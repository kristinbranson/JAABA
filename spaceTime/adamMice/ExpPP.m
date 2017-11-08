classdef ExpPP
% Experiment postprocessor
  
properties (Constant)
  PARAMETER_DESCRIPTION_FILE = 'ExpPP.params.txt';
  PARAMETER_DESCRIPTIONS = ExpPP.loadParameterDescriptions();
      
  ALLBEHS = {'Atmouth' 'Badgrab' 'Botharm' 'Chew' 'Grab' 'Handopen' 'Lefthand' 'Lift' 'Reach' 'Sup'};
  BASICBEHAVIORS = {'Lift' 'Handopen' 'Grab' 'Sup' 'Atmouth' 'Chew'};
  NBASICBEHAVIORS = numel(ExpPP.BASICBEHAVIORS);
  
  DEFAULT_MINBOUTSIZE = cell2struct(repmat({3},1,ExpPP.NBASICBEHAVIORS),ExpPP.BASICBEHAVIORS,2);
  
  FRMRATE = 180; % frames/sec. TODO this is not constant?
  TRXFNAME = 'trx.mat';
  
  FLDS_EXPORT_MD = {'id' 'upper' 'mouse' 'date' 'exp'};  
  
  FIELDNAME_CUSTOMGROUPLIST = 'customgrpNames';
  CUSTOMGROUPNAME_PAT = 'customgrp_%s';
  CUSTOMGROUP_NOTPREFIX = 'not_';
  
  TRIALSUCCESSCOLOR = struct(...
    'successful_1grab',[0 0.55 0], ...
    'successful_2plusgrab',[0 0 1], ...
    'unsuccessful',[0 0 0]);
end

methods (Static) % data/metadata load/export
  
  function data = loadexps(dirs,jabOrScoreFiles,varargin)
    % data = loadexps(dirs,jabOrScoreFiles,varargin)
    %
    % dirs: char or cellstr.
    % jabOrScoreFiles: cellstr of filenames. This should contain either a)
    % full paths to jabfiles, or b) scorefile names (eg
    % scores_Liftm119.mat), for the mouse+behaviors of interest. At the
    % moment, this set must include all "basic" behaviors (Lift, Grab,
    % Atmouth, Chew). If empty, you will be prompted to select Jabfiles.
    % s: data structure, one element per trial/experiment. 
    %
    % optional PVs:
    % - dorecurse (scalar logical). If true, treat dirs as root dirs.
    % - mouse (positive integer). If not provided, will try to infer mouse
    % from jabfiles.
    % - condition (char). Mouse condition/manipulation.
    % - includelabels (scalar logical). If true, manual labels from jabfiles
    % will be loaded along with score data.
    %
    % Note: dorecurse==false does not guarantee that
    % numel(data)==numel(dirs). data may have fewer elements than dirs, if
    % one or more experiments in dirs is missing scorefiles or is otherwise
    % invalid.
    
    if ischar(dirs)
      dirs = cellstr(dirs);
    end   
    assert(iscellstr(dirs),'Expected ''dirs'' to be a char or cellstr.');
    
    % uiGetJabFiles if necessary
    if isempty(jabOrScoreFiles)
      tmpstr = civilizedStringFromCellArrayOfStrings(ExpPP.BASICBEHAVIORS);
      tmpstr = sprintf('Select jab files for: %s',tmpstr);
      [tfsuccess,jabOrScoreFiles] = ExpPP.uiGetJabFiles('promptstr',tmpstr);
      if ~tfsuccess
        data = [];
        return;
      end
    end
    
    % Get:
    % * scorefiles: Nx1 cellstr of score files
    % * behnames: Nx1 cellstr of behaviornames for each scorefile
    tfJabs = ~isempty(regexp(jabOrScoreFiles{1},'.jab$','once'));
    if tfJabs
      [scorefiles,behnames] = cellfun(@Macguffin.jabfileScoresBehs,jabOrScoreFiles,'uni',0);      
      scorefiles = cellfun(@(x)x(:),scorefiles,'uni',0);
      behnames = cellfun(@(x)x(:),behnames,'uni',0);
      scorefiles = cat(1,scorefiles{:});
      behnames = cat(1,behnames{:});
    else
      scorefiles = jabOrScoreFiles;
      behnames = cellfun(@(x){x},scorefiles,'uni',0); % best guess for behaviors: just use scorefilename
    end   
    assert(numel(scorefiles)==numel(behnames));
    % Now:
    % Convert behnames to abstract behnames, as possible
    for iScoreFile = 1:numel(scorefiles)
      bnames = behnames{iScoreFile};
      tf = regexpmatch(bnames,ExpPP.ALLBEHS,'caseinsens',true);
      if nnz(tf)==0
        warning('ExpPP:loadexps','Unrecognized score/behavior ''%s''.',bnames);
      elseif nnz(tf)==1
        behnames{iScoreFile} = ExpPP.ALLBEHS{tf};
      else
        warning('ExpPP:loadexps','Score/behavior ''%s'' matches multiple abstract behaviors.',bnames);
      end
    end
    
    [dorecurse,mouse,condition,includelabels] = myparse(varargin,...
      'dorecurse',false,...
      'mouse',[],...
      'condition',[],...
      'includelabels',false);
    validateattributes(dorecurse,{'logical' 'numeric'},{'scalar'},'','dorecurse');
    if ~isempty(mouse)
      validateattributes(mouse,{'numeric'},{'positive' 'integer'},'','mouse');
    end
    validateattributes(includelabels,{'numeric' 'logical'},{'scalar'});
    
    % get the expdirs
    isExpDirFcn = @(d)ExpPP.isExpDir(d,scorefiles);
    if dorecurse
      founddirs = cell(0,1);
      foundtf = false(0,1);
      for i = 1:numel(dirs)
        tmpmap = depthFirstSearch(dirs{i},isExpDirFcn);
        keys = tmpmap.keys;
        vals = cell2mat(tmpmap.values);
        founddirs = [founddirs;keys(:)]; %#ok<AGROW>
        foundtf = [foundtf;vals(:)]; %#ok<AGROW>
      end
    else
      founddirs = dirs;
      foundtf = cellfun(isExpDirFcn,dirs);
      nLoad = nnz(foundtf);
      nFound = numel(founddirs);
      fprintf(1,'Loading %d/%d experiments specified.\n',nLoad,nFound);
    end
    assert(isequal(size(founddirs),size(foundtf)));
    expdirs = founddirs(foundtf);
    if isempty(expdirs)
      warningNoTrace('ExpPP:noDirs','No loadable experiments found.');
    end
        
    % load exps
    if tfJabs
      jabs = jabOrScoreFiles;
    else
      jabs = [];
    end
    if includelabels
      assert(tfJabs && isscalar(jabOrScoreFiles),'''includelabels'' option only supported for scalar jabfile input.');
      Q = loadAnonymous(jabOrScoreFiles{1});
      
      nExp = numel(expdirs);
      data = cell(nExp,1);
      for i = 1:nExp
        tf = strcmp(expdirs{i},Q.expDirNames);
        switch nnz(tf)
          case 0,    labels = Labels.labels(1); % expdir not in jab; use empty labels structure
          case 1,    labels = Q.labels(tf);
          otherwise, assert(false,'Repeated expdirs in jabfile ''%s''.',jabOrScoreFiles{1});
        end
        data{i} = ExpPP.loadsingleexp(expdirs{i},scorefiles,behnames,...
          'condition',condition,'jabs',jabs,'labels',labels);
      end
      data = cell2mat(data);
    else     
      data = cellfun(@(d)ExpPP.loadsingleexp(d,scorefiles,behnames,...
        'condition',condition,'jabs',jabs),expdirs);
    end
    data = data(:);
    
    % check mouse
    % TODO: guess this is pointless unless actual filtering is done?
    if ~isempty(mouse)
      for i = 1:numel(data)
        if ~isequal(data(i).mouse,mouse)
          warning('ExpPP:mouseMismatch',...
            'Experiment %s: experiment may have been run on wrong mouse.',data(i).id);
        end
      end
    end
  end
  
  function [tfsuccess,jabfiles] = uiGetJabFiles(varargin)
    % tfsuccess: if false, user canceled GUI
    
    promptstr = myparse(varargin,...
      'promptstr','Select jab files');
    
    jabpath = ExpPP.loadConfigVal('jabpath');
    if isempty(jabpath)      
      jabpath = pwd;
    end
    jabfiles = uipickfiles('Prompt',promptstr,...
      'FilterSpec',jabpath,'Type',{'*.jab','JAB-files'},'SmartAdd',true);
    if ~iscell(jabfiles) || isempty(jabfiles), % canceled
      tfsuccess = false;
      jabfiles = [];
    else
      tfsuccess = true;
      jabpath = fileparts(jabfiles{1});
      ExpPP.saveConfigVal('jabpath',jabpath);      
    end
  end
  
  function tf = isExpDir(d,scorefnames)
    % - a dir is an expdir if it contains all desired scorefiles
    % - if a dir contains some but not all scorefiles, a warning is thrown
    % - in the future, can try to handle such "partial" dirs
    
    tfscr = cellfun(@(x)exist(fullfile(d,x),'file'),scorefnames);    
    if all(tfscr)
      tf = true;
    elseif any(tfscr)
      warningNoTrace('ExpPP:partialExpDir',...
        '%s score files missing: %s',d,...
        civilizedStringFromCellArrayOfStrings(scorefnames(~tfscr)));
      tf = false;
    else
      tf = false;
    end    
  end
  
  function s = loadsingleexp(exp,scorefnames,abstractbehnames,varargin)
    % scorefnames: cellstr. short scorefilenames, eg 'scores_Grabm76'
    % abstractbehnames: cellstr of abstract behavior names for scorefnames, 
    %   eg {'Grab' 'Lift'}. behnames{i} labels the scores in
    %   scorefnames{i}. behnames{i} will typically differ from the
    %   .behaviorName field in the scorefnames{i} scorefile.
    % 
    % Optional PVs:
    %   - condition
    %   - jabs
    %   - labels. scalar labels structure for this experiment. If provided,
    %   labeled bouts are added to s, for behaviors matching the specific
    %   (not abstract) behaviornames specified in scorefnames{i}. Note,
    %   currently only a scalar 'labels' value is supported, ie the labels
    %   are expected to come from a multi-classifier jab that encompasses
    %   the behaviors in scorefnames.
    %
    % s: scalar struct containing all loaded metadata/data for experiment
    % exp.
    %
    % Currently, exp must contain i) a valid trx file, and ii) valid
    % scorefiles for all elements of scorefnames. 
    
    % AL 20150203. The restriction that exp must contain valid scorefiles
    % for all scorefnames may seem overly restrictive. Note however that 
    % if the optional 'labels' argument is provided, all scorefiles must
    % exist in order to specify the concrete behaviornames to look for in
    % labels.
    
    assert(exist(exp,'dir')>0,'Dir not found: ''%s''.',exp);
    if ischar(scorefnames)
      scorefnames = {scorefnames};
    end
    assert(iscellstr(scorefnames) && iscellstr(abstractbehnames) && numel(scorefnames)==numel(abstractbehnames),...
      'Invalid input arguments scorefnames/behnames.');
    assert(numel(abstractbehnames)==numel(unique(abstractbehnames)),...
      'Input/specified behavior names must all be distinct.');
        
    [condition,jabs,labels] = myparse(varargin,...
      'condition',[],...
      'jabs',[],...
      'labels',[]);
    
    % Start with metadata
    sMD = ExpPP.parseexpname(exp);
    sMD.condition = condition;
    sMD.jabs = jabs;
        
    trxfn = fullfile(exp,ExpPP.TRXFNAME);
    trx = load(trxfn,'trx');
    trx = trx.trx;
    sMD.trxt0 = min([trx.firstframe]);
    sMD.trxt1 = max([trx.endframe]);

    % Load scores
    scoredata = cellfun(@(x)ExpPP.loadandcheckscores(fullfile(exp,x)),scorefnames,'uni',0);
    Nscorefiles = numel(scorefnames);
    
    % - Check scoredata behaviornames against abstract names  
    for i = 1:Nscorefiles
      sdata = scoredata{i};
      if isfield(sdata,'behaviorName')
        % new-style scorefiles store behaviorname; check against abstract names
        abname = abstractbehnames{i};
        bname = sdata.behaviorName;
        if isempty(regexpi(bname,abname,'once'))
          warningNoTrace('ExpPP:scoreBehMismatch',...
            'Possible behavior mismatch. Concrete behavior: ''%s''. Abstract behavior: ''%s''.',...
            bname,abname);
        end
      end
    end
    
    tfLabels = ~isempty(labels);
    if tfLabels
      assert(isscalar(labels) && isstruct(labels),'Only scalar ''labels'' argument supported.');
      Labels.verifyLabels(labels);
      
      if isempty(labels.flies)
        % empty labels structure, no labels
        lblt0s = nan(1,0);
        lblt1s = nan(1,0);
        lblnames = cell(1,0);
      else     
        assert(labels.flies==1);      
        assert(labels.off==0);
      
        lblt0s = labels.t0s{1};
        lblt1s = labels.t1s{1};
        lblnames = labels.names{1};
        
        lbltts = labels.timelinetimestamp{1};
        lblttsfns = fieldnames(lbltts);
        assert(numel(lblttsfns)>1,'Currently expect multiclassifier labels.');    
      end
      
      nLblBouts = numel(lblt0s);
    end
    
    assert(isequal(Nscorefiles,numel(scoredata),numel(abstractbehnames)));
    for i = 1:Nscorefiles
      % scores (predictions, machine-generated)
      scr = struct();
      scr.tfscorefileexists = true;
      scr.scorefilename = scorefnames{i};
      scr = structmerge(scr,scoredata{i});
      
      % trx check
      t1 = min(sMD.trxt1,numel(scr.postprocessed));
      if t1 < sMD.trxt1
        fprintf('Number of frames in scores file %s = %d < number of frames in trajectory file %s = %d\n',...
          scr.scorefilename,t1,trxfn,sMD.trxt1);
      end
            
      scr = structprefixfieldnames(scr,[abstractbehnames{i} '_']);
      sMD = structmerge(sMD,scr);
      
      % manual labels
      if tfLabels
        scrdata = scoredata{i};
        beh = scrdata.behaviorName;
        nobeh = Labels.noBehaviorName(beh);
        if exist('lblttsfns','var')>0 && ~any(strcmp(beh,lblttsfns))
          warning('ExpPP:behMissing',...
            'While loading ''%s'', behavior ''%s'' not found in labels.timelinetimestamp.',...
            exp,beh);
          % Should be impossible now for beh/nobeh to show up in lblnames below
        end

        assert(scrdata.tStart==1);
        labl = struct();
        % label: 1/0/-1 for beh/unk/nobeh
        labl.label = zeros(size(scrdata.postprocessed));
        labl.t0sPos = nan(1,0);
        labl.t1sPos = nan(1,0);
        labl.t0sNeg = nan(1,0);
        labl.t1sNeg = nan(1,0);
        for iBout = 1:nLblBouts
          switch lblnames{iBout}
            case beh
              idx = lblt0s(iBout):lblt1s(iBout)-1;
              labl.label(idx) = 1;
              labl.t0sPos(1,end+1) = lblt0s(iBout);
              labl.t1sPos(1,end+1) = lblt1s(iBout);
            case nobeh
              idx = lblt0s(iBout):lblt1s(iBout)-1;
              labl.label(idx) = -1;
              labl.t0sNeg(1,end+1) = lblt0s(iBout);
              labl.t1sNeg(1,end+1) = lblt1s(iBout);
          end
        end
        
        labl = structprefixfieldnames(labl,[abstractbehnames{i} '_labl_']);
        sMD = structmerge(sMD,labl);
      end
    end
    
    s = sMD;
  end
  
  function s = loadandcheckscores(scrfile)
    % s = loadandcheckscores(scrfile)
    % s: struct array, score data in scrfile (.allScores field)
    assert(exist(scrfile,'file')==2,'File not found: ''%s''.',scrfile);
    
    s = load(scrfile);
    if iscell(s.allScores)
      assert(isscalar(s.allScores));
      s.allScores = s.allScores{1};
    end
    assert(isstruct(s.allScores));    
    s = structmerge(s,s.allScores);
    s = rmfield(s,'allScores');
        
    for i = 1:numel(s)
      FLDS = {'scores' 't0s' 't1s' 'postprocessed'};
      for f = FLDS,f=f{1}; %#ok<FXSET>
        assert(isequal(s(i).(f){:}));
        s(i).(f) = s(i).(f){1};
      end
      FLDS = {'tStart' 'tEnd'};
      for f = FLDS,f=f{1}; %#ok<FXSET>
        assert(isscalar(unique(s(i).(f))));
        s(i).(f) = s(i).(f)(1);
      end
    end
    
    % probably too strict
    %assert(s.tStart==1 && s.tEnd==ExpPP.NUMSAMPS && numel(s.postprocessed)==ExpPP.NUMSAMPS);
  end
  
  function sMD = parseexpname(expname)
    % parse expname for metadata.
    
    sMD = struct();    
    sMD.expfull = expname;
    
    % upper/date/exp/file/trial
    [set,sMD.exp] = myfileparts(expname);
    TRIALPAT = '_v(?<trl>[0-9]{2,3})';
    tmp = regexp(sMD.exp,TRIALPAT,'names');
    if ~isempty(tmp)
      sMD.trial = str2double(tmp.trl);
    else
      sMD.trial = nan;
    end 
    
    [sMD.upper,sMD.date] = myfileparts(set);
    if strcmp(sMD.date,'Good')
      % special case; extra intervening dir
      [sMD.upper,sMD.date] = myfileparts(sMD.upper);
    end
    DATEPAT1 = '[0-9]{1,2}_[0-9]{1,2}_[0-9]{1,2}';
    if ~isempty(regexp(sMD.date,DATEPAT1,'once'))
      sMD.date = datestr(sMD.date,'yyyymmdd');
    end
    MOUSEPAT1 = '_[Mm]ouse([0-9]+)';
    MOUSEPAT2 = '[mM]([0-9]+)[/_]';
    %MOUSEPAT2 = '_[mM]([0-9]+)/';
    tmp1 = regexp(set,MOUSEPAT1,'tokens');
    tmp2 = regexp(set,MOUSEPAT2,'tokens');
    if ~isempty(tmp1)
      sMD.mouse = tmp1{1}{1};
    elseif ~isempty(tmp2)
      sMD.mouse = tmp2{1}{1};
    else
      warningNoTrace('No mouse parsed from %s.',set); %#ok<WNTAG>
      sMD.mouse = '<unknown>';
    end
    
    %sMD.id = sprintf('%s#%s#%s',sMD.mouse,sMD.date,sMD.exp);
    IDPATS = {'/tier2/hantman' 'Y:' ''};
    for pat = IDPATS, pat=pat{1}; %#ok<FXSET>
      if strncmp(expname,pat,numel(pat))
        sMD.id = expname(numel(pat)+1:end);
        break;
      end
    end
    sMD.id = regexprep(sMD.id,'\\','/');
  end
  
  function export(data,fname,varargin)
    % export(data,fname,p1,v1,...)
    % 
    % data: data struct.
    % fname: file to export into. Leave empty to be prompted.
    % Optional PVs: 
    %   - 'format'. One of {'csv' 'xls'}, defaults to 'csv'.
    %   - 'expppstatprefix', defaults to 'auto'
    % 
    % Note: we export some additional fields not present in the data 
    % structure purely for convenience and eg interfacing with Uri.
    
    [format,expppstatprefix] = myparse(varargin,...
      'format','csv',...
      'expppstatprefix','auto');
    assert(ismember(format,{'csv' 'xls'}),'Invalid ''format'' specification.');
    if ~ispc && strcmp(format,'xls')
      warning('ExpPP:xlsNonPc','Export to xls supported only on Windows. Using csv format.');
      format = 'csv';
    end
    
    % confirm/verify filename
    if ~exist('fname','var') || isempty(fname)
      exportpath = ExpPP.loadConfigVal('exportpath');
      if isempty(exportpath)
        exportpath = pwd;
      end      
      switch format
        case 'csv', ext = 'csv';
        case 'xls', ext = 'xlsx';
      end
      filterspec = fullfile(exportpath,['*.' ext]);
      
      [f,p] = uiputfile(filterspec,'Export data');
      if isequal(f,0) || isequal(p,0)
        return;
      end
      fname = fullfile(p,f);
    else
      if exist(fname,'file')
        warning('ExpPP:overwrite','File ''%s'' exists and will be overwritten.',fname);
      end
    end
    
    %%% Compute additional 'convenience' fields
    flds = fieldnames(data);
    baseFlds = cell(0,1); % create convenience fields for these base fields. Base fields should take frames indices for values.

    % Determine behaviors present
    pat = sprintf('^%s_([a-zA-Z]+)_0$',expppstatprefix);
    tmptoks = regexp(flds,pat,'tokens');
    tmptoks = cat(2,tmptoks{:});
    behs = cat(1,tmptoks{:});
    
    % Fields like auto_Lift_0, auto_Lift_0_bl
    for b = behs(:)', b = b{1}; %#ok<FXSET>
      pat = sprintf('^%s_%s_[01]',expppstatprefix,b);
      baseFlds = [baseFlds; flds(regexpmatch(flds,pat))]; %#ok<AGROW>
    end
    
    % special case
    baseFlds{end+1,1} = sprintf('%s_del_Lift1_Lift0',expppstatprefix);
    
    % GS* fields
    baseFlds = [baseFlds; flds(regexpmatch(flds,sprintf('^%s_.*GS0[01]_',expppstatprefix)))];
    %baseFlds = [baseFlds; flds(regexpmatch(flds,sprintf('%s_del_GS0[01]_',expppstatprefix)))];
    baseFlds = [baseFlds; flds(regexpmatch(flds,sprintf('^%s_.*GSSS_',expppstatprefix)))];
    %baseFlds = [baseFlds; flds(regexpmatch(flds,sprintf('%s_del_GSSS_',expppstatprefix)))];
    
    fldGrabSuccessType = sprintf('%s_Grab_successtype',expppstatprefix);
    Nexp = numel(data);
    for i = 1:Nexp
      trialtype = data(i).(fldGrabSuccessType);
      for fldfrm = baseFlds(:)', fldfrm = fldfrm{1}; %#ok<FXSET>
        frm = data(i).(fldfrm);
        fldcol1 = sprintf('%s_successfultrial',fldfrm); 
        fldcol2 = sprintf('%s_unsuccessfultrial',fldfrm);
        fldcol3 = sprintf('%s_successfultrial_1grab',fldfrm); 
        fldcol4 = sprintf('%s_successfultrial_2plusgrab',fldfrm);
        switch trialtype
          case 'successful_1grab'
            data(i).(fldcol1) = frm;
            data(i).(fldcol2) = [];
            data(i).(fldcol3) = frm;
            data(i).(fldcol4) = [];
          case 'successful_2plusgrab'
            data(i).(fldcol1) = frm;
            data(i).(fldcol2) = [];
            data(i).(fldcol3) = [];
            data(i).(fldcol4) = frm;
          case 'unsuccessful'
            data(i).(fldcol1) = [];
            data(i).(fldcol2) = frm;
            data(i).(fldcol3) = [];
            data(i).(fldcol4) = [];
          otherwise
            assert(false,'Unknown value of ''%s''.',fldGrabSuccessType);
        end
      end
    end
    
    % form reduced data structure
    flds = fieldnames(data);

    fldsMD = ExpPP.FLDS_EXPORT_MD(:);
    tfparam = regexpmatch(flds,'^param_');
    tfspecial = ismember(flds,{sprintf('%s_Grab_success',expppstatprefix);...
                               sprintf('%s_Grab_successtype',expppstatprefix)}) | ...
                regexpmatch(flds,'^customgrp_');
    tfbasic = regexpmatch(flds,sprintf('^%s_[a-zA-Z]+_[01]$',expppstatprefix)) | ...
              regexpmatch(flds,sprintf('^%s_[a-zA-Z]+_[01]_bl$',expppstatprefix)) | ...
              regexpmatch(flds,sprintf('^%s_[a-zA-Z]+_num$',expppstatprefix)) | ...
              regexpmatch(flds,sprintf('^%s_laser',expppstatprefix));
    tfremain = regexpmatch(flds,sprintf('^%s_GS0[01]_',expppstatprefix)) | ...      
               regexpmatch(flds,sprintf('^%s_GSSS_',expppstatprefix)) | ...
               regexpmatch(flds,sprintf('^%s_del_',expppstatprefix)) | ...
               regexpmatch(flds,'successfultrial');
            
    flds2keep = [fldsMD;flds(tfparam);flds(tfspecial);flds(tfbasic);flds(tfremain)]; % order intentional    
    % for now, exclude all labl_* fields
%     tfLabl = regexpmatch(flds2keep,'^labl_');
%     flds2keep(tfLabl,:) = [];
    
    sout = rmfield(data,setdiff(flds,flds2keep));
    sout = orderfields(sout,flds2keep);

    structexport(sout,fname,'format',format);
    
    ExpPP.saveConfigVal('exportpath',fileparts(fname));
  end
  
  function exportStats(stats,fname)
    % confirm/verify filename
    if ~exist('fname','var')
      exportpath = ExpPP.loadConfigVal('exportpath');
      if isempty(exportpath)
        exportpath = pwd;
      end
      filterspec = fullfile(exportpath,'*.csv');
      
      [f,p] = uiputfile(filterspec,'Export statistics');
      if isequal(f,0) || isequal(p,0)
        return;
      end
      fname = fullfile(p,f);
    else
      if exist(fname,'file')
        warning('ExpPP:overwrite','File ''%s'' exists and will be overwritten.',fname);
      end
    end
    
    fh = fopen(fname,'w');
    fprintf(fh,'field,set,N,mean,median,std\n');
    flds = fieldnames(stats);
    for f = flds(:)', f = f{1}; %#ok<FXSET>
      %fprintf(fh,'%s\n',f);
      sets = fieldnames(stats.(f));
      for s = sets(:)', s = s{1}; %#ok<FXSET>
        el = stats.(f).(s);
        fprintf(fh,'%s,%s,%d,%.3f,%.3f,%.3f\n',f,s,el.N,el.mean,el.median,el.std);
      end
    end
    
    fclose(fh);     
  end
  
  function data = mergeData(varargin)
    % data = ExpPP.mergeData(d1,d2,d3,...)
    
    dMerge = varargin;
    dMerge = dMerge(~cellfun(@isempty,dMerge));
    if isempty(dMerge)
      data = [];
      return;
    end
    
    nMerge = numel(dMerge);
    cgrpsMerge = cell(nMerge,1);
    for i = 1:nMerge
      [dMerge{i},cgrpsMerge{i}] = ExpPP.ensureCustomGroupInit(dMerge{i});
    end
    cgrpsAll = cat(2,cgrpsMerge{:});
    cgrpsAll = unique(cgrpsAll);
    
    for i = 1:nMerge
      dm = dMerge{i};
      
      % Add missing custom groups 
      cgrpsAdd = setdiff(cgrpsAll,cgrpsMerge{i});
      for iGrp = 1:numel(cgrpsAdd)
        dm = ExpPP.createCustomGroup(dm,false(size(dm)),cgrpsAdd{iGrp});
      end
      
      if i==1
        data = dm;
      else
        % check experiments seen
        %expfullSeen = {data.expfull}';
        %expSeen = {data.exp}';
        idSeen = {data.id}';
        assert(isequal(numel(data),numel(unique(idSeen))));
        
        %expfull = {dm.expfull}';
        id = {dm.id}';
        tfRpt = ismember(id,idSeen);
        if any(tfRpt)
          idIgnore = id(tfRpt);
          warning('ExpPP:mergeData','Ignoring repeated experiments found in input %d: %s\n',...
            i,civilizedStringFromCellArrayOfStrings(idIgnore));
        end
        dm = dm(~tfRpt);
        
        dm = orderfields(dm,data);
        data = [data(:); dm(:)];
      end
    end
  end  
    
end

methods (Static) % custom groups
  
  function [data,currentCustomGroups] = ensureCustomGroupInit(data)
    % Ensure custom-group state is initialized on data; return current
    % custom group list
    assert(~isempty(data),'Currently not supported for empty data.');
    if ~isfield(data,ExpPP.FIELDNAME_CUSTOMGROUPLIST)
      [data.(ExpPP.FIELDNAME_CUSTOMGROUPLIST)] = deal(cell(1,0)); % row vec
    end
    currentCustomGroups = data(1).(ExpPP.FIELDNAME_CUSTOMGROUPLIST);
  end
  
  function data = createCustomGroup(data,tf,grpname)
    validateattributes(tf,{'logical'},{'vector','numel' numel(data)});
    assert(ischar(grpname));
        
    [data,cgrps] = ExpPP.ensureCustomGroupInit(data);    
    if ismember(grpname,cgrps)
      error('ExpPP:existingGroup','Custom group with name ''%s'' already exists.',grpname);
    end
    
    fldtf = sprintf(ExpPP.CUSTOMGROUPNAME_PAT,grpname);
    if ~isvarname(fldtf)
      error('ExpPP:invalidFieldname',...
        'Custom group name ''%s'' leads to invalid structure fieldname. Please use a simpler group name.',...
        grpname);
    end
    cgrps{1,end+1} = grpname;
    [data.(ExpPP.FIELDNAME_CUSTOMGROUPLIST)] = deal(cgrps);
    tf = num2cell(tf);
    [data.(fldtf)] = deal(tf{:});
  end
  
  function data = removeCustomGroup(data,grpname)
    [data,cgrps] = ExpPP.ensureCustomGroupInit(data);
    if ~ismember(grpname,cgrps)
      error('ExpPP:nosuchGroup',...
        'Custom group ''%s'' does not exist in data.',grpname);
    end
    fldtf = sprintf(ExpPP.CUSTOMGROUPNAME_PAT,grpname);
    assert(isfield(data,fldtf))
    
    cgrps = setdiff(cgrps,grpname);
    cgrps = cgrps(:)'; % ensure cgrps is a row vec
    [data.(ExpPP.FIELDNAME_CUSTOMGROUPLIST)] = deal(cgrps);
    data = rmfield(data,fldtf);
  end
  
  function tf = findCustomGroup(data,grpname)
    assert(ischar(grpname));
    assert(~isempty(data) && any(strcmp(grpname,data(1).(ExpPP.FIELDNAME_CUSTOMGROUPLIST))),...
      'Group ''%s'' does not exist in data.',grpname);
    
    fldtf = sprintf(ExpPP.CUSTOMGROUPNAME_PAT,grpname);
    tf = [data.(fldtf)];
  end
  
  function data = addNotGroups(data)
    % For not particularly good reasons, a client wants explicit not-groups
    % for each custom group

    [data,cgrps] = ExpPP.ensureCustomGroupInit(data);
    Ngrp = numel(cgrps);
    for i = 1:Ngrp
      cgrpNot = [ExpPP.CUSTOMGROUP_NOTPREFIX cgrps{i}];
      if ismember(cgrpNot,cgrps)
        warning('ExpPP:groups',...
          'Complement group ''%s'' already exists. Not creating.',cgrpNot);
        continue;
      end
      tf = ExpPP.findCustomGroup(data,cgrps{i});
      data = ExpPP.createCustomGroup(data,~tf,cgrpNot);
    end
  end
  
  function data = cleanNotGroups(data)
    PFIX = ExpPP.CUSTOMGROUP_NOTPREFIX;
    
    [data,cgrps] = ExpPP.ensureCustomGroupInit(data);
    Ngrp = numel(cgrps);
    for i = Ngrp:-1:1
      if strncmp(cgrps{i},PFIX,numel(PFIX))
        notgrp = cgrps{i};
        grp = notgrp(numel(PFIX)+1:end);
        tf = strcmp(grp,cgrps);
        assert(nnz(tf)<=1);
        if any(tf)
          tfG = ExpPP.findCustomGroup(data,grp);
          tfNotG = ExpPP.findCustomGroup(data,notgrp);
          if all(xor(tfG,tfNotG))
            % precise complements
            data = ExpPP.removeCustomGroup(data,notgrp);
          end
        end
      end
    end 
  end
  
  function data = convertCustomGroupsToStrs(data)
    % Convert custom group vals from 0/1 to <grpname>/not<grpname> 
    
    [data,cgrps] = ExpPP.ensureCustomGroupInit(data);
    Ngrp = numel(cgrps);
    for i = 1:Ngrp
      g = cgrps{i};
      tf = ExpPP.findCustomGroup(data,g);
      newg = cell(size(tf));
      newg(tf) = {g};
      newg(~tf) = {['~' g]};
      
      gfld = sprintf(ExpPP.CUSTOMGROUPNAME_PAT,g);
      [data.(gfld)] = deal(newg{:});
    end
  end
  
  function data = convertLogicalFieldToStrs(data,fld01,str0,str1)
    % Convert a field with logical dichotomous data to string values:
    % 0 -> str0, eg 'notCNO'
    % 1 -> str1, eg 'CNO'
    
    for i = 1:numel(data)
      v = data(i).(fld01);
      assert(islogical(v));
      if v
        data(i).(fld01) = str1;
      else
        data(i).(fld01) = str0;
      end
    end
  end
  
  function [data,gPair] = addPairwiseGroups(data,gFlds)
    % Form/add all pairwise grouping variables.
    %
    % gFlds: vector cellstr, list of grouping vars. data must take on 
    % categorical data over these fields.
    % gPair: vector cellstr, list of pairwise-grouping vars added. The
    % values for these grouping vars are char categoricals.
    
    gFlds = gFlds(:);
    
    % check gFlds
    tf = false(size(gFlds));
    valsCellStr = cell(size(gFlds));
    for i = 1:numel(gFlds)
      gF = gFlds{i};
      [tf(i),msg,vals] = ExpPP.checkCategoricalVar(data,gF);
      if tf(i)
        valsCellStr{i} = ExpPP.ensureStringCategorical(vals,gF);
      else
        warning('ExpPP:invalidCatVar',...
          'Skipping invalid categorical variable ''%s'': %s',gF,msg);
      end
    end
    gFlds = gFlds(tf);
    valsCellStr = valsCellStr(tf);   
        
    % Add pairwise
    nGrp = numel(gFlds);
    gPair = cell(0,1);
    for i1 = 1:nGrp
      g1 = gFlds{i1};
      v1 = valsCellStr{i1};      
      for i2 = i1+1:nGrp
        g2 = gFlds{i2};
        v2 = valsCellStr{i2};
      
        gComb = [g1 '__' g2];
        if isfield(data,gComb) % unlikely
          warning('ExpPP:groupCollision',...
            'Field ''%s'' for combined-group already exists in data. Not adding combined-custom group.',...
            gComb);
          continue;
        end
        
        assert(isequal(size(v1),size(v2),size(data)));
        for iData = 1:numel(data)
          data(iData).(gComb) = sprintf('%s & %s',v1{iData},v2{iData});
        end
        
        gPair{end+1,1} = gComb; %#ok<AGROW>
      end
    end
  end
  
  function [tf,msg,g] = checkCategoricalVar(data,gFld)
    % Checks that a field of a data struct is a char, double, or logical
    % categorical.
    %
    % tf: if true, gFld represents an acceptable categorical field over data.
    % msg: message detail, used only if tf==false
    % g: either a cellstr, logical, or numeric vector containing
    %   categorical values. Same size as data.
    
    tf = false;
    msg = '';
    g = {data.(gFld)}';
    
    if any(cellfun(@isempty,g))
      msg = 'One or more empty values';
      return;
    end
    
    cls = cellfun(@class,g,'uni',0);
    cls = unique(cls);
    if ~isscalar(cls)
      msg = 'Values with differing types';
      return;
    end
    cls = cls{1};
    
    switch cls
      case 'char'
        % none
      case {'logical' 'double'}
        g = cell2mat(g);
      otherwise
        msg = sprintf('Unknown class: %s',cls);
        return;
    end
    
    assert(isequal(size(g),size(data)));
    tf = true;
  end
  
  function gstr = ensureStringCategorical(g,gName)
    % Convert a categorical vector to cellstr.
    % 
    % g: vector of categorical values: either double, logical, or cellstr
    % gName: fieldname, used when g is logical of "effectively logical"
    %
    % gstr: cellstr
    
    if iscellstr(g)
      gstr = g;
    elseif islogical(g) || isnumeric(g) && all(ismember(g(:),[0 1]))
      tf = logical(g);
      gstr = cell(size(g));
      gstr(tf) = {gName};
      gstr(~tf) = {['~' gName]};
    elseif isnumeric(g)
      gstr = arrayfun(@num2str,g,'uni',0);
    else
      assert(false,'Unknown type for categorical variable ''%s''.',gName);
    end
  end
  
end

methods (Static) % parameters 
  
  function tbl = loadParameterDescriptions
    % Manually read parameter file (want to use readtable here but this 
    % requires R2013b or later)
    fh = fopen(ExpPP.PARAMETER_DESCRIPTION_FILE);
    if isequal(fh,-1)
      error('ExpPP:errRead','Cannot read postprocessing parameter file.');
    end
    out = textscan(fh,'%s%s%s%s%s','delimiter',',','collectoutput',true);
    fclose(fh);
    out = out{1};
    assert(isequal(out(1,:),{'name' 'type' 'nanok' 'default' 'desc'}));
    out = out(2:end,:);
    tbl = struct();
    tbl.name = out(:,1);
    tbl.type = out(:,2);
    tbl.nanok = str2double(out(:,3));
    tbl.default = str2double(out(:,4));
    tbl.desc = out(:,5);
  end
  
  function s = formParamStruct(names,values)
    % eg: names = {'fillsize' 'minboutsize#Lift'}; vals = {1 2};
    %     => s.fillsize = 1; s.minboutsize.Lift = 2;
    
    assert(iscellstr(names) && iscell(values));
    assert(isequal(numel(names),numel(values)));
    Nflds = numel(names);
    
    paramDescs = ExpPP.PARAMETER_DESCRIPTIONS;

    s = struct();
    for i = 1:Nflds
      prm = names{i};
      val = values{i};
      
      iParamDesc = strcmp(prm,paramDescs.name);
      switch nnz(iParamDesc)
        case 0
          warningNoTrace('ExpPP:formParamStruct',...
            'Unrecognized parameter name ''%s'' will be ignored.',prm);
          continue;
        case 1
          % none
        otherwise
          assert(false);
      end
      
      % check val
      if isnan(val)
        if paramDescs.nanok(iParamDesc)
          % none
        else
          warningNoTrace('ExpPP:nanVal',...
            'Parameter ''%s'' cannot be NaN. Using factory-default value of %d.\n',...
            prm,paramDescs.default(iParamDesc));
          val = paramDescs.default(iParamDesc);
        end
      else
        switch paramDescs.type{iParamDesc}
          case 'nat'
            if floor(val)~=val
              warningNoTrace('ExpPP:fracVal',...
                'Parameter ''%s'' should be an integer. Rounding to %d.\n',prm,round(val));
              val = round(val);
            end
            if val<0
              warningNoTrace('ExpPP:negVal',...
                'Parameter ''%s'' cannot be negative. Using 0.\n',prm);
              val = 0;
            end
          otherwise
            assert(false,'All types currently nat');
        end
      end
      
      % assign
      subparams = regexp(prm,'#','split');
      switch numel(subparams)
        case 1
          s.(subparams{1}) = val;
        case 2
          s.(subparams{1}).(subparams{2}) = val;
        otherwise
          assert(false,'Unexpected parameter.');
      end
    end

  end
  
  function expparammap = loadExperimentSpecificParams(fname)
    % expparammap = loadExperimentSpecificParams(fname)
    %
    % fname: optional. filename of csv containing
    % experiment-specific params. This file must have a header row of
    % fieldnames. The field 'exp' is mandatory; this field uniquely
    % identifies experiments by name. All other fields should be a subset
    % of parameters with names as listed in ExpPP.params.txt.
    %
    % expparammap: containers.Map. Keys are experiment names (field 'exp');
    % values are a struct containing parameters in the format used by
    % ExpPP.computeStats.
    
    if exist('fname','var')==0 || isempty(fname);
      expparampath = ExpPP.loadConfigVal('expparampath');
      if isempty(expparampath)
        expparampath = pwd;
      end
      [fname,pname] = uigetfile({'*.csv'},'Select parameter file',expparampath);
      if isequal(fname,0)
        return;
      end
      fname = fullfile(pname,fname);
    else
      if exist(fname,'file')==0
        error('ExpPP:fileNotFound','File ''%s'' not found.',fname);
      end
    end
    
    [pth,~,ext] = fileparts(fname);
    switch ext
      case '.csv'
        fh = fopen(fname);
        if isequal(fh,-1)
          error('ExpPP:errRead','Cannot read file ''%s''.',fname);
        end
        onc = onCleanup(@()fclose(fh));
        try
          line1 = fgetl(fh);
          flds = regexp(line1,',','split');
          Nflds = numel(flds);

          dat = textscan(fh,repmat('%s',1,Nflds),'delimiter',',','collectoutput',true);
          dat = dat{1};
        catch ME
          error('ExpPP:readFile','Error reading parameter file %s: %s',fname,ME.message);
        end
      otherwise
        error('ExpPP:fileFormat','Experiment parameter file must be a csv file.');
    end
    
    tfExp = strcmp('exp',flds);
    if nnz(tfExp)~=1
      error('ExpPP:noExp','Parameter file must include exactly one field ''exp'' which uniquely identifies experiments.');
    end
    
    exps = dat(:,tfExp);
    if numel(exps)~=numel(unique(exps))
      error('ExpPP:noExp','Parameter file contains duplicate values of ''exp'' field.');
    end
    
    assert(size(dat,2)==Nflds,'Expected data size to have %d columns.',Nflds);
    Nrows = size(dat,1);    
    
    fldsNotExp = flds(~tfExp);
    expparammap = containers.Map();
    for i = 1:Nrows
      key = dat{i,tfExp};
      vals = dat(i,~tfExp);
      vals = cellfun(@str2double,vals);
      vals = num2cell(vals);
      s = ExpPP.formParamStruct(fldsNotExp,vals);
      expparammap(key) = s;
    end
    
    fprintf(1,'Loaded %d experiment-specific parameter(s) for %d experiments.\n',...
      Nflds-1,Nrows);
    
    ExpPP.saveConfigVal('expparampath',pth);
  end 

end


methods (Static) % compute stats, massage bouts
  
  function [data,stats] = computeStats(data,varargin)
    % [s,stats] = computeStats(s,varargin)
    
    [computeonlabels,... 
     expparammap,...
     paramdefault.toneframe,...
     paramdefault.fillsize,...
     paramdefault.minboutsize,...
     paramdefault.grabseq_overlap,...
     paramdefault.grabseq_minboutlengths,...
     paramdefault.grabseq_maxdelayframes,...
     paramdefault.grabseq_usefirsthandopen,...
     paramdefault.trialsuccess_subsequentgrabframes,...
     paramdefault.trialsuccess_t0,...
     paramdefault.trialsuccess_t1,...
     paramdefault.trialsuccess_atmouthminframes,...
     paramdefault.trialsuccess_atmouthminframesrange,...
     paramdefault.trialsuccess_allowsuccessfulgrabafterperchreturn,...
     paramdefault.laser_inuse,...
     paramdefault.laser_on,...
     paramdefault.laser_progressmaxdelay,...     
     statparam.framelimits] = myparse(varargin,...
      'computeonlabels',false,... % boolean
      'expparammap',[],... % experiment-specific parameter map
      'toneframe',0,...
      'fillsize',1,...
      'minboutsize',ExpPP.DEFAULT_MINBOUTSIZE,...
      'grabseq_overlap',10,...
      'grabseq_minboutlengths',cell2struct(num2cell(zeros(ExpPP.NBASICBEHAVIORS,1)),ExpPP.BASICBEHAVIORS(:),1),... 
      'grabseq_maxdelayframes',ExpPP.FRMRATE,... % 1 second
      'grabseq_usefirsthandopen',false,...
      'trialsuccess_subsequentgrabframes',100,...
      'trialsuccess_t0',nan,...
      'trialsuccess_t1',nan,...
      'trialsuccess_atmouthminframes',20,...
      'trialsuccess_atmouthminframesrange',200,...
      'trialsuccess_allowsuccessfulgrabafterperchreturn',true,...
      'laser_inuse',0,...
      'laser_on',nan,...
      'laser_progressmaxdelay',nan,...
      'statframelimits',[-inf inf]);    
    
    tfExpParams = ~isempty(expparammap);
    if tfExpParams
      assert(isa(expparammap,'containers.Map'),'Expected ''expparams'' to be a containers.Map object.');
      keys = expparammap.keys;
      vals = expparammap.values;
      assert(all(cellfun(@isstruct,vals)),'Expected ''expparams'' to be a map from experiment-names to parameter structs.');
      
      exps = {data.exp}';
      assert(numel(unique(exps))==numel(exps),'Data contains repeated experiment names.');
      tf = ismember(exps,keys);
      if ~any(tf)
        warningNoTrace('ExpPP:expparams',...
          'Experiment-specific processing parameters were provided, but experiment names do not match any experiments in dataset.');
      elseif ~all(tf)
        warningNoTrace('ExpPP:expparams','Only %d/%d experiments in dataset have experiment-specific parameters. The remaining experiments will use default parameters.',...
          nnz(tf),numel(data));
      end
    end
    validateattributes(statparam.framelimits,{'numeric'},{'numel' 2});

    % check that all basic behaviors were loaded
    for i = 1:numel(ExpPP.BASICBEHAVIORS)
      % use presence/absence of field '<behavior>_t0s' as proxy for whether
      % behavior was loaded
      fld = sprintf('%s_t0s',ExpPP.BASICBEHAVIORS{i}); 
      if ~isfield(data,fld)
        error('ExpPP:behaviorMissing',...
          'Behavior ''%s'' appears to be missing from data.',ExpPP.BASICBEHAVIORS{i});
      end
    end
    for i = 1:numel(data)
      
      % Get/verify the processing parameters for this exp
      exp = data(i).exp;
      param = paramdefault;
      param.exp_specific_params = tfExpParams && expparammap.isKey(exp);
      if param.exp_specific_params
        expparam = expparammap(exp);
        flds = fieldnames(param);
        fldsE = fieldnames(expparam);
        fldsUnrecog = setdiff(fldsE,flds);
        if ~isempty(fldsUnrecog)
          warning('ExpPP:extraFields','Unused/unrecognized experiment-specific parameters: %s',...
            civilizedStringFromCellArrayOfStrings(fldsUnrecog));
        end
        for f = flds(:)',f=f{1}; %#ok<FXSET>
          if isfield(expparam,f)
            param.(f) = expparam.(f);
          end
        end
      end
      param = ExpPP.verifyComputeStatParams(param);
      
      for b = ExpPP.BASICBEHAVIORS, b = b{1}; %#ok<FXSET>
        %%% auto bouts preprocessing

        fldpp = [b '_postprocessed'];
        ypp = data(i).(fldpp);

        [ypp,nfill] = ExpPP.fillInHoles(ypp,param.fillsize);
        [i0s,i1s] = ExpPP.significantBouts(ypp,param.minboutsize.(b));
        
        % keep only after toneFrm
        tf = i0s>=param.toneframe;
        i0s = i0s(tf);
        i1s = i1s(tf);
                
        nfillfld = sprintf('auto_%s_nfill',b);
        t0fld = sprintf('auto_%s_t0s',b);
        t1fld = sprintf('auto_%s_t1s',b);        
        data(i).(nfillfld) = nfill;
        data(i).(t0fld) = i0s;
        data(i).(t1fld) = i1s;        
        
        %%% labl bouts
        if computeonlabels
          t0origFld = sprintf('%s_labl_t0sPos',b);
          t0newFld = sprintf('labl_%s_t0s',b);
          t1origFld = sprintf('%s_labl_t1sPos',b);
          t1newFld = sprintf('labl_%s_t1s',b);
          data(i).(t0newFld) = data(i).(t0origFld);
          data(i).(t1newFld) = data(i).(t1origFld);
        end
      end
      
      %%% auto stats
      lift = ExpPP.hlpBout(data(i).auto_Lift_t0s,data(i).auto_Lift_t1s);
      hand = ExpPP.hlpBout(data(i).auto_Handopen_t0s,data(i).auto_Handopen_t1s);
      grab = ExpPP.hlpBout(data(i).auto_Grab_t0s,data(i).auto_Grab_t1s);
      sup = ExpPP.hlpBout(data(i).auto_Sup_t0s,data(i).auto_Sup_t1s);
      atmouth = ExpPP.hlpBout(data(i).auto_Atmouth_t0s,data(i).auto_Atmouth_t1s);
      chew = ExpPP.hlpBout(data(i).auto_Chew_t0s,data(i).auto_Chew_t1s);      
      
      NLiftPP = numel(data(i).Lift_postprocessed);
      NAtmouthPP = numel(data(i).Atmouth_postprocessed);
      expName = data(i).exp;
      
      [lift0,lift0len,hand0,hand0len,grab0,grab0len,...
       sup0,sup0len,atmouth0,atmouth0len,chew0,chew0len,...
       lift1,lift1len,lift_num,...
       GS00,GS01,trialsuccess,grab_num,GSSS,laserEvent] = ...
        ExpPP.grabSeqAnalysis(param,lift,hand,grab,sup,atmouth,chew,NLiftPP,NAtmouthPP,expName);      
      
      %%% labl stats
      if computeonlabels
        lift_M = ExpPP.hlpBout(data(i).labl_Lift_t0s,data(i).labl_Lift_t1s);
        hand_M = ExpPP.hlpBout(data(i).labl_Handopen_t0s,data(i).labl_Handopen_t1s);
        grab_M = ExpPP.hlpBout(data(i).labl_Grab_t0s,data(i).labl_Grab_t1s);
        sup_M = ExpPP.hlpBout(data(i).labl_Sup_t0s,data(i).labl_Sup_t1s);
        atmouth_M = ExpPP.hlpBout(data(i).labl_Atmouth_t0s,data(i).labl_Atmouth_t1s);
        chew_M = ExpPP.hlpBout(data(i).labl_Chew_t0s,data(i).labl_Chew_t1s);
        
        assert(NLiftPP==numel(data(i).Lift_labl_label));
        assert(NAtmouthPP==numel(data(i).Atmouth_labl_label));
        
        % parameters for labl grabSeqAnls
        assert(isstruct(param));
        param_M = param;
        param_M.toneframe = 0;
        param_M.fillsize = 0;
        param_M.minboutsize.Lift = 0;
        param_M.minboutsize.Handopen = 0;
        param_M.minboutsize.Grab = 0;
        param_M.minboutsize.Sup = 0;
        param_M.minboutsize.Atmouth = 0;
        param_M.minboutsize.Chew = 0;
        param_M.grabseq_minboutlengths.Lift = 0;
        param_M.grabseq_minboutlengths.Handopen = 0;
        param_M.grabseq_minboutlengths.Grab = 0;
        param_M.grabseq_minboutlengths.Sup = 0;
        param_M.grabseq_minboutlengths.Atmouth = 0;
        param_M.grabseq_minboutlengths.Chew = 0;
        param_M.trialsuccess_atmouthminframes = 1;
        
        [lift0_M,lift0len_M,hand0_M,hand0len_M,grab0_M,grab0len_M,...
          sup0_M,sup0len_M,atmouth0_M,atmouth0len_M,chew0_M,chew0len_M,...
          lift1_M,lift1len_M,lift_num_M,...
          GS00_M,GS01_M,trialsuccess_M,grab_num_M,GSSS_M,laserEvent_M] = ...
          ExpPP.grabSeqAnalysis(param_M,lift_M,hand_M,grab_M,sup_M,atmouth_M,chew_M,NLiftPP,NAtmouthPP,expName);
      end
      
      %%% Assign results to struct fields
      
      % assign auto-params to struct
      paramfns = fieldnames(param);
      for fn = paramfns(:)', fn = fn{1}; %#ok<FXSET>
        val = param.(fn);
        if isstruct(val)
          % special case: should recurse
          valfieldnames = fieldnames(val);
          for valfn = valfieldnames(:)', valfn = valfn{1}; %#ok<FXSET>
            newfn = ['param_' fn '__' valfn];
            data(i).(newfn) = val.(valfn);
          end
        else
          newfn = ['param_' fn];
          data(i).(newfn) = val;
        end
      end
      
      % don't bother with labl-params for now
      
      nstDataIAssign('auto',...
        lift0,lift0len,hand0,hand0len,grab0,grab0len,...
        sup0,sup0len,atmouth0,atmouth0len,chew0,chew0len,...
        lift1,lift1len,lift_num,...
        GS00,GS01,trialsuccess,grab_num,GSSS,laserEvent);
      
      if computeonlabels
        nstDataIAssign('labl',...
           lift0_M,lift0len_M,hand0_M,hand0len_M,grab0_M,grab0len_M,...
           sup0_M,sup0len_M,atmouth0_M,atmouth0len_M,chew0_M,chew0len_M,...
           lift1_M,lift1len_M,lift_num_M,...
          GS00_M,GS01_M,trialsuccess_M,grab_num_M,GSSS_M,laserEvent_M);
      end
      
    end
    
    function nstDataIAssign(pfix,...
        zz_lift0,zz_lift0len,zz_hand0,zz_hand0len,zz_grab0,zz_grab0len,...
        zz_sup0,zz_sup0len,zz_atmouth0,zz_atmouth0len,zz_chew0,zz_chew0len,...
        zz_lift1,zz_lift1len,zz_lift_num,...
        zz_GS00,zz_GS01,zz_trialsuccess,zz_grab_num,zz_GSSS,zz_laserEvent)
      
      data(i).([pfix '_Lift_0']) = zz_lift0;
      data(i).([pfix '_Lift_0_bl']) = zz_lift0len;
      data(i).([pfix '_Handopen_0']) = zz_hand0;
      data(i).([pfix '_Handopen_0_bl']) = zz_hand0len;
      data(i).([pfix '_Grab_0']) = zz_grab0;
      data(i).([pfix '_Grab_0_bl']) = zz_grab0len;
      data(i).([pfix '_Sup_0']) = zz_sup0;
      data(i).([pfix '_Sup_0_bl']) = zz_sup0len;
      data(i).([pfix '_Atmouth_0']) = zz_atmouth0;
      data(i).([pfix '_Atmouth_0_bl']) = zz_atmouth0len;
      data(i).([pfix '_Chew_0']) = zz_chew0;
      data(i).([pfix '_Chew_0_bl']) = zz_chew0len;
      
      data(i).([pfix '_Lift_1']) = zz_lift1;
      data(i).([pfix '_Lift_1_bl']) = zz_lift1len;
      data(i).([pfix '_Lift_num']) = zz_lift_num;
      data(i).([pfix '_del_Lift1_Lift0']) = zz_lift1-zz_lift0;
      
      % GS00
      % The following assignments work even if no Grab Seq was found for
      % GS00, because the default GS00 is created with all nan prop vals
      data(i).([pfix '_GS00_Lift_0']) = zz_GS00.lift(1);
      data(i).([pfix '_GS00_Handopen_0']) = zz_GS00.handopen(1);
      data(i).([pfix '_GS00_Grab_0']) = zz_GS00.grab(1);
      
      data(i).([pfix '_GS00_Lift_bl']) = diff(zz_GS00.lift);
      data(i).([pfix '_GS00_Handopen_bl']) = diff(zz_GS00.handopen);
      data(i).([pfix '_GS00_Grab_bl']) = diff(zz_GS00.grab);
      
      data(i).([pfix '_del_GS00_Handopen_Lift']) = zz_GS00.handopen(1)-zz_GS00.lift(1);
      data(i).([pfix '_del_GS00_Grab_Lift']) = zz_GS00.grab(1)-zz_GS00.lift(1);
      data(i).([pfix '_del_GS00_Grab_Handopen']) = zz_GS00.grab(1)-zz_GS00.handopen(1);
      
      % GS01
      % Same comment as for GS00
      data(i).([pfix '_GS01_Lift_0']) = zz_GS01.lift(1);
      data(i).([pfix '_GS01_Handopen_0']) = zz_GS01.handopen(1);
      data(i).([pfix '_GS01_Grab_0']) = zz_GS01.grab(1);            
      
      data(i).([pfix '_del_GS01_Grab_GS00_Grab']) = zz_GS01.grab(1)-zz_GS00.grab(1);

      % GSSS
      data(i).([pfix '_Grab_success']) = zz_trialsuccess.tf;
      data(i).([pfix '_Grab_successtype']) = zz_trialsuccess.type;
      data(i).([pfix '_Grab_num']) = zz_grab_num;
            
      % The following assignments work even if there is no successful Grab
      % Seq, because the default GSSS is created with all nan prop vals.
      
      data(i).([pfix '_GSSS_Lift_0']) = zz_GSSS.lift(1);
      data(i).([pfix '_GSSS_Handopen_0']) = zz_GSSS.handopen(1);
      data(i).([pfix '_GSSS_Grab_0']) = zz_GSSS.grab(1);
      data(i).([pfix '_GSSS_Sup_0']) = zz_GSSS.sup(1);
      data(i).([pfix '_GSSS_Atmouth_0']) = zz_GSSS.atmouth(1);
      data(i).([pfix '_GSSS_Chew_0']) = zz_GSSS.chew(1);
      data(i).([pfix '_GSSS_Lift_bl']) = diff(zz_GSSS.lift);
      data(i).([pfix '_GSSS_Handopen_bl']) = diff(zz_GSSS.handopen);
      data(i).([pfix '_GSSS_Grab_bl']) = diff(zz_GSSS.grab);
      data(i).([pfix '_GSSS_Sup_bl']) = diff(zz_GSSS.sup);
      data(i).([pfix '_GSSS_Atmouth_bl']) = diff(zz_GSSS.atmouth);
      data(i).([pfix '_GSSS_Chew_bl']) = diff(zz_GSSS.chew);
      
      data(i).([pfix '_del_GSSS_Handopen_Lift']) = zz_GSSS.handopen(1)-zz_GSSS.lift(1);
      data(i).([pfix '_del_GSSS_Grab_Lift']) = zz_GSSS.grab(1)-zz_GSSS.lift(1);
      data(i).([pfix '_del_GSSS_Grab_Handopen']) = zz_GSSS.grab(1)-zz_GSSS.handopen(1);
      data(i).([pfix '_del_GSSS_Sup_Grab']) = zz_GSSS.sup(1)-zz_GSSS.grab(1);
      data(i).([pfix '_del_GSSS_Atmouth_Lift']) = zz_GSSS.atmouth(1)-zz_GSSS.lift(1);
      data(i).([pfix '_del_GSSS_Atmouth_Grab']) = zz_GSSS.atmouth(1)-zz_GSSS.grab(1);
      data(i).([pfix '_del_GSSS_Atmouth_Sup']) = zz_GSSS.atmouth(1)-zz_GSSS.sup(1);
      data(i).([pfix '_del_GSSS_Chew_Atmouth']) = zz_GSSS.chew(1)-zz_GSSS.atmouth(1);
      
      % LaserEvent
      data(i).([pfix '_laser_intrpt_nbeh']) = zz_laserEvent.intrptNbeh;
      data(i).([pfix '_laser_intrpt_behName']) = zz_laserEvent.intrptBeh0Name;
      data(i).([pfix '_laser_intrpt_beht0']) = zz_laserEvent.intrptBeh0t0;
      data(i).([pfix '_laser_intrpt_behtdel']) = zz_laserEvent.intrptBeh0tdel;
      data(i).([pfix '_laser_intrpt_isPhantom']) = zz_laserEvent.intrptBeh0IsPhantom;
      data(i).([pfix '_laser_intrpt_progress']) = zz_laserEvent.intrptTFProgress;
      data(i).([pfix '_laser_intrpt_behProgName']) = zz_laserEvent.intrptBehProgName;
      data(i).([pfix '_laser_intrpt_behProgt0']) = zz_laserEvent.intrptBehProgt0;
    end      
      
    function nstAssignFLDS(zzPfix)
      FLDS.([zzPfix '_Lift_0']) = lims;
      FLDS.([zzPfix '_Handopen_0']) = lims;
      FLDS.([zzPfix '_Grab_0']) = lims;
      FLDS.([zzPfix '_Sup_0']) = lims;
      FLDS.([zzPfix '_Atmouth_0']) = lims;
      FLDS.([zzPfix '_Chew_0']) = lims;
      FLDS.([zzPfix '_Lift_0_bl']) = @(x)lclInRange(x.([zzPfix '_Lift_0']),lims);
      FLDS.([zzPfix '_Handopen_0_bl']) = @(x)lclInRange(x.([zzPfix '_Handopen_0']),lims);
      FLDS.([zzPfix '_Grab_0_bl']) = @(x)lclInRange(x.([zzPfix '_Grab_0']),lims);
      FLDS.([zzPfix '_Sup_0_bl']) = @(x)lclInRange(x.([zzPfix '_Sup_0']),lims);
      FLDS.([zzPfix '_Atmouth_0_bl']) = @(x)lclInRange(x.([zzPfix '_Atmouth_0']),lims);
      FLDS.([zzPfix '_Chew_0_bl']) = @(x)lclInRange(x.([zzPfix '_Chew_0']),lims);
      FLDS.([zzPfix '_Lift_1']) = lims;
      FLDS.([zzPfix '_Lift_1_bl']) = @(x)lclInRange(x.([zzPfix '_Lift_1']),lims);
      
      FLDS.([zzPfix '_GS00_Lift_bl']) = @(x)lclInRange(x.([zzPfix '_GS00_Lift_0']),lims);
      FLDS.([zzPfix '_GS00_Handopen_bl']) = @(x)lclInRange(x.([zzPfix '_GS00_Handopen_0']),lims);
      FLDS.([zzPfix '_GS00_Grab_bl']) = @(x)lclInRange(x.([zzPfix '_GS00_Grab_0']),lims);
      
      FLDS.([zzPfix '_GSSS_Lift_bl']) = @(x)lclInRange(x.([zzPfix '_GSSS_Lift_0']),lims);
      FLDS.([zzPfix '_GSSS_Handopen_bl']) = @(x)lclInRange(x.([zzPfix '_GSSS_Handopen_0']),lims);
      FLDS.([zzPfix '_GSSS_Grab_bl']) = @(x)lclInRange(x.([zzPfix '_GSSS_Grab_0']),lims);
      FLDS.([zzPfix '_GSSS_Sup_bl']) = @(x)lclInRange(x.([zzPfix '_GSSS_Sup_0']),lims);
      FLDS.([zzPfix '_GSSS_Atmouth_bl']) = @(x)lclInRange(x.([zzPfix '_GSSS_Atmouth_0']),lims);
      FLDS.([zzPfix '_GSSS_Chew_bl']) = @(x)lclInRange(x.([zzPfix '_GSSS_Chew_0']),lims);
      
      FLDS.([zzPfix '_del_Lift1_Lift0']) = @(x)lclInRange(x.([zzPfix '_Lift_1']),lims) && lclInRange(x.([zzPfix '_Lift_0']),lims);
      
      FLDS.([zzPfix '_del_GS00_Handopen_Lift']) = @(x)lclInRange(x.([zzPfix '_GS00_Lift_0']),lims) && lclInRange(x.([zzPfix '_GS00_Handopen_0']),lims);
      FLDS.([zzPfix '_del_GS00_Grab_Lift']) = @(x)lclInRange(x.([zzPfix '_GS00_Lift_0']),lims) && lclInRange(x.([zzPfix '_GS00_Grab_0']),lims);
      FLDS.([zzPfix '_del_GS00_Grab_Handopen']) = @(x)lclInRange(x.([zzPfix '_GS00_Grab_0']),lims) && lclInRange(x.([zzPfix '_GS00_Handopen_0']),lims);
      
      FLDS.([zzPfix '_del_GS01_Grab_GS00_Grab']) = @(x)lclInRange(x.([zzPfix '_GS01_Grab_0']),lims) && lclInRange(x.([zzPfix '_GS00_Grab_0']),lims);
      
      FLDS.([zzPfix '_del_GSSS_Handopen_Lift']) = @(x)lclInRange(x.([zzPfix '_GSSS_Lift_0']),lims) && lclInRange(x.([zzPfix '_GSSS_Handopen_0']),lims);
      FLDS.([zzPfix '_del_GSSS_Grab_Lift']) = @(x)lclInRange(x.([zzPfix '_GSSS_Grab_0']),lims) && lclInRange(x.([zzPfix '_GSSS_Lift_0']),lims);
      FLDS.([zzPfix '_del_GSSS_Grab_Handopen']) = @(x)lclInRange(x.([zzPfix '_GSSS_Grab_0']),lims) && lclInRange(x.([zzPfix '_GSSS_Handopen_0']),lims);
      FLDS.([zzPfix '_del_GSSS_Sup_Grab']) = @(x)lclInRange(x.([zzPfix '_GSSS_Sup_0']),lims) && lclInRange(x.([zzPfix '_GSSS_Grab_0']),lims);
      FLDS.([zzPfix '_del_GSSS_Atmouth_Lift']) = @(x)lclInRange(x.([zzPfix '_GSSS_Atmouth_0']),lims) && lclInRange(x.([zzPfix '_GSSS_Lift_0']),lims);
      FLDS.([zzPfix '_del_GSSS_Atmouth_Grab']) = @(x)lclInRange(x.([zzPfix '_GSSS_Atmouth_0']),lims) && lclInRange(x.([zzPfix '_GSSS_Grab_0']),lims);
      FLDS.([zzPfix '_del_GSSS_Atmouth_Sup']) = @(x)lclInRange(x.([zzPfix '_GSSS_Atmouth_0']),lims) && lclInRange(x.([zzPfix '_GSSS_Sup_0']),lims);
      FLDS.([zzPfix '_del_GSSS_Chew_Atmouth']) = @(x)lclInRange(x.([zzPfix '_GSSS_Chew_0']),lims) && lclInRange(x.([zzPfix '_GSSS_Atmouth_0']),lims);      
    end
    
     if nargout>1
      % compute descriptive stats
            
      lims = statparam.framelimits;
      % FLDS structure
      % Fieldname: Field in data struct for which to compute stats.
      % Val: inlier limits (1x2 numeric array), or predicate fcn handle. 
      %   If a 1x2 numerical vec, then the value of the field must lie in 
      % the specified range to be considered an inlier. If a fcn handle, 
      % then the fcn handle is evaluated on each data element and should 
      % return true for inliers and false for outliers.
      %
      % The point is that fields which represent frames/times, such as
      % auto_Lift_0 or "first Lift", have inlier status specified by
      % inclusion in a simple range. For fields like auto_Lift_0_bl, or
      % "first Lift bout length", inlier status should probably be governed
      % by the inlier status of auto_Lift_0... 
      %
      % The framelimits thing is probably not totally thought out but ok
      % for now.
      
      FLDS = struct();
      nstAssignFLDS('auto');
      nstAssignFLDS('labl');
      
      % SETS structure
      % Fields: groups under which to compute stats
      % Values: Predicate fcns for group membership
      %
      % Default groups/sets
      SETS.alldata = @(x)true;
      SETS.auto_Grab_successful = @(x)x.auto_Grab_success;
      SETS.auto_Grab_unsuccessful = @(x)~x.auto_Grab_success;
      SETS.auto_Grab_successful_1grab = @(x)strcmp(x.auto_Grab_successtype,'successful_1grab');
      SETS.auto_Grab_successful_2plusgrab = @(x)strcmp(x.auto_Grab_successtype,'successful_2plusgrab');
      if computeonlabels
        SETS.labl_Grab_successful = @(x)x.labl_Grab_success;
        SETS.labl_Grab_unsuccessful = @(x)~x.labl_Grab_success;
        SETS.labl_Grab_successful_1grab = @(x)strcmp(x.labl_Grab_successtype,'successful_1grab');
        SETS.labl_Grab_successful_2plusgrab = @(x)strcmp(x.labl_Grab_successtype,'successful_2plusgrab');
      end
      % Custom groups/sets
      [data,customgrps] = ExpPP.ensureCustomGroupInit(data);
      Ncustomgrps = numel(customgrps);
      for iGrp = 1:Ncustomgrps
        grpname = customgrps{iGrp};
        fldtf = sprintf(ExpPP.CUSTOMGROUPNAME_PAT,grpname);
        SETS.(fldtf) = @(x)x.(fldtf);
      end
      
      setNames = fieldnames(SETS);
      Nset = numel(setNames);
      setMat = false(numel(data),Nset);
      for iSet = 1:Nset
        setMat(:,iSet) = arrayfun(SETS.(setNames{iSet}),data);
      end
      
      stats = ExpPP.computeDescriptiveStatistics(data,FLDS,setMat,setNames);
    end
  end
  
  function param = verifyComputeStatParams(param)
    validateattributes(param.toneframe,{'numeric'},{'integer' 'scalar'});
    validateattributes(param.fillsize,{'numeric'},{'nonnegative' 'integer' 'scalar'});
    assert(isstruct(param.minboutsize)&&all(isfield(param.minboutsize,ExpPP.BASICBEHAVIORS)),'Parameter ''minboutsize'' must encompass all basic behaviors.');
    validateattributes(param.grabseq_overlap,{'numeric'},{'nonnegative' 'integer' 'scalar'});
    
    assert(isstruct(param.grabseq_minboutlengths)&&all(isfield(param.grabseq_minboutlengths,ExpPP.BASICBEHAVIORS)),'Parameter ''grabseq_minboutlengths'' must encompass all basic behaviors.');
    unrecogFlds = setdiff(fieldnames(param.grabseq_minboutlengths),ExpPP.BASICBEHAVIORS);
    if ~isempty(unrecogFlds)
      warnstr = civilizedStringFromCellArrayOfStrings(unrecogFlds);
      warning('ExpPP:unrecognizedFields',...
        'Unrecognized behaviors in min-bout-lengths structure: %s',...
        warnstr);
    end
    validateattributes(param.grabseq_maxdelayframes,{'numeric'},{'nonnegative' 'integer' 'scalar'});
    validateattributes(param.grabseq_usefirsthandopen,{'numeric' 'logical'},{'scalar'});
    
    validateattributes(param.trialsuccess_subsequentgrabframes,{'numeric'},{'nonnegative' 'integer' 'scalar'});
    validateattributes(param.trialsuccess_t0,{'numeric'},{'scalar'});
    validateattributes(param.trialsuccess_t1,{'numeric'},{'scalar'});
    validateattributes(param.trialsuccess_allowsuccessfulgrabafterperchreturn,{'numeric' 'logical'},{'scalar'});
    
    validateattributes(param.laser_inuse,{'numeric' 'logical'},{'scalar'});
    validateattributes(param.laser_on,{'numeric'},{'integer' 'scalar'});
    validateattributes(param.laser_progressmaxdelay,{'numeric'},{'scalar'});
    
    if isnan(param.trialsuccess_t0)
      param.trialsuccess_t0 = -inf;
    end
    if isnan(param.trialsuccess_t1)
      param.trialsuccess_t1 = inf;
    end    
  end
      
  function [lift0,lift0len,hand0,hand0len,grab0,grab0len,...
            sup0,sup0len,atmouth0,atmouth0len,chew0,chew0len,...
            lift1,lift1len,lift_num,...
            GS00,GS01,trialsuccess,grab_num,GSSS,laserEvent] = ...
      grabSeqAnalysis(param,lift,hand,grab,sup,atmouth,chew,NLiftPP,NAtmouthPP,expName)
    
    % First <behavior>: first bout after tone
    [lift0,lift0len] = ExpPP.firstBout(lift);
    [hand0,hand0len] = ExpPP.firstBout(hand);
    [grab0,grab0len] = ExpPP.firstBout(grab);
    [sup0,sup0len] = ExpPP.firstBout(sup);
    [atmouth0,atmouth0len] = ExpPP.firstBout(atmouth);
    [chew0,chew0len] = ExpPP.firstBout(chew);
    
    % Global Last Lift/Num Lift
    ref = nanmin([hand0 grab0 sup0 atmouth0 chew0 NLiftPP+1]);
    [lift1,lift1len] = ExpPP.lastBoutBeforeRef(lift,ref);
    lift_num = ExpPP.numBoutStartsBetweenRefs(lift.t0s,param.toneframe,ref); % lift.t0s should always be >= param.toneframe anyway
    % Alternative to anove condition: if ~isnan(hand0), use hand0 as ref;
    % otherwise if ~isnan(grab0); etc
    
    % Setup for
    % * HG seq-finding
    % * Successful grab seq-finding
    Iatmouth = ExpPP.boutEndpoints2IndicatorVec(NAtmouthPP,atmouth.t0s,atmouth.t1s);
    Iatmouthtf = Iatmouth>0;
    
    % GSSS: Grab Sequence, Successful
    %
    % - A successful grab seq contains all behaviors (Lift-Hand-Grab-Sup-Atmouth-Chew),
    % in order
    % - There must be a minimum number of frames of Atmouth within a
    % prescribed range. The point here is to exclude sequences where the
    % atmouth is very brief, _unless_ there are several brief atmouths in
    % rapid succession.
    % - if trialsuccess_t0/t1 are specified, the start of the Atmouth bout
    % must occur in range [trialsuccess_t0,trialsuccess_t1)
    %
    % For each Chew bout (starting front-to-back)
    %   Look for immediately preceding Atmouth (with atmouth constraints,
    %   ie duration and trialsuccess_t0/t1), can overlap slightly etc
    %     If found, look for immediately preceding Sup, can overlap slightly etc
    %       ...
    %         If found all the way up to and including Lift, we are done
    %         and have found the first successful Grab seq.
    % Otherwise no successful GS found
    
    trialsuccess = struct();
    trialsuccess.tf = false;
    GSSS = ExpPPGrabSeq;
    
    tfchews = chew.dts >= param.grabseq_minboutlengths.Chew;
    ichews = find(tfchews(:)');
    for ichew = ichews
      % loop over chews in order (front-to-back) as we are looking for
      % first successful Chew. If/once we find it we will break.
      
      tfatmos = ExpPP.boutsOfB_RoughlyBeforeBout(atmouth.t0s,chew.t0s(ichew),chew.t1s(ichew),param.grabseq_overlap) & ...
        atmouth.dts >= param.grabseq_minboutlengths.Atmouth & ...
        atmouth.t0s >= param.trialsuccess_t0 & ...
        atmouth.t0s < param.trialsuccess_t1 & ...
        ExpPP.hlpAtmouthBoutConstraint(atmouth.t0s,Iatmouthtf,...
        param.trialsuccess_atmouthminframesrange,param.trialsuccess_atmouthminframes);
      % require 20/200 frames starting from atmouth bout to be atmouth
      iatmo = find(tfatmos(:)',1,'last');
      if ~isempty(iatmo)
        tfsup = ExpPP.boutsOfB_RoughlyBeforeBout(sup.t0s,atmouth.t0s(iatmo),atmouth.t1s(iatmo),param.grabseq_overlap) & ...
          sup.dts >= param.grabseq_minboutlengths.Sup;
        isup = find(tfsup(:)',1,'last');
        if ~isempty(isup)
          tfgrab = ExpPP.boutsOfB_RoughlyBeforeBout(grab.t0s,sup.t0s(isup),sup.t1s(isup),param.grabseq_overlap) & ...
            grab.dts >= param.grabseq_minboutlengths.Grab;
          igrab = find(tfgrab(:)',1,'last');
          if ~isempty(igrab)
            tfhand = ExpPP.boutsOfB_RoughlyBeforeBout(hand.t0s,grab.t0s(igrab),grab.t1s(igrab),param.grabseq_overlap) & ...
              hand.dts >= param.grabseq_minboutlengths.Handopen;
            ihand = find(tfhand(:)',1,'last');
            if ~isempty(ihand)
              tflift = ExpPP.boutsOfB_RoughlyBeforeBout(lift.t0s,hand.t0s(ihand),hand.t1s(ihand),param.grabseq_overlap) & ...
                lift.dts >= param.grabseq_minboutlengths.Lift;
              ilift = find(tflift(:)',1,'last');
              if ~isempty(ilift)
                % successful seq!
                
                if param.grabseq_usefirsthandopen
                  % go back and look for multiple Hand bouts in range
                  % [lift.t0s(ilift),hand.t0s(ihand)]. If there is more
                  % than one bout, take the first.
                  
                  % due to fudgefac, hand.t0s(ihand) could occur before
                  % lift.t0s(ilift). Note, we could go a step further and
                  % account for the fudgefac by using something like
                  % lift.t0s(ilift)-param.grabseq_overlap.
                  
                  if igrab>1
                    % Previous Grab serves as additional constraint
                    ref0 = min(max(lift.t0s(ilift),grab.t0s(igrab-1)+1),hand.t0s(ihand));
                  else
                    ref0 = min(lift.t0s(ilift),hand.t0s(ihand));
                  end
                  ref1 = hand.t0s(ihand)+1; % +1 due to boutsBetweenRefs API
                  tftmp = ExpPP.boutsBetweenRefs(hand.t0s,ref0,ref1);
                  assert(nnz(tftmp)>=1);
                  if nnz(tftmp)>1
                    warningNoTrace('ExpPP:multiHand','Multiple Handopen bouts in GSSS: %s.',expName);
                  end
                  ihanduse = find(tftmp,1);
                else
                  ihanduse = ihand;
                end
                
                GSSS.lift = [lift.t0s(ilift) lift.t1s(ilift)];
                %                   GSSS.handopen = [hand.t0s(ihand) hand.t1s(ihand)];
                GSSS.handopen = [hand.t0s(ihanduse) hand.t1s(ihanduse)];
                GSSS.grab = [grab.t0s(igrab) grab.t1s(igrab)];
                GSSS.sup = [sup.t0s(isup) sup.t1s(isup)];
                GSSS.atmouth = [atmouth.t0s(iatmo) atmouth.t1s(iatmo)];
                GSSS.chew = [chew.t0s(ichew) chew.t1s(ichew)];
                
                trialsuccess.tf = true;
                break;
              end
            end
          end
        end
      end
    end
    
    % GS00: Grab (Hand-Grab) Sequence, First
    % - An HG seq contains a Hand then Grab. The Hand is defined to be the last
    % Handbout roughly before the Grab.
    % - In future, HG may optionally contain a Sup, Atmouth, and Chew
    % - Each Hand/Grab bout my participate in at most one HG seq
    %
    % Note, a GS00 will always be found if a GSSS is present, but not
    % vice versa.
    %
    % For now, here are possibilities for bouts in GS00 overlapping with
    % bouts in GSSS:
    % - Shared Hand and Grab. OK.
    % - Shared Hand, different Grab. OK.
    % - Different Hand, same Grab. Impossible, the Hand is always the last etc.
    
    GS00 = ExpPPGrabSeq;
    GS01 = ExpPPGrabSeq;
    state = 'found_none'; % {found_none, found_onegrab, found_twograbs}
    tfgrabs = grab.dts >= param.grabseq_minboutlengths.Grab;
    igrabs = find(tfgrabs(:)');
    return2perch = nan;
    for igrab = igrabs
      % loop over grabs in order (front-to-back) as we are looking for
      % first HG seq, then second.
      
      tfhand = ExpPP.boutsOfB_RoughlyBeforeBout(hand.t0s,grab.t0s(igrab),grab.t1s(igrab),param.grabseq_overlap) & ...
        hand.dts >= param.grabseq_minboutlengths.Handopen;
      ihand = find(tfhand(:)',1,'last');
      if ~isempty(ihand)
        % currently we require a Lift as well
        tflift = ExpPP.boutsOfB_RoughlyBeforeBout(lift.t0s,hand.t0s(ihand),hand.t1s(ihand),param.grabseq_overlap) & ...
          lift.dts >= param.grabseq_minboutlengths.Lift;
        ilift = find(tflift(:)',1,'last');
        if ~isempty(ilift)
          
          if param.grabseq_usefirsthandopen
            % See comment in GSSS-finding
            if igrab>1
              % Previous Grab serves as additional constraint
              ref0 = min(max(lift.t0s(ilift),grab.t0s(igrab-1)+1),hand.t0s(ihand));
            else
              ref0 = min(lift.t0s(ilift),hand.t0s(ihand));
            end
            ref1 = hand.t0s(ihand)+1; % +1 due to boutsBetweenRefs API
            tftmp = ExpPP.boutsBetweenRefs(hand.t0s,ref0,ref1);
            assert(nnz(tftmp)>=1);
            if nnz(tftmp)>1
              warningNoTrace('ExpPP:multiHand','Multiple Handopen bouts in GS00/01: %s.',expName);
            end
            ihanduse = find(tftmp,1);
          else
            ihanduse = ihand;
          end
          
          switch state
            case 'found_none'
              assert(all(isnan(GS00.grab)));
              
              GS00.lift = [lift.t0s(ilift) lift.t1s(ilift)];
              GS00.handopen = [hand.t0s(ihanduse) hand.t1s(ihanduse)];
              GS00.grab = [grab.t0s(igrab) grab.t1s(igrab)];              
  
              tfReturnLift = ExpPP.boutsBetweenRefs(lift.t0s,GS00.grab(1)+1,inf);
              if any(tfReturnLift)
                return2perch = lift.t0s(tfReturnLift);
                return2perch = return2perch(1); % first lift that occurs after start of GS00.grab
              else
                % return2perch initialized to nan
              end
              
              state = 'found_onegrab';
            case 'found_onegrab'
              assert(all(isnan(GS01.grab)));
              
              if hand.t0s(ihanduse)==GS00.handopen(1)
                % Another grab found with same preceding handopen; this is
                % not a new Grab seq. 
              else
                GS01.lift = [lift.t0s(ilift) lift.t1s(ilift)];
                GS01.handopen = [hand.t0s(ihanduse) hand.t1s(ihanduse)];
                GS01.grab = [grab.t0s(igrab) grab.t1s(igrab)];
                
                state = 'found_twograb'; %#ok<NASGU>
                break;
              end
              
            otherwise
              assert(false);
          end
        end
      end
    end
    
    % Grab_num, Grab success, Grab successtype
    if trialsuccess.tf
      % Conceptually we want to count the Grabs that occur between lift0
      % and the successfulGrabSeq (sGS). The natural reference point in
      % sGS is the sGS.Grab itself. Previously we had the notion of
      % counting to the sGS.Atmouth. Typically, there will be no
      % difference between these two counts, as the sGS.Grab will
      % be the last Grab before sGS.Atmouth.
      %
      % Edge cases:
      % 1. sGS.Grab occurs after sGS.Atmouth, due to allowed overlap. In
      % this case, using sGS.Grab is desired for the second ref point.
      % 2. sGS.Grab occurs before sGS.Atmouth, but is not the last Grab
      % before sGS.Atmouth, because of the placement of the intervening
      % Sup. In other words, imagine the sequence Grab-Sup-Grab-Atmouth.
      
      % - GSSS is present => GS00 is present.
      % - The Grab in GS00 should definitely be counted.
      % - Due to overlap, GS00.grab(1) could theoretically occur before
      % lift0.
      assert(~isnan(lift0) && ~isnan(GS00.grab(1)));
      ref1 = min([lift0 GS00.grab(1)]);
      %ref2 = max(successfulGrabSeq.grab(2),successfulGrabSeq.atmouth(1));
      grab_num = ExpPP.numBoutStartsBetweenRefs(grab.t0s,ref1,GSSS.grab(2));
    else
      if isnan(lift0) % No lifts
        grab_num = ExpPP.numBoutStartsBetweenRefs(grab.t0s,param.toneframe,inf); % all grab.t0s should be after param.toneframe anyway
      else
        % if present, we want to include grab in GS00.
        ref1 = nanmin([lift0 GS00.grab(1)]);
        grab_num = ExpPP.numBoutStartsBetweenRefs(grab.t0s,ref1,inf);
      end
    end
    % massage trialsuccess.tf based on
    % param.trialsuccess_allowsuccessfulgrabafterperchreturn and return2perch
    if trialsuccess.tf && ~param.trialsuccess_allowsuccessfulgrabafterperchreturn ...
        && ~isnan(return2perch)
      assert(return2perch>GS00.grab(1));
      if GSSS.grab(1)>return2perch
        % This implies GSSS.Grab(1)>GS00.Grab(1)
        % The successful grab is occuring after the return2perch.
        trialsuccess.tf = false;
      end
    end
    % successfrm, successtype
    if trialsuccess.tf
      assert(grab_num>0,'Trial successful, expect at least one Grab.');
      if grab_num==1
        trialsuccess.type = 'successful_1grab'; % TODO: trialsuccess.type is a 3-way enum
      else
        trialsuccess.type = 'successful_2plusgrab';
      end
    else
      trialsuccess.type = 'unsuccessful';
    end
    
    laserEvent = ExpPPLaserEvent;
    if param.laser_inuse
      laserEvent.laserOn = param.laser_on;
      
      allbouts = ExpPP.aggregateBouts('lift',lift,'hand',hand,'grab',grab,...
        'sup',sup,'atmouth',atmouth,'chew',chew);
            
      % find all real- and phantom-interupted bouts
      allt0s = getstructarrayfield(allbouts,'t0');
      allt1s = getstructarrayfield(allbouts,'t1');
      allbehs = {allbouts.name}';
      assert(isequal(numel(allt0s),numel(allt1s),numel(allbehs)));
      tfInterrupt = allt0s<=laserEvent.laserOn & laserEvent.laserOn<allt1s;
      tfPhantomInterrupt = allt0s<=laserEvent.laserOn;
      laserEvent.intrptNbeh = nnz(tfInterrupt);
      
      if any(tfInterrupt) || any(tfPhantomInterrupt)
        % at least one bout interrupted, analyze further
        
        if any(tfInterrupt)
          iBoutInterrupt = find(tfInterrupt);
          laserEvent.intrptBeh0IsPhantom = false;
        else
          assert(any(tfPhantomInterrupt));
          iBoutInterrupt = find(tfPhantomInterrupt);
          laserEvent.intrptBeh0IsPhantom = true;
        end

        % Reduce to a single interruption if necessary
        if numel(iBoutInterrupt)>1
          % AL: don't throw the warning for now since this codepath applies
          % to phantom bouts as well
          
%           warningNoTrace('ExpPP:grabSeqAnalysis',...
%             'Exp %s: Laser on at frame=%d interrupts multiple (%d) behaviors. Analysis will use last/latest interrupted behavior.',...
%             expName,laserEvent.laserOn,laserEvent.intrptNbeh);
          
          t0Interrupt = allt0s(iBoutInterrupt);
          [~,idxtmp] = sort(t0Interrupt(:),1,'descend');
          iBoutInterrupt = iBoutInterrupt(idxtmp(1));
        end
        
        assert(isscalar(iBoutInterrupt));
        boutIntrpt = allbouts(iBoutInterrupt);
        laserEvent.intrptBeh0Name = boutIntrpt.name;
        laserEvent.intrptBeh0t0 = boutIntrpt.t0;
        
        % Look for 'Progression' within window [laserOn,laserOn+maxDelay):
        % * Interrupted behavior (beh0) must end
        % * Another behavior (behProgress), different from beh0, must start
        % AG QUESTION: does beh0 have to end? if so, does it have to end before behProg starts?
        
        laserWinEnd = laserEvent.laserOn + param.laser_progressmaxdelay;
        tfBeh0Ended = boutIntrpt.t1 < laserWinEnd;
        tfProgBeh = laserEvent.laserOn<allt0s & allt0s<laserWinEnd & ...
                    ~strcmp(allbehs,laserEvent.intrptBeh0Name); 
        iBoutProg = find(tfProgBeh);
        nBoutProg = numel(iBoutProg);
        
        laserEvent.intrptTFProgress = tfBeh0Ended && nBoutProg>0;
        if laserEvent.intrptTFProgress
          if nBoutProg>1
            %warning('ExpPP:grabSeqAnalysis','Laser on at frame=%d, multiple progression behaviors. Analysis will use earliest progression behavior.',...
            %  laserEvent.laserOn);
            % Use first/earliest progression behavior
            t0Started = allt0s(iBoutProg);
            [~,idxtmp] = sort(t0Started);
            iBoutProg = iBoutProg(idxtmp(1));
          end
          assert(isscalar(iBoutProg));
          
          boutProg = allbouts(iBoutProg);
          laserEvent.intrptBehProgName = boutProg.name;
          laserEvent.intrptBehProgt0 = boutProg.t0;
        end
      end      
    end
  end
  
  function allbouts = aggregateBouts(varargin)
    % Alternative data structure for bouts
    % aggregateBouts(type1,bouts1,type2,bouts2,...)
    
    assert(mod(numel(varargin),2)==0);

    allbouts = struct('name',cell(0,1),'t0',[],'t1',[]);
    for iArg = 1:2:numel(varargin)
      type = varargin{iArg};
      bouts = varargin{iArg+1};
      for iBout = 1:numel(bouts.t0s)
        allbouts(end+1,1).name = type; %#ok<AGROW>
        allbouts(end).t0 = bouts.t0s(iBout);
        allbouts(end).t1 = bouts.t1s(iBout);
      end
    end
  end
  
  function s = hlpBout(t0s,t1s)
    s.t0s = t0s;
    s.t1s = t1s;
    s.dts = t1s-t0s;
  end
 
  function tfB = boutsOfB_RoughlyBeforeBout(boutB_t0s,boutA_t0,boutA_t1,dt_overlap)
    % tfB: logical, same size as boutB_t0s.

    assert(~any(isnan(boutB_t0s)));
    assert(~isnan(boutA_t0)&&~isnan(boutA_t1));
    
    ref = min(boutA_t0+dt_overlap,boutA_t1); 
    % if dt_overlap=0, then boutB_t0 must occur strictly before boutA_t0
    % if dt_overlap=1, then bouts B and A are allowed to occur simultaneously
    tfB = boutB_t0s < ref;    
  end
  
  function tf = hlpAtmouthBoutConstraint(atmouth_t0s,Iatmouthtf,range,minframes)
    % tf: logical, same size as atmouth_t0s. 
    
    tf = false(size(atmouth_t0s));
    for i = 1:numel(atmouth_t0s)
      intervalIdx = atmouth_t0s(i):min(atmouth_t0s(i)+range,numel(Iatmouthtf));
      tf(i) = nnz(Iatmouthtf(intervalIdx))>=minframes;
    end
  end
   
%   function tfB = boutsOfB_RoughlyAfterBoutOfA(boutB_t0s,boutB_t1s,boutA_t0,boutA_t1,dt_before,dt_after)
%     % boutA_t0, boutA_t1: [t0 t1] of bout of A
%     % tfB: logical, same size as boutB_t0s/t1s.
%     
%     assert(numel(boutB_t0s)==numel(boutB_t1s));
% 
%     tfB = false(size(boutB_t0s));
%     Nb = numel(boutB_t0s);
%     for i = 1:Nb
%       bt0 = boutB_t0s(i);
%       bt1 = boutB_t1s(i);
%       if boutA_t0<=bt0 && bt0<boutA_t1+dt_after
%         % B starts not too long after A
%         tfB(i) = true;
%       elseif boutA_t0-dt_before<=bt0 && bt0<boutA_t0 && bt1>boutA_t0
%         % B starts slightly before start of A, but extends past start of A
%         tfB(i) = true;
%       else
%         % none, tfB(i) remains false
%       end
%     end    
%   end
    
  function stats = computeDescriptiveStatistics(data,flds,setmat,setnames,varargin)
    % Compute actual descriptive stats, as opposed to the questionably-named
    % computeStats (which maybe is computeMarks)
    % 
    % data: Nexp x 1 struct array, data struct
    % flds: struct. fields: fields of data struct for which to compute 
    % stats. values: inlier limits, or predicate fcn
    % setmat: Nexp x Nset logical array. The ith col of setmat is the
    % indicator vec for membership in set i.
    % setnames: Nset x 1 cellstr. Names/labels for columns of setmat.
    % optional PVs:
    %
    % stats: Structure containing computed descriptive stats:
    % stats.fld1.set1.y: data for fld1, set1.
    % stats.fld1.set1.N: etc
    % stats.fld1.set1.mean: ...
    
    Nexp = numel(data);
    assert(iscolumn(data));
    Nset = size(setmat,2);
    assert(isstruct(flds) && all(ismember(fieldnames(flds),fieldnames(data))),'Invalid input argument ''flds''.');
    assert(islogical(setmat) && size(setmat,1)==Nexp,'Invalid input argument ''setmat''.');
    assert(iscellstr(setnames) && numel(setnames)==Nset,'Invalid input argument ''setnames''.');    

    fieldstruct = flds;
    flds = fieldnames(fieldstruct);
    Nfld = numel(flds);    
    
    stats = struct();
    for iFld = 1:Nfld
      f = flds{iFld};
      
      % compute inlier indicator vec
      inlier = fieldstruct.(f);
      if isnumeric(inlier)
        assert(isequal(size(inlier),[1 2]),'Invalid inlier limit specification.');
        %yall = getstructarrayfield(data,f,'numericscalar',true);
        yall = [data.(f)]'; % unsafe version faster
        assert(isequal(size(yall),size(data)));
        
        tfInlier = inlier(1)<= yall & yall <=inlier(2);
      elseif isa(inlier,'function_handle')
        tfInlier = arrayfun(inlier,data);
        assert(islogical(tfInlier) && isequal(size(tfInlier),size(data)));        
      else
        assert(false,'Invalid inlier specification.');
      end
      
      for iSet = 1:Nset
        tfSet = setmat(:,iSet);
        setnm = setnames{iSet};
        
%         yset = getstructarrayfield(data(tfSet),f,'numericscalar',true);
%         ysetinlier = getstructarrayfield(data(tfSet & tfInlier),f,'numericscalar',true);
        % AL: trying to optimize performance 
        yset = [data(tfSet).(f)]';
        ysetinlier = [data(tfSet & tfInlier).(f)]';
        assert(numel(yset)==nnz(tfSet));
        assert(numel(ysetinlier)==nnz(tfSet&tfInlier));
                
        tmp = struct();
        tmp.yset = yset;
        tmp.inlier = inlier;
        tmp.ysetinlier = ysetinlier;
        tmp.N = numel(ysetinlier);
        tmp.mean = mean(ysetinlier);
        tmp.median = median(ysetinlier);
        tmp.std = std(ysetinlier);
        
        stats.(f).(setnm) = tmp;
      end
    end
  end
   
  function [ypp,nfill] = fillInHoles(ypp,sz)
    % Fill in holes of size sz or smaller in postprocessed signal (ypp).
    
    assert(all(ypp>=0));
    
    % A hole is a bout in the inverted signal away from the edge.
    [bout0,bout1] = get_interval_ends(ypp==0);
    d = bout1-bout0;
    N = numel(ypp);
    tfHole = d<=sz & bout0~=1 & bout1~=N+1;
    btIdx = find(tfHole);
    nfill = numel(btIdx);
    for b = btIdx(:)'
      idx = bout0(b):bout1(b)-1;
      assert(all(ypp(idx)==0));
      ypp(idx) = 1;
    end
  end
 
  function [t,len] = firstBout(bouts)
    if isempty(bouts.t0s)
      t = nan;
      len = nan;
    else
      t = bouts.t0s(1);
      len = bouts.dts(1);
    end
  end
  
  function [t0,len] = lastBoutBeforeRef(bouts,ref)
    [t0,idx] = ExpPP.lastBoutStartBeforeRef(bouts.t0s,ref);
    if isempty(idx)
      len = nan;
    else
      len = bouts.dts(idx);
    end
  end
  
  function [t,idx] = lastBoutStartBeforeRef(t0s,ref)
    % start before if t0<ref
    
    assert(all(diff(t0s)>0),'expected monotonic');
    
    if isempty(t0s) || isnan(ref)
      % either no bouts, or no reference point
      t = nan;
      idx = [];
    else
      tf = t0s<ref;
      idx = find(tf,1,'last');
      if isempty(idx)
        % Bouts and ref exist, but no bouts before ref
        t = nan;
      else
        t = t0s(idx);
      end
    end
  end
  
  function tf = boutsBetweenRefs(t0s,ref0,ref1)
    % t0s: vector of bout starts
    % ref0/ref1: reference points, can be NaN.
    % tf: logical, same size as t0s. if tf(i) is true, t0s(i) falls in
    % [ref0,ref1), that is the ith bout starts in that interval. 
    % NOTE: If either of ref0 or ref1 is nan, no bouts will be be
    % considered "between".

    if isnan(ref0) || isnan(ref1)
      tf = false(size(t0s));
    else
      tf = ref0<=t0s & t0s<ref1;
    end
  end
    
  function n = numBoutStartsBetweenRefs(t0s,ref0,ref1)
    % starts between if t0 in [ref0,ref1)
    if isnan(ref0) || isnan(ref1)
      n = nan;
    else
      tf = ref0<=t0s & t0s<ref1;
      n = nnz(tf);
    end
  end
  
  function [i0s,i1s] = significantBouts(ypp,boutSzMin)
    [bout0,bout1] = get_interval_ends(ypp);
    sz = bout1-bout0;
    assert(all(sz>=0));
    tf = sz>=boutSzMin;    
    i0s = bout0(tf);
    i1s = bout1(tf);
  end
  
  function i = first(ypp)
    i = find(ypp,1,'first');
    if isempty(i)
      i = nan;
    end
  end
  
  function i = last(ypp)
    i = find(ypp,1,'last');
    if isempty(i)
      i = nan;
    end
  end
  
  function i = lastBeforeRef(ypp,i0)
    ypp = ypp(1:i0);
    i = ExpPP.last(ypp);
  end
  
  function Ibout = boutEndpoints2IndicatorVec(N,t0s,t1s)
    % Ibout: indicator vector of length N. Ibout(i) is 0 if frame i is not
    % located in a bout, or j if frame i is contained in bout j.
    Ibout = zeros(1,N);
    assert(numel(t0s)==numel(t1s));
    for i = 1:numel(t0s)
      idx = t0s(i):t1s(i)-1;
      Ibout(idx) = i;
    end
  end
  
  function i = firstOfBoutContainingKth(ypp,kmax)
    % KB's idea: start of bout that contains kth frame of (predicted) activity
    % This computes this stats for k = 1:kmax.

    if ~exist('kmax','var')
      kmax = 25;
    end
      
    assert(false);
    [Ibout,yppstarts] = ExpPP.boutind(ypp>0);
    frm = find(ypp,kmax,'first'); % first (up-to) kmax frames of predicted activity
    Nfrmfound = numel(frm);
    i = nan(kmax,1);
    for k = 1:kmax
      if k<=Nfrmfound
        f = frm(k); % kth frame of predicted activity
        ib = Ibout(f);
        assert(ib>0); % must be during a bout
        kbstat = yppstarts(ib);
      else % k > Nfrmfound; there is no k'th frame of predicted activity
        kbstat = nan;
      end
      i(k) = kbstat;
    end
  end
    
end

methods (Static) % misc
  
  function [tf,missingbehs] = doNamesSpanAllBasicBehaviors(names)
    % fnames: cellstr of eg behavior names, jab filenames or score filenames.
    % tf: true if all elements of ExpPP.BASICBEHAVIORS are represented (via
    % regexpmatch) in fnames.
    % missingbehs: cellstr. if tf is false, missingbehs contains missing
    % behaviors. 
    
    assert(iscellstr(names));
    missingbehs = cell(1,0);
    for beh = ExpPP.BASICBEHAVIORS(:)', beh=beh{1}; %#ok<FXSET>
      tfreg = regexpmatch(names,beh,'caseinsens',true);
      if nnz(tfreg)==0
        missingbehs{end+1} = beh; %#ok<AGROW>
      elseif nnz(tfreg)==1
        % none
      else
        warning('ExpPP:multipleFiles','More than one filename matches behavior ''%s''',beh);
      end
    end
    
    tf = isempty(missingbehs);
  end
  
  function value = loadConfigVal(rcfield)
    % Load configuration value from .rc file for ExpPP. 
    % rcfield: name of configuration field
    
    mpath = fileparts(mfilename('fullpath'));
    rcfname = fullfile(mpath,'.exppprc.mat');
    rc = struct;
    if exist(rcfname,'file'),
      try
        rc = load(rcfname);
      catch ME,
        fprintf(2,'Error loading saved configuration in %s: %s\n',rcfname,ME.getReport());
      end
    end
    if isfield(rc,rcfield)
      value = rc.(rcfield);
    else
      value = [];
    end      
  end
  
  function saveConfigVal(rcfield,val) 
    % See loadConfigVal.
    
    mpath = fileparts(mfilename('fullpath'));
    rcfname = fullfile(mpath,'.exppprc.mat');
    tmp.(rcfield) = val; %#ok<STRNU>
    if exist(rcfname,'file')
      save(rcfname,'-struct','tmp','-append');
    else
      save(rcfname,'-struct','tmp');
    end
  end
  
end
end

function tf = lclInRange(x,lims)
assert(isscalar(x)&&isnumeric(x));
tf = lims(1)<=x && x<=lims(2);
end
