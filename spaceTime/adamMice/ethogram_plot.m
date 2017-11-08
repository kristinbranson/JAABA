function [rootexpdirs,jabfiles,T0,T1] = ethogram_plot(rootexpdirs,jabfiles,nframesplot,varargin)
% ethogram_plot(rootexpdirs,jabfiles,nframesplot,varargin)
% Plots experiment ethograms.
%
% rootexpdirs: cellstr, root directories for experiments to include in 
% plot. Eg:
% rootexpdirs = {
%   '/groups/branson/home/kabram/matlab/adamTracking/videos/20130628/Laser On'
%   '/groups/branson/home/kabram/matlab/adamTracking/videos/20130628/Laser Off'
%   };
% jabfiles: cellstr, jab files for classifiers of interest. Eg:
% jabfiles = {
%   '/groups/branson/home/kabram/matlab/adamTracking/data/reach.jab'
%   '/groups/branson/home/kabram/matlab/adamTracking/data/grab.jab'
%   };
% nframesplot: Show this many frames in ethogram.
% varargin: optional P-V pairs:
%  - 'xlsexport', scalar logical.
%  - 'xlsfile', string specifying exported filename for xlsexport.
%  - 'automarks', scalar logical. If true, include automatic predictions
%     for "first lift", etc. In this case, the specified jabfiles will need 
%     to include Lift/Handopen/Grab/Sup/Atmouth/Chew.
%  - 'exptags', Either 'readfromjabs' OR (for explicit specification) cell 
%     of row cellstrs. If 'readfromjabs', ethogram_plot will look in the
%     jabfiles to find experiment tags. In the case of explicit 
%     specification, exptags must have same size as rootexpdirs; exptags{i} 
%     is a cellstr specifying tags for all experiments under rootexpdirs{i}. 
%     
% Example: ethogram_plot(expdirs,jabfiles,nframesplot,'xlsexport',true)
%
% ethogram_plot will run classifyMovie on experiments that do not have the
% appropriate scorefiles, or whose scorefiles are out of date with respect
% to their classifiers.
% 
% Double clicking: If you double click on any patch, the corresponding
% experiment and the behavior detector will be opened up in JAABA.
% Shift clicking: Shift-clicking zooms in on an experiment.

[xlsexport,xlsfile,doforce,automarks,smartorder,scoresnotjabs,dontreclassify,exptags,ethcoreargexport] =...
  myparse(varargin,...
  'xlsexport',false,...
  'xlsfile',[],...
  'doforce',false,...
  'automarks',false,...
  'smartorder',true,...
  'scoresnotjabs',false,...
  'dontreclassify',false,...
  'exptags',[],...
  'ethcoreargexport',false);

if nargin < 1,
  rootexpdirs = {};
end
if nargin < 2,
  jabfiles = {};
end
if nargin < 3
  nframesplot = [];
elseif ~isempty(nframesplot)
  assert(isnumeric(nframesplot));
  assert(isscalar(nframesplot) || numel(nframesplot)==2,'Invalid nframesplot');
end

if isempty(rootexpdirs),  
  rootexpdirs = uipickfiles('Prompt','Select root directories containing experiments','DirsOnly',true);
  if ~iscell(rootexpdirs) || isempty(rootexpdirs),
    return;
  end    
end

% ensure jabfiles
if isempty(jabfiles)
  [tfsuccess,jabfiles] = ExpPP.uiGetJabFiles;
  if ~tfsuccess
    return;
  end
end

% determine exptag status
if ischar(exptags) && strcmp(exptags,'readfromjabs')
  tfExpTagsJab = true;  
elseif ~isempty(exptags)
  % It would be odd for a user to initially specified exptags, but not
  % initially specify rootexpdirs. Nonetheless, if this occurs and the user
  % select rootexpdirs (when prompted) that correspond to their exptags,
  % then everything will work and more power to them.
  ExperimentTags.verifyTags(exptags,numel(rootexpdirs));
  tfExpTagsJab = false;
else
  exptags = ExperimentTags.expTags(rootexpdirs);
  tfExpTagsJab = false;
end

% compile expdirs, expdirtags
expdirs = {};
expdirtags = {};
for ndx = 1:numel(rootexpdirs)
  dd = findAllSubDirs(rootexpdirs{ndx});
  nn = numel(dd);
  expdirs(end+1:end+nn) = dd;
  if ~tfExpTagsJab
    expdirtags(end+1:end+nn) = exptags(ndx);
  end
end
nexps = numel(expdirs);
if nexps == 0,
  fprintf('No experiments found.\n');
  T0 = nan;
  T1 = nan;
  return;
end

% compile expdirtags, read-from-jabfile case
if tfExpTagsJab
  % compile exp->tag map
  exp2tagmap = containers.Map();
  for iJab = 1:numel(jabfiles)
    Q = loadAnonymous(jabfiles{iJab});
    jabexpdirs = Q.expDirNames;
    jabexptags = Q.expDirTags;
    if isempty(jabexptags)
      % legacy Jab, no expDirTags
      jabexptags = ExperimentTags.expTags(jabexpdirs);
    else
      jabexptags = ExperimentTags.cleanLegacyNotTags(jabexptags);
    end
    for iExp = 1:numel(jabexpdirs)
      lowerexpname = lower(jabexpdirs{iExp});
      if exp2tagmap.isKey(lowerexpname)
        warning('ethogram_plot:expTags',...
          'Experiment ''%s'' is present in more than one Jabfile. Using tags from one Jab only.',...
          lowerexpname);
      end
      exp2tagmap(lowerexpname) = jabexptags{iExp};
    end
  end
  
  % look up tags for each element of expdir
  expdirtags = ExperimentTags.expTags(expdirs);
  for iExp = 1:nexps
    lowerexpname = lower(expdirs{iExp});
    if exp2tagmap.isKey(lowerexpname);
      expdirtags{iExp} = exp2tagmap(lowerexpname);
    else
      warning('ethogram_plot:expTags',...
          'Experiment ''%s'' not present in Jabfile(s). No tags will be associated.',...
          lowerexpname);
    end
  end
end

assert(numel(expdirtags)==nexps);

% Reclassify expdirs as appropriate
if ~dontreclassify && ~scoresnotjabs
  for iExp = 1:nexps
    for iJab = 1:numel(jabfiles)
      classifyMovie(expdirs{iExp},jabfiles{iJab},'doforce',doforce,'verbose',1);
    end
  end
end

[scorefns,behaviornames,jabsforscorefiles,trxfilename] = compileJabInfo(jabfiles,scoresnotjabs);
assert(isequal(numel(scorefns),numel(behaviornames),numel(jabsforscorefiles)));

if smartorder
  % order the scorefiles/behaviors/jabnames
  [behaviornames,idx] = orderbypats(behaviornames,ExpPP.BASICBEHAVIORS);
  scorefns = scorefns(idx);
  jabsforscorefiles = jabsforscorefiles(idx);
end

% some experiments may be missing scorefiles at this point
tfGoodExp = cellfun(@(x)ExpPP.isExpDir(x,scorefns),expdirs);
if ~all(tfGoodExp)
  tmpexp = expdirs(~tfGoodExp);
  tmpstr = sprintf('%s\n',tmpexp{:});
  warningNoTrace('ethogram_plot:expsDiscarded','Excluding %d experiments with missing scorefiles: %s',...
    nnz(~tfGoodExp),tmpstr);  
end
expdirs = expdirs(tfGoodExp);
expdirtags = expdirtags(tfGoodExp);

[labels,T0,T1] = ...
  formLabelMatrix(expdirs,trxfilename,scorefns,behaviornames,nframesplot);

hfig = figure;

[tfbehspresent,missingbehs] = ExpPP.doNamesSpanAllBasicBehaviors(behaviornames);
reallydoautomarks = automarks && tfbehspresent;
if automarks && ~tfbehspresent
  % user said automarks, but basic behaviors missing
  warning('ethogram_plot:missingBehaviors',...
      'The following behaviors are missing: %s.',...
      civilizedStringFromCellArrayOfStrings(missingbehs));
end
if reallydoautomarks
  dauto = ExpPP.loadexps(expdirs,jabfiles,'includelabels',true);
  assert(numel(dauto)==numel(expdirs)); % ExpPP.loadexps does not guarantee this, but here all expdirs should be valid
  dauto = ExpPP.ensureCustomGroupInit(dauto);
  allUniqueTags = ExperimentTags.allUniqueTags(expdirtags);
  allUniqueTagsValid = matlab.lang.makeValidName(allUniqueTags);
  if numel(unique(allUniqueTagsValid))~=numel(allUniqueTagsValid)
    error('ethogram_plot:tag','Repeated tags.'); % very rare
  end
  Ntag = numel(allUniqueTags);
  for iTag = 1:Ntag
    tagRaw = allUniqueTags{iTag};
    cgrpname = allUniqueTagsValid{iTag};
    tfTag = ExperimentTags.findTag(expdirtags,tagRaw);
    dauto = ExpPP.createCustomGroup(dauto,tfTag,cgrpname);
  end
  ethplotcoreargs = {hfig,labels,T0,T1,expdirs,behaviornames,jabsforscorefiles,'doautomarks',true,'automarkdata',dauto,'scoresnotjabs',scoresnotjabs};
else
  ethplotcoreargs = {hfig,labels,T0,T1,expdirs,behaviornames,jabsforscorefiles,'scoresnotjabs',scoresnotjabs};
end
if ethcoreargexport
  fname = sprintf('ethcoreargs.%s.mat',datestr(now,'yyyymmddTHHMMSS'));
  save(fname,'ethplotcoreargs');  
end
[boutmat,line_names] = ethogram_plot_core(ethplotcoreargs{:});

if xlsexport
  if isempty(xlsfile)
    [xlsfile,xlspath] = uiputfile('bouts.txt','Export to file');
    if isequal(xlsfile,0)
      xlsfile = [];
    else
      xlsfile = fullfile(xlspath,xlsfile);
    end
  elseif exist(xlsfile,'file')
    warning('ethogram_plot:fileExists','Overwriting file ''%s''.',xlsfile);    
  end      
  
  if ~isempty(xlsfile)
    Ethogram.boutExport(xlsfile,boutmat,line_names,behaviornames);
    fprintf(1,'Exported to tab-delimited file ''%s''.\n',xlsfile);
  end
end

function [scorefns,behaviornames,jabs,trxfilename] = compileJabInfo(jabfiles,scoresnotjabs)
% Compile scorefiles/behaviors from jabfiles.
%
% jabfiles: cellstr of jabfilenames, or cellstr of scorefilenames (archaic)
% if scoresnotjabs is true
%
% scorefns: cellstr of (short) score filenames
% behaviornames: cellstr of behaviornames, one per scorefn
% jabs: cellstr of jabfiles, one per scorefn. eg, jabs{i} is the jabfile
% containing scorefns{i} and behaviornames{i}. jabs may contain repeated
% entries when multi-classifier jabs are present.
% trxfilename: char, unique trxfilename in jabs (all jabs are expected to
% be consistent)
  
if scoresnotjabs
  scorefns = jabfiles;
  behaviornames = cell(size(scorefns));
  for i = 1:numel(scorefns)
    tok = regexp(scorefns{i},'scores_(.+).mat','tokens');
    assert(~isempty(tok),'Unable to parse score filename.');
    behaviornames{i} = tok{1}{1};
  end
  jabs = scorefns; % just use scorefilenames here; no jabs avail
  trxfilename = 'trx.mat';  
else
  % jabfiles specified.
  Njab = numel(jabfiles);
  scorefns = cell(0,1);
  behaviornames = cell(0,1);
  jabs = cell(0,1);
  trxfilename = cell(Njab,1);
  for iJab = 1:Njab
    Q = loadAnonymous(jabfiles{iJab});
    
    sfns = Q.file.scorefilename;
    if ischar(sfns)
      sfns = {sfns};
    end
    behs = Labels.verifyBehaviorNames(Q.behaviors.names);
    assert(isequal(numel(sfns),numel(behs)),...
      'Inconsistent scores/behaviors in JAB file: %s',jabfiles{iJab});
    Nbeh = numel(behs);
    scorefns = cat(1,scorefns,sfns(:));
    behaviornames = cat(1,behaviornames,behs(:));
    jabs = cat(1,jabs,repmat(jabfiles(iJab),Nbeh,1));
    
    trxfilename{iJab} = Q.file.trxfilename; 
  end
  
  trxfilename = unique(trxfilename);
  assert(isscalar(trxfilename),'Jabs contain inconsistent trxfilenames.');
  trxfilename = trxfilename{1};  
end

function [labels,T0,T1] = formLabelMatrix(expdirs,trxfilename,scorefns,behaviornames,nframesplot)

% behaviornames just used to check against behaviornames in scorefiles

nexps = numel(expdirs);
nbehaviors = numel(behaviornames);
assert(numel(scorefns)==nbehaviors);
labels = false(nexps,0,nbehaviors);
tStarts = nan(1,nexps);
tEnds = nan(1,nexps);

for i = 1:nexps
  trxfile = fullfile(expdirs{i},trxfilename);
  load(trxfile,'trx');
  tStarts(i) = min([trx.firstframe]);
  tEnds(i) = max([trx.endframe]);

  % load score data  
  allScores = cell(nbehaviors,1);
  for j = 1:nbehaviors,
    scorefile = fullfile(expdirs{i},scorefns{j});
    allScores{j} = load(scorefile);
    if isfield(allScores{j},'behaviorName') && ...
       ~strcmp(allScores{j}.behaviorName,behaviornames{j})
        warning('ethogram_plot:behaviorMismatch',...
          'Behavior label mismatch in experiment %s. Scorefile: %s. Expected: %s.',...
          expdirs{i},allScores{j}.behaviorName,behaviornames{j});
    end
    allScores{j} = allScores{j}.allScores;
    if iscell(allScores{j})
      assert(isscalar(allScores{j}));
      allScores{j} = allScores{j}{1};
    end
  end
  
  for j = 1:nbehaviors
    t0 = tStarts(i);
    t1 = min(tEnds(i),numel(allScores{j}.postprocessed{1}));
    if t1 < tEnds(i),
      fprintf('Number of frames for behavior idx%d = %d < number of frames in trajectory file %s = %d\n',...
        j,t1,trxfile,tEnds(i));
      tEnds(i) = t1;
    end    
    
    if t1 > size(labels,2),
      labels = cat(2,labels,nan(nexps,t1-size(labels,2),nbehaviors));
    end
    
    labels(i,t0:t1,j) = allScores{j}.postprocessed{1}(t0:t1);
  end
end

T0 = min(tStarts);
T1 = max(tEnds)-1;

if ~isempty(nframesplot)
  if numel(nframesplot) == 1,
    T1 = min(T1,T0+nframesplot-1);
  else
    assert(numel(nframesplot)==2);
    T0 = nframesplot(1);
    T1 = min(T1,nframesplot(2));
  end
end
%nframesplot = T1-T0+1;

labels = labels(:,T0:T1,:);


function [cellstr,idx] = orderbypats(cellstr,pats)
% sort cellstr by pats (patterns) as best as possible

Npat = numel(pats);
loc = nan(size(cellstr));
for i = 1:Npat
  tf = regexpmatch(cellstr,pats{i},'caseinsens',true);
  loc(tf) = i;
end
loc(isnan(loc)) = inf; % elements of cellstr that didn't match go at the end

% loc(i) gives the desired relative position of cellstr{i}

[~,idx] = sort(loc);
cellstr = cellstr(idx);
