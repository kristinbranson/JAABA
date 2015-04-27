classdef Jab < handle
  % .jab files
  
  methods (Static)
    
    function J = load(jabfile)
      % J = load(jabfile)
      %
      % J: Modernized Macguffin object

      assert(ischar(jabfile) && exist(jabfile,'file')==2,...
        'Cannot find jabfile ''%s''.',jabfile);

      J = loadAnonymous(jabfile);
      if isstruct(J)
        J = Macguffin(J);
      end
      J.modernize();
    end
    
    function save(J,jabfile)
      saveAnonymous(jabfile,J);
    end
    
    function savePrompt(J,jabfile)
      if exist('jabfile','var')==0
        jabfile = [];
      end
      [fname,pname] = uiputfile('*.jab','Save',jabfile);
      fname = fullfile(pname,fname);
      Jab.save(J,fname);
    end
    
    function clearLabels(jabfile,realbeh)
      % Clear all labels for realbeh and No-realbeh       
      
      J = Jab.load(jabfile);
      
      jabrealbehs = J.behaviors.names(1:J.behaviors.nbeh);
      tf = strcmpi(realbeh,jabrealbehs);
      assert(any(tf),'Specified behavior ''%s'' not present in jabfile.',realbeh);
      assert(nnz(tf)==1);
      realbeh = jabrealbehs{tf}; % case could be different
      
      tfMultiCls = J.behaviors.nbeh>1;
      if tfMultiCls
        nobeh = Labels.noBehaviorName(realbeh);
      else
        nobeh = 'None';
      end
      
      J.labels = Labels.clearLabels(J.labels,realbeh,realbeh);
      J.labels = Labels.clearLabels(J.labels,nobeh,realbeh);
      
      J.savePrompt(J,jabfile);
    end
    
    function tf = isMultiClassifier(jabfile)
      % tf = isMultiClassifier(jabfile)
      J = Jab.load(jabfile);
      tf = J.isMultiClassifier();
    end
    
    function rmClassifier(jabfile,clsnames)
      % rmClassifier(jabfile,clsnames)
      %
      % jab: char, jab filename
      % clsname: char or cellstr
      
      J = Jab.load(jabfile);
      Jbeh = J.behaviors;
      Jclsnames = Jbeh.names(1:Jbeh.nbeh);
      
      if exist('clsnames','var')>0
        if ischar(clsnames)
          clsnames = cellstr(clsnames);
        end
        assert(iscellstr(clsnames),'Expected ''clsname'' argument to be a char or cellstr.');                
      else
        [sel,ok] = listdlg('PromptString','Select classifiers to remove:',...
          'ListString',Jclsnames,...
          'ListSize',[240 240]);
        if ~ok
          return;
        end
        clsnames = Jclsnames(sel);
      end
      
      [tf,iClsRm] = ismember(clsnames,Jclsnames);
      if ~all(tf)
        clsmissing = clsnames(~tf);
        error('Jab:rmClassifier','Classifier(s) %s not present in jabfile.',...
          civilizedStringFromCellArrayOfStrings(clsmissing));
      end
      
      nClsRm = numel(iClsRm);
      nClsNew = Jbeh.nbeh - nClsRm;
      if nClsNew==0
        error('Jab:rmClassifier','Cannot remove all classifiers from a .jab file.');
      end
      
      iLblsRm = [iClsRm iClsRm+Jbeh.nbeh];
      origBehs = J.behaviors.names(iClsRm);
      origNoBehs = J.behaviors.names(iClsRm+Jbeh.nbeh);
      J.behaviors.names(iLblsRm) = [];
      J.behaviors.labelcolors(iLblsRm,:) = [];
      J.behaviors.nbeh = nClsNew;
      J.file.scorefilename(iClsRm) = [];
      J.windowFeaturesParams(iClsRm) = [];
      for iCls = 1:nClsRm      
        J.labels = Labels.removeClassifier(J.labels,origBehs{iCls},origNoBehs{iCls});
        J.gtLabels = Labels.removeClassifier(J.gtLabels,origBehs{iCls},origNoBehs{iCls});
      end      
      J.classifierStuff(iClsRm) = [];
      
      % Quirk, no-beh name for single-classifier projs
      if nClsNew==1
        oldNoBehavior = J.behaviors.names{2};
        newNoBehavior = 'None';
        J.behaviors.names{2} = newNoBehavior;
        
        J.labels = Labels.renameBehaviorRaw(J.labels,oldNoBehavior,newNoBehavior);
        J.gtLabels = Labels.renameBehaviorRaw(J.gtLabels,oldNoBehavior,newNoBehavior);
        
        % ALXXX to be fixed
        warning('Jab:rmClassifier',...
          'Going from multi-classifier to single-classifier project. Label importance may not be set correctly for GT mode.');
      end
      
      Jab.savePrompt(J,jabfile);
    end
     
    function tf = formatChangedBetweenVersions(v1,v2)
      % v2 should be later-than-or-equal-to-v1.
      %
      % This method introduced in version 0.6.0. Jabs saved in prior
      % versions have version prop enclosed in double-cell array. We avoid
      % the complications for now and just return true if v2 is 0.6.0 and
      % v1 is not. This will need to be updated with subsequent versions.

      tf = strcmp(v2,'0.6.0') && ~strcmp(v1,'0.6.0');
    end  
    
    function merge(jabfiles,jabout)
      % jabMerge(jabfiles,jabout)
      % jabfiles: optional. cellstr of jab filenames
      % jabout: optional. output jab filename

      % ALTODO: verify behavior, eg this uses ExpPP.loadConfigVal

      
      if ~exist('jabfiles','var') || isempty(jabfiles)
        [tfsuccess,jabfiles] = ExpPP.uiGetJabFiles('promptstr','Select jab files to combine');
        if ~tfsuccess
          return;
        end
      end
      
      if ischar(jabfiles)
        jabfiles = cellstr(jabfiles);
      end
      assert(iscellstr(jabfiles),'Expected ''jabfiles'' to be a cellstr of jab filenames.');
      
      % ALTODO: create a new verify method
      %Macguffin.jabVerifyAug2014(jabfiles);      
      
      Q = cellfun(@loadAnonymous,jabfiles,'uni',0);
      Q = cat(1,Q{:});
      Qmerge = Macguffin(Q);
      
      % come up with a proposed name for combined jab
      if ~exist('jabout','var') || isempty(jabout)
        MAXFILENAMELENGTH = 45;
        behnames = Labels.verifyBehaviorNames(Qmerge.behaviors.names);
        combjabname = '';
        for i = 1:numel(behnames)
          combjabname = [combjabname behnames{i} '.']; %#ok<AGROW>
          if numel(combjabname)>MAXFILENAMELENGTH
            combjabname = [combjabname 'etc.']; %#ok<AGROW>
            break;
          end
        end
        combjabname = [combjabname 'jab'];

        jabpath = ExpPP.loadConfigVal('jabpath');
        if isempty(jabpath)
          jabpath = pwd;
        end    
        [filename,pathname] = ...
          uiputfile({'*.jab','JAABA files (*.jab)'},'Save combined jabfile',fullfile(jabpath,combjabname));
        if ~ischar(filename),
          % user hit cancel
          return;
        end
        
        jabout = fullfile(pathname,filename);
      else
        if exist(jabout,'file')
          uiwait(warndlg(sprintf('Output file ''%s'' exists and will be overwritten.',jabout)));
        end
      end
        
      tmp.x = Qmerge; %#ok<STRNU>
      save(jabout,'-struct','tmp');
      ExpPP.saveConfigVal('jabpath',fileparts(jabout));
    end    
    
    function [scorefiles,behaviors] = jabfileScoresBehs(jabname)
      % [scorefiles,behaviors] = jabfileScoresBehs(jabname)
      % scorefile: cellstr of score files
      % behaviors: cellstr of "real" behaviors (doesn't include No_<beh>)
      
      x = loadAnonymous(jabname);
      scorefiles = x.file.scorefilename;      
      if ischar(scorefiles)
        scorefiles = cellstr(scorefiles);
      end
      behaviors = Labels.verifyBehaviorNames(x.behaviors.names);      
      assert(numel(scorefiles)==numel(behaviors));
    end    
    
    
    
  end
  
end
  
  
  
  