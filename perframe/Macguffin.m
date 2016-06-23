classdef Macguffin < handle
  % This holds all the stuff that gets saved to the .jab file.  In fact,
  % the .jab file is a .mat file holding one variable, x, which is a
  % Macguffin.
  
  properties
    featureLexiconName
    featureLexicon
    behaviors
    file
    trxGraphicParams
    scoreFeatures
    sublexiconPFNames
    windowFeaturesParams
    labels
    gtLabels
    expDirNames
    gtExpDirNames
    classifierStuff
    gtSuggestions
    extra=struct()  % a structure that stores additional information
    version=''
      % 0.5.0 : Original
      % 0.5.1 : Supports nextra_markers, flies_extra_markersize, 
      %           flies_extra_marker, and flies_extra_linestyle fields in
      %           trxGraphicParams.
      % 0.6.0 : Multiclassifier and ST updates.
  end
    
  % -----------------------------------------------------------------------
  methods (Access=private)
    % ---------------------------------------------------------------------
    function initFromFeatureLexiconName(self,featureLexiconName)
      [featureLexicon,animalType]=featureLexiconFromFeatureLexiconName(featureLexiconName);
      self.featureLexiconName=featureLexiconName;
      self.behaviors.type=animalType;
      self.behaviors.names = {};  % no behaviors defined yet
      self.file.perframedir = 'perframe';
      self.file.clipsdir = 'clips';
      %self.windowfeatures = struct;
      self.behaviors.labelcolors = [JLabelGUIData.COLOR_BEH_DEFAULT JLabelGUIData.COLOR_NOBEH_DEFAULT];
      self.behaviors.unknowncolor = [0,0,0];
      self.trxGraphicParams=trxGraphicParamsFromAnimalType(animalType);
      self.file.scorefilename = '';
      self.file.trxfilename = '';
      self.file.moviefilename = '';
      self.file.movieindexfilename = 0;
      self.scoreFeatures = Macguffin.ScoreFeatures;
      self.featureLexicon=featureLexicon;
      featureLexiconPFNames = fieldnames(featureLexicon.perframe);
      self.sublexiconPFNames = featureLexiconPFNames;
      self.windowFeaturesParams={struct()};  % scalar struct with no fields
        % This is valid b/c no per-frame features have been enabled yet
      self.labels = Labels.labels(0);
      self.gtLabels = Labels.labels(0);
      self.expDirNames={};
      self.gtExpDirNames={};
      self.classifierStuff=ClassifierStuff();
      
      self.extra.perframe.params.fov = pi;
      self.extra.perframe.params.thetafil = [0.0625,0.2500,0.3750,0.2500,0.0625];
      self.extra.perframe.params.nbodylengths_near = 2.5;
      self.extra.perframe.params.max_dnose2ell_anglerange = 127;
      self.extra.perframe.params.nroi = 0;
      self.extra.perframe.params.nflies_close = [];
      
      self.extra.perframe.landmarkParams.arena_center_mm_x = 0;
      self.extra.perframe.landmarkParams.arena_center_mm_y = 0;
      self.extra.perframe.landmarkParams.arena_radius_mm = 60;
      self.extra.perframe.landmarkParams.arena_type = 'circle';

      self.extra.usePastOnly = false;

      self.addversion
      
    end  % method
    
    
    % ---------------------------------------------------------------------
    function initFromJLabelData(self,jld)
      self.featureLexiconName=jld.featureLexiconName;
      self.featureLexicon=jld.featureLexicon;
      self.scoreFeatures=jld.scoreFeatures;
      subdialectPFNames=jld.allperframefns;
      nScoreFeaturess=length(jld.scoreFeatures);
      sublexiconPFNames=subdialectPFNames(1:end-nScoreFeaturess); %#ok<PROP>
      self.sublexiconPFNames=sublexiconPFNames; %#ok<PROP>
      self.behaviors.type=jld.targettype;
      self.behaviors.names=jld.labelnames;
      self.behaviors.labelcolors=jld.labelcolors;
      self.behaviors.unknowncolor=jld.unknowncolor;
      % TODO AL: .behaviors.names seems more like label names. Resolve
      % timelines vs behaviors naming
      self.behaviors.nbeh = Labels.determineNumTimelines(jld.labelnames); 
      self.file.moviefilename=jld.moviefilename;
      self.file.movieindexfilename=jld.movieindexfilename;
      self.file.trxfilename=jld.trxfilename;
      self.file.scorefilename=jld.scorefilename;
      self.file.clipsdir=jld.clipsdir;                  
      self.file.perframedir=jld.perframedir;
      self.file.stfeatures=jld.stfeatures;
      self.trxGraphicParams=jld.trxGraphicParams;
      if isnonempty(jld.landmark_params) ,
        self.extra.perframe.landmarkParams=jld.landmark_params;
      end
      if isnonempty(jld.perframe_params) ,
        self.extra.perframe.params=jld.perframe_params;
      end
      % Get the labels, put them in self
      if jld.gtMode ,
        self.gtExpDirNames=jld.expdirs;
        self.expDirNames=jld.otherModeLabelsEtc.expDirNames;
      else
        self.expDirNames=jld.expdirs;
        self.gtExpDirNames=jld.otherModeLabelsEtc.expDirNames;
      end
      [self.labels,self.gtLabels]=jld.getLabelsAndGTLabels();
      % Get the window feature params, put in self
      self.windowFeaturesParams=jld.windowfeaturesparams;
      % Put the classifier in self
      self.classifierStuff=jld.getClassifierStuff();
      self.version = jld.version;

      if isprop(jld,'usePastOnly'),
        self.extra.usePastOnly = jld.usePastOnly;
      end
      
      % save GTSuggestions by default
      self.gtSuggestions.GTSuggestionMode = jld.GTSuggestionMode;
      self.gtSuggestions.randomGTSuggestions = jld.randomGTSuggestions;
      self.gtSuggestions.thresholdGTSuggestions = jld.thresholdGTSuggestions ;
      self.gtSuggestions.loadedGTSuggestions = jld.loadedGTSuggestions;
      self.gtSuggestions.balancedGTSuggestions = jld.balancedGTSuggestions ;


    end  % method

    
    % ---------------------------------------------------------------------
    function initFromBasicDataStruct(self,basicDataStruct)
      % Init from the featureLexiconName to start
      self.initFromFeatureLexiconName(basicDataStruct.featureLexiconName);
      % Copy over specific fields to start
      if isfield(basicDataStruct,'extra')
        self.extra=basicDataStruct.extra;
      end
      % Now overwrite other properties defined in the struct
      fieldNames=fieldnames(basicDataStruct);
      for iField=1:length(fieldNames)
        fieldName=fieldNames{iField};
        if ismember(fieldName,{'featureLexiconName' 'extra'})
          % already copied over
          continue
        elseif isprop(self,fieldName)
          self.(fieldName)=basicDataStruct.(fieldName);
        else
          % non-property fields get put in the extra property
          self.extra.(fieldName)=basicDataStruct.(fieldName);
        end
      end  % for loop
      self.addversion
    end  % method

    
    % ---------------------------------------------------------------------
    function initFromOldStyleProjectAndClassifier(self, ...
                                                  projectParams, ...
                                                  classifierParams, ...
                                                  gtExpDirNames, ...
                                                  gtLabels, ...
                                                  scoreFeatureMatFileNames, ...
                                                  scoreFeatureJabFileNames)
    
      % process args
      if ~exist('scoreFeatureMatFileNames','var') || ~exist('scoreFeatureJabFileNames','var')
        scoreFeatureMatFileNames=cell(0,1);
        scoreFeatureJabFileNames=cell(0,1);
      end
      
      % get the featureLexicon from the relevant file
      featureLexiconFileNameRel=projectParams.file.featureconfigfile;
      pathToMisc=fileparts(mfilename('fullpath'));
      pathToJaaba=fileparts(pathToMisc);
      featureLexiconFileNameAbs= ...
        fullfile(pathToJaaba,'perframe',featureLexiconFileNameRel);
      featureLexicon = ReadXMLParams(featureLexiconFileNameAbs);

      % Try to match up the relative feature lexicon file name with a
      % feature lexicon name
      [featureLexiconNameList, ...
       featureLexiconFileNameRelList] = ...
        getFeatureLexiconListsFromXML();
      isSameLexicon=strcmp(featureLexiconFileNameRel,featureLexiconFileNameRelList);
      i=find(isSameLexicon);  
      if isempty(i) ,
        featureLexiconName='custom';
        %featureLexiconAnimalType='';
      else
        featureLexiconName=featureLexiconNameList{i};
        %featureLexiconAnimalType=featureLexiconAnimalTypeList{i};
      end

      % Convert featureparamlist, if present, to new format
      % If not present, the sublexicon is identical to the lexicon
      if isfield(projectParams,'featureparamlist') ,
        sublexiconPFNames=fieldnames(projectParams.featureparamlist);
      else
        sublexiconPFNames=fieldnames(featureLexicon.perframe);
      end

      % copy the project params over
      self.behaviors=projectParams.behaviors;
      fileParams=projectParams.file;
      % delete some fields that we don't need/want
      if isfield(fileParams,'featureconfigfile') ,
        fileParams=rmfield(fileParams,'featureconfigfile');
      end
      self.file=fileParams;
      if isfield(projectParams,'trx')
        trxGraphicParamsRaw=projectParams.trx;
      elseif isfield(projectParams,'plot')  && ...
             isfield(projectParams.plot,'trx')
        trxGraphicParamsRaw=projectParams.plot.trx;
      else
        trxGraphicParamsRaw=trxGraphicParamsFromAnimalType(animalType);
      end
      self.trxGraphicParams=cookTrxGraphicParams(trxGraphicParamsRaw);
      self.featureLexiconName=featureLexiconName;
      self.featureLexicon=featureLexicon;
      self.sublexiconPFNames=sublexiconPFNames;
      self.scoreFeatures= ...
        replaceScoreFeatureClassifierFileNames(projectParams.scoresinput, ...
                                               scoreFeatureMatFileNames, ...
                                               scoreFeatureJabFileNames);
%       if isempty(projectParams.behaviors.names) ,
%         self.behaviorName='';
%       else
%         self.behaviorName=projectParams.behaviors.names{1};
%       end
      if isfield(projectParams,'perframe') && isfield(projectParams.perframe,'landmark_params')
        self.extra.perframe.landmarkParams=projectParams.perframe.landmark_params;
      else
        self.extra.perframe.landmarkParams=[];
      end
      
      if isfield(projectParams,'perframe') && isfield(projectParams.perframe,'params')
        self.extra.perframe.params=projectParams.perframe.params;
      else
        self.extra.perframe.params=[];
      end

      % copy the experiment dirs, labels over
      if isempty(classifierParams) ,
        % if no classifierParams provided
        self.appendEmptyLabelsAndDefaultClassifier(projectParams);
      else
        % The usual case: caller provided a non-empty classifierParams
        self.appendClassifierAndLabels(projectParams,classifierParams);
      end

      % Make sure the GT experiment dir names are absolute paths
      nGTExpDirs=length(gtExpDirNames);
      gtExpDirAbsPathNames=cell(nGTExpDirs,1);
      for i=1:nGTExpDirs
        gtExpDirName=gtExpDirNames{i};
        if isFileNameAbsolute(gtExpDirName) ,
          gtExpDirAbsPathNames{i}=gtExpDirName;
        else
          gtExpDirAbsPathNames{i}=fullfile(pwd(),gtExpDirName);
        end
      end

      % append the GT labels
      self.gtExpDirNames=gtExpDirAbsPathNames;
      self.gtLabels=cookLabels(gtLabels);
      self.addversion;
    end  % method
    
    function initFromMacguffins(self,m)
      % init from 1+ Macguffins
      % Combines multiple Macguffins into single multi-classifier Macguffin
                    
      assert(isa(m,'Macguffin'));
      assert(~isempty(m),'Input Macguffin cannot be empty,');
      
      m.modernize();
      
      flds = fieldnames(m);
            
      % Special case early return 
      if isscalar(m)
        for f = flds(:)',f=f{1}; %#ok<FXSET>
          self.(f) = m.(f);
        end
        return;
      end
      
      %%% Check all fields that must be the same for all input jabs 
      % (otherwise we don't know how to merge)
      FIELDS_ALLOWED_TO_DIFFER = {'behaviors' 'file' 'trxGraphicParams' ...
        'windowFeaturesParams' 'labels' 'expDirNames' 'classifierStuff' 'version'};
      fieldsMustBeSame = setdiff(flds,FIELDS_ALLOWED_TO_DIFFER);
      for f = fieldsMustBeSame(:)', f=f{1}; %#ok<FXSET>
        vals = {m.(f)};
        if ~isequaln(vals{:})
          error('Macguffin:merge',...
            'Currently unsupported merge: input objects differ in field ''%s''.',f);
        end
        self.(f) = vals{1};
      end      

      %%% Behaviors
      Nobj = numel(m);
      allbehs = cell(1,0);
      allnobehs = cell(1,0);
      nbehs = nan(Nobj,1);
      allsfns = cell(1,0);
      allwfp = cell(1,0);
      allcs = ClassifierStuff.empty(0,1);
      for i = 1:Nobj
        bhvrs = m(i).behaviors;
        if ~strcmp(bhvrs.type,m(1).behaviors.type)
          error('Macguffin:incompatibleObjs',...
            'Unsupported merge: input objects differ in field ''%s''.','behaviors');
        end
        [behs,nobehs] = Labels.verifyBehaviorNames(bhvrs.names);
        
        allbehs = [allbehs behs]; %#ok<AGROW>
        allnobehs = [allnobehs nobehs]; %#ok<AGROW>
        nbehs(i) = numel(behs); % record for comparison to classifierStuff
        
        sfn = m(i).file.scorefilename;
        assert(numel(sfn)==nbehs(i));
        allsfns = [allsfns sfn]; %#ok<AGROW>
        
        wfp = m(i).windowFeaturesParams;
        assert(numel(wfp)==nbehs(i));
        allwfp = [allwfp wfp]; %#ok<AGROW>
        
        cs = m(i).classifierStuff(:);
        assert(numel(cs)==nbehs(i));
        allcs = [allcs; cs(:)]; %#ok<AGROW>
      end
      % Print some diagnostics
      allbehsUn = unique(allbehs);
      allbehsUnCnt = cellfun(@(x)nnz(strcmp(x,allbehs)),allbehsUn);
      nAllbehsUn = numel(allbehsUn);
      fprintf(1,'New jab contents:\n');      
      for iBeh = 1:numel(allbehsUn)
        fprintf(1,'  Classifier %s (%d objs)\n',allbehsUn{iBeh},allbehsUnCnt(iBeh));
      end
      % Set combined .behaviors 
      self.behaviors.type = m(1).behaviors.type;
      self.behaviors.names = Labels.behnames2labelnames(allbehsUn);
      % Note, labelcolors are reset and not carried over
      self.behaviors.labelcolors = Labels.cropOrAugmentLabelColors(zeros(1,0),numel(self.behaviors.names),'darkened');
      %self.behaviors.labelcolors = Labels.addNoBehColors(self.behaviors.labelcolors);
      self.behaviors.labelcolors = reshape(self.behaviors.labelcolors,3,[])';
      self.behaviors.unknowncolor = m(1).behaviors.unknowncolor;
      self.behaviors.nbeh = nAllbehsUn;
      
      %%% Files
      % all fields but scorefilename must be equal
      files = {m.file};
      filesNoSfn = cellfun(@(x)rmfield(x,'scorefilename'),files,'uni',0);
      if ~isequaln(filesNoSfn{:})
        error('Macguffin:incompatibleObjs',...
          'Currently unsupported merge: input objects differ in field ''.file''.');
      end
      self.file = files{1}; % self.file.scorefilename set below
      % Check for inconsistent scorefilenames for repeated behaviors
      sfn = cell(nAllbehsUn,1);
      for iBeh = 1:nAllbehsUn
        beh = allbehsUn{iBeh};
        defsfn = ScoreFile.defaultScoreFilename(beh);
        idx = strcmp(beh,allbehs);
        if ~all(strcmp(defsfn,allsfns(idx)))
          warningNoTrace('Macguffin:merge','Scorefile name for classifier ''%s'' reset to ''%s''.',...
            beh,defsfn);
        end
        sfn{iBeh} = defsfn;
      end
      self.file.scorefilename = sfn(:)';
      
      
      %%% trxGraphicParams
      vals = {m.trxGraphicParams};
      if isequaln(vals{:})
        self.trxGraphicParams = vals{1};
      else
        warningNoTrace('Macguffin:merge',...
          'Projects differ in field .trxGraphicParams. This property will be re-initialized from the animal type.');
        self.trxGraphicParams = trxGraphicParamsFromAnimalType(self.behaviors.type);
      end
            
      %%% expdirnames
      % - compile combined expdir list
      arrayfun(@(x)assert(isrow(x.expDirNames),...
        'Expected .expDirNames to be a rowvector.'),m);
      allexpdirnames = cat(2,m.expDirNames);
      allexpdirnamesUn = unique(allexpdirnames);
      allexpdirnamesUnCnt = cellfun(@(x)nnz(strcmp(x,allexpdirnames)),allexpdirnamesUn);
      nAllExp = numel(allexpdirnamesUn);
      nRptExp = nnz(allexpdirnamesUnCnt>1);
      fprintf(1,'Merged object will contain %d experiments in all.\n',nAllExp);
      fprintf(1,'%d/%d experiments are present in more than one input object.\n',nRptExp,nAllExp);      
      self.expDirNames = allexpdirnamesUn;
      
      %%% Labels
      alllabels = cell(Nobj,1); % contains modified .labels struct for each MacGuffin
      for i = 1:Nobj
        lbls = m(i).labels;
        if nAllbehsUn>1 
          [behs,nobehs] = Labels.verifyBehaviorNames(m(i).behaviors.names);
          if numel(behs)==1 && strcmpi(nobehs,'none')
            % We are now a multiclassifier jab, but input obj i is a single
            % classifier jab. Convert 'None' to No_<beh> in labels
            lbls = Labels.renameBehaviorRaw(lbls,nobehs{1},Labels.noBehaviorName(behs{1}));
          end
        end
        assert(numel(lbls)==numel(m(i).expDirNames));
        alllabels{i} = lbls(:)';
      end
      alllabels = cat(2,alllabels{:});
      assert(numel(alllabels)==numel(allexpdirnames));
      % compile/combine labels
      self.labels = Labels.compileLabels(allexpdirnamesUn,alllabels,allexpdirnames,self.behaviors.names);
        
      %%% windowFeaturesParams
      assert(numel(allbehs)==numel(allwfp));
      wfpNew = cell(1,nAllbehsUn);
      for iBeh = 1:nAllbehsUn
        beh = allbehsUn{iBeh};
        idx = strcmp(beh,allbehs);
        wfpBeh = allwfp(idx); % all wfps for this behavior
        if numel(wfpBeh)>1 && ~all(isequaln(wfpBeh{:}))
          error('Macguffin:incompatibleObjs',...
            'Unsupported merge: Differing .windowFeaturesParams encountered for behavior ''%s''.',beh);
        end
        wfpNew{iBeh} = wfpBeh{1};
      end
      self.windowFeaturesParams = wfpNew;
    
      %%% classifierStuff
      % Consider moving into ClassifierStuff
      assert(numel(allbehs)==numel(allcs));
      csNewAll = [];
      for iBeh = 1:nAllbehsUn
        beh = allbehsUn{iBeh};
        idx = strcmp(beh,allbehs);
        csBeh = allcs(idx); % all classifierstuffs for this behavior        
        csNew = ClassifierStuff();

        CSFIELDS_SAME = {'type' 'postProcessParams' 'trainingParams' 'featureNames'};
        for f = CSFIELDS_SAME,f=f{1}; %#ok<FXSET>
          val = {csBeh.(f)};
          if numel(csBeh)>1 && ~isequaln(val{:})
            error('Macguffin:incompatibleObjs',...
              'Unsupported merge: ClassifierStuffs for behavior ''%s'' differ in field ''%s''.',...
              beh,f);
          end
          csNew.(f) = val{1};
        end
        
        CSFIELDS_CLEAR = {'params' 'timeStamp' 'windowdata' 'savewindowdata','selFeatures','predictOnlyCurrentFly'};
        for f = CSFIELDS_CLEAR,f=f{1}; %#ok<FXSET>
          val = {csBeh.(f)};
          val{1,end+1} = csNew.(f); %#ok<AGROW> csNew contains default/new value of field f
          if ~isequaln(val{:})
            warningNoTrace('Macguffin:clearField',...
              'Field ''%s'' in classifierStuff for behavior ''%s'' will be reset.',...
              f,beh);
          end
        end
        
        CSFIELD_PICKONE = {'confThresholds' 'scoreNorm'};        
        for f = CSFIELD_PICKONE,f=f{1}; %#ok<FXSET>
          val = {csBeh.(f)};
          if numel(val)>1 && ~isequaln(val{:})
            warningNoTrace('Macguffin:oneField',...
              'Differing values present for classifierStuff field ''%s'' for behavior ''%s''. One value will be selected.',...
              f,beh);
          end
          csNew.(f) = val{1};
        end
        
        csNewAll = [csNewAll;csNew]; %#ok<AGROW>
      end
      self.classifierStuff = csNewAll;
            
      self.addversion;
    end
            
    
    % ---------------------------------------------------------------------
    function appendClassifierAndLabels(self,projectParams,classifierParams)
      % this is used when converting old-style files to .jab files
      self.expDirNames=classifierParams.expdirs;
      if isfield(classifierParams,'labels')
        self.labels=cookLabels(classifierParams.labels);
      else
        % Look for labels in the experiment directories
        self.labels=getLabelsFromExpDirs(self.expDirNames,projectParams.file.labelfilename);
      end
      if isfield(classifierParams,'windowfeaturesparams')
        self.windowFeaturesParams=classifierParams.windowfeaturesparams;
      elseif isfield(projectParams.windowfeatures,'windowfeaturesparams')
        self.windowFeaturesParams=projectParams.windowfeatures.windowfeaturesparams;
      end
      % if isfield(classifierParams,'gt_labels')
      %   gtLabels=classifierParams.gt_labels;
      % else
      %   nExps=length(classifierParams.expdirs);
      %   gtLabels=struct();
      %   for i=1:nExps
      %     gtLabels(i).t0s={};
      %     gtLabels(i).t1s={};
      %     gtLabels(i).names={};
      %     gtLabels(i).flies=[];
      %     gtLabels(i).off=[];
      %     gtLabels(i).timestamp={};
      %     gtLabels(i).imp_t0s={};
      %     gtLabels(i).imp_t1s={};
      %   end  
      % end
      % everythingParams.gtLabels=gtLabels;

      % copy the classifier params proper over
      classifierStuff=ClassifierStuff();
      %classifierStuff.animalType=projectParams.behaviors.type;
      %if ~isempty(projectParams.behaviors.names)                             
      %  classifierStuff.behaviorName=projectParams.behaviors.names{1};
      %else
      %  classifierStuff.behaviorName='';
      %end
      classifierStuff.type=classifierParams.classifiertype;
      classifierStuff.params=classifierParams.classifier;
      classifierStuff.confThresholds=classifierParams.confThresholds;
      classifierStuff.scoreNorm=classifierParams.scoreNorm;
      classifierStuff.postProcessParams=classifierParams.postprocessparams;
      classifierStuff.trainingParams=classifierParams.classifier_params;
      classifierStuff.timeStamp=classifierParams.classifierTS;
      self.classifierStuff=classifierStuff;                              
    end  % method
    
    
    % ---------------------------------------------------------------------
    function appendEmptyLabelsAndDefaultClassifier(self,projectParams)
      % this is used when converting old-style files to .jab files
      self.expDirNames={};
      self.labels=struct([]);
      if isfield(projectParams.windowfeatures,'windowfeaturesparams')
        self.windowFeaturesParams=projectParams.windowfeatures.windowfeaturesparams;
      end
      classifierStuff=ClassifierStuff();
      %animalType=projectParams.behaviors.type;
%       classifierStuff.type='boosting';  % e.g., 'boosting'
%       classifierStuff.params=struct([]);
%       classifierStuff.confThresholds=[0 0];
%       classifierStuff.scoreNorm=[];
%       classifierStuff.postProcessParams=struct([]);
%       classifierStuff.trainingParams= ...
%         struct('iter',100, ...
%                'iter_updates',10, ...
%                'numSample',2500, ...
%                'numBins',30, ...
%                'CVfolds',7, ...
%                'baseClassifierTypes','Decision Stumps', ...
%                'baseClassifierSelected',1);
%       classifierStuff.timeStamp=[];
      self.classifierStuff=classifierStuff;                            
    end  % method   
  end  % private methods
  
  
  % -----------------------------------------------------------------------
  methods
    % ---------------------------------------------------------------------
    function self=Macguffin(varargin)
      if length(varargin)==1 && isa(varargin{1},'Macguffin')
        self.initFromMacguffins(varargin{1});
      elseif length(varargin)==1 && ischar(varargin{1})
        self.initFromFeatureLexiconName(varargin{1});
      elseif length(varargin)==1 && isequal(class(varargin{1}),'JLabelData')
        jld=varargin{1};
        self.initFromJLabelData(jld);
      elseif length(varargin)==1 && isstruct(varargin{1})
        basicDataStruct=varargin{1};
        self.initFromBasicDataStruct(basicDataStruct);
      elseif length(varargin)==4 || length(varargin)==6
        projectParams=varargin{1};
        classifierParams=varargin{2};
        gtExpDirNames=varargin{3};
        gtLabels=varargin{4};
        if length(varargin)==6
           scoreFeatureMatFileNames=varargin{5};
           scoreFeatureJabFileNames=varargin{6};
        else
           scoreFeatureMatFileNames=cell(0,1);
           scoreFeatureJabFileNames=cell(0,1);
        end
        self.initFromOldStyleProjectAndClassifier(projectParams, ...
                                                  classifierParams, ...
                                                  gtExpDirNames, ...
                                                  gtLabels, ...
                                                  scoreFeatureMatFileNames, ...
                                                  scoreFeatureJabFileNames);
      else
        error('Macguffin:badArgumentsToConstructor', ...
              'The arguments to Macguffin() are no good');
      end
    end  % constructor method

    % ---------------------------------------------------------------------
    function tf=isMultiClassifier(self)
      % tf=isMultiClassifier(self)
      tf = isfield(self.behaviors,'nbeh') && self.behaviors.nbeh>1;
    end
    
    % ---------------------------------------------------------------------
    function result=getMainBehaviorName(self)
      % The name of the "main" behavior, if present.  The "main" behavior
      % is the first one that is not "None" or "none" or some such.
      % This is not a dependent property because other code sometimes turns
      % a Macguffin into a struct when no behaviors are defined, and doing
      % that calls all the get. methods for the dependent parameters, which
      % throws an error.
      mainBehaviorIsDefined=false;
      if isfield(self.behaviors,'names')
        behaviorNames=self.behaviors.names;
        isNone=strcmpi('none',behaviorNames);
        realBehaviorNames=behaviorNames(~isNone);
        if isnonempty(realBehaviorNames)
          mainBehaviorIsDefined=true;
          result=realBehaviorNames{1};
        end
      end
      if ~mainBehaviorIsDefined
        error('Macguffin:mainBehaviorNotDefined', ...
              'Main behavior is not defined');
      end
    end  % method
    
    % ---------------------------------------------------------------------
    function setMainBehaviorName(self,behaviorName)
      self.behaviors.names={behaviorName 'None'};
    end  % method

    % ---------------------------------------------------------------------
    function setScoreFileName(self,scoreFileName)
      self.file.scorefilename = scoreFileName;
    end
    
    % ---------------------------------------------------------------------
    function setTrxFileName(self,trxFileName)
      self.file.trxfilename = trxFileName;
    end
    
    % ---------------------------------------------------------------------
    function setMovieFileName(self,movieFileName)
      self.file.moviefilename = movieFileName;
    end

    function setMovieIndexFileName(self,movieIndexFileName)
      self.file.movieindexfilename = movieIndexFileName;
    end

    
    function addversion(self)
        vid = fopen('version.txt','r');
        vv = textscan(vid,'%s');
        fclose(vid);
        self.version = vv{1}{1};
    end
    
    function modernize(self,dowarn)
      if ~exist('dowarn','var')
        dowarn = false;
      end
      
      for i = 1:numel(self)        
        obj = self(i);
        
        if isscalar(obj.behaviors.names)
          % Legacy (pre v0.5.2?)
          assert(~strcmpi(obj.behaviors.names,'none'));
          obj.behaviors.names{end+1} = 'None';
        end
          
        if isfield(obj.behaviors,'nbeh')
          assert(obj.behaviors.nbeh==numel(obj.behaviors.names)/2);
        else
          obj.behaviors.nbeh = numel(obj.behaviors.names)/2;
        end
        
        FILE_OBSOLETE_FIELDS = {'gt_labelfilename' 'labelfilename' 'rootoutputdir' 'windowfilename'};
        fileflds = fieldnames(obj.file);
        obj.file = rmfield(obj.file,intersect(FILE_OBSOLETE_FIELDS,fileflds));
        
        if isempty(obj.scoreFeatures)
          % Standardize empty scoreFeatures shape
          obj.scoreFeatures = Macguffin.ScoreFeatures;
        end

        realbehs = Labels.verifyBehaviorNames(obj.behaviors.names);
        [obj.labels,tfmodlbl] = Labels.modernizeLabels(obj.labels,realbehs);
        [obj.gtLabels,tfmodGTlbl] = Labels.modernizeLabels(obj.gtLabels,realbehs);
        tfmodcls = obj.classifierStuff.modernize();
        % tfmodtags = isequal(obj.expDirTags,[]);
        % if tfmodtags
        %   obj.expDirTags = ExperimentTags.expTags(obj.expDirNames);
        % end
        
        if ischar(obj.file.scorefilename)
          obj.file.scorefilename = {obj.file.scorefilename};
        end          
        if isstruct(obj.windowFeaturesParams)
          obj.windowFeaturesParams = {obj.windowFeaturesParams};
        end
        nCls = numel(obj.classifierStuff);
        assert(isequal(nCls,numel(obj.file.scorefilename),...
                       numel(obj.windowFeaturesParams)));
                             
        tfmod = tfmodlbl || tfmodGTlbl || tfmodcls; % || tfmodtags;
        if dowarn && tfmod
          warningNoTrace('Macguffin:modernized','Jab contents modernized. Opening and resaving jabfiles in JAABA will update them and eliminate this warning.');
        end
      end
    end
        
  end  % methods
  
  methods (Static) % ScoreFeatures
  
    function sf = ScoreFeatures
      % ScoreFeatures constructor
      % 
      % Currently we default to a 0x0 struct but a 1x0 or 0x1 probably
      % makes more sense
          
      sf = struct('classifierfile',{},'ts',{},'scorefilename',{});
    end
    
  end
  
end
