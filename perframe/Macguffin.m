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
    labelGraphicParams
    scoreFeatures
    sublexiconPFNames
    windowFeaturesParams
    labels
    gtLabels
    expDirNames
    gtExpDirNames
    classifierStuff
    extra=struct()  % a structure that stores additional information
    version=''
      % 0.5.0 : Original
      % 0.5.1 : Supports nextra_markers, flies_extra_markersize, 
      %           flies_extra_marker, and flies_extra_linestyle fields in
      %           trxGraphicParams.
  end  % properties
    
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
      self.behaviors.labelcolors = [0.7,0,0,0,0,0.7];
      self.behaviors.unknowncolor = [0,0,0];
      self.trxGraphicParams=trxGraphicParamsFromAnimalType(animalType);
      self.labelGraphicParams.colormap = 'line';
      self.labelGraphicParams.linewidth = 3;
      self.file.scorefilename = '';
      self.file.trxfilename = '';
      self.file.moviefilename = '';
      self.scoreFeatures = struct('classifierfile',{},'ts',{},'scorefilename',{});
      self.featureLexicon=featureLexicon;
      featureLexiconPFNames = fieldnames(featureLexicon.perframe);
      self.sublexiconPFNames = featureLexiconPFNames;
      self.windowFeaturesParams={struct()};  % scalar struct with no fields
        % This is valid b/c no per-frame features have been enabled yet
      self.labels=struct('t0s',{}, ...
                         't1s',{}, ...
                         'names',{}, ...
                         'flies',{}, ...
                         'off',{}, ...
                         'timestamp',{});
      self.gtLabels=struct('t0s',{}, ...
                           't1s',{}, ...
                           'names',{}, ...
                           'flies',{}, ...
                           'off',{}, ...
                           'timestamp',{});
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
      self.file.trxfilename=jld.trxfilename;
      self.file.scorefilename=jld.scorefilename;
      self.file.clipsdir=jld.clipsdir;                  
      self.file.perframedir=jld.perframedir;                  
      self.labelGraphicParams=jld.labelGraphicParams;
      self.trxGraphicParams=jld.trxGraphicParams;
      if isnonempty(jld.landmark_params) ,
        self.extra.perframe.landmarkParams=jld.landmark_params;
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
      if isfield(projectParams,'labels')
        self.labelGraphicParams=projectParams.labels;
      else
        self.labelGraphicParams=projectParams.plot.labels;
      end
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
      % ALTODO: Caveat emptor, unverified functionality
      
      assert(false,'ALTODO MERGEST update');
              
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
      
      FIELDSTHATMIGHTDIFFER = {'behaviors' 'file' 'labels' 'expDirNames' 'classifierStuff'};
      fieldsMustBeSame = setdiff(flds,FIELDSTHATMIGHTDIFFER);       
      for f = fieldsMustBeSame(:)', f=f{1}; %#ok<FXSET>
        files = {m.(f)};
        if ~isequaln(files{:}),
          error('Macguffin:incompatibleObjs','Input Macguffins differ unexpectedly in field ''%s''.',f);
        end
        self.(f) = files{1};
      end      
      
      Nobj = numel(m);
      
      %%% files
      % all fields but scorefilename must be equal
      files = {m.file};
      filesNoSfn = cellfun(@(x)rmfield(x,'scorefilename'),files,'uni',0);
      if ~isequaln(filesNoSfn{:})
        error('Macguffin:incompatibleObjs','Input Macguffins differ unexpectedly in field ''.file''.');
      end
      % combined .file.scorefilename is concatenated cellstr
      self.file = files{1};
      scoreFileNames = cellfun(@(x)cellstr(x.scorefilename),files,'uni',0); % x.scorefilename might be char for single-classifier jab
      scoreFileNames = cellfun(@(x)x(:)',scoreFileNames,'uni',0);
      scoreFileNames = cat(2,scoreFileNames{:});
      if numel(unique(scoreFileNames))~=numel(scoreFileNames)
        error('Macguffin:repeatedScoreFile','Macguffins contain duplicate scorefile names.');
      end
      self.file.scorefilename = scoreFileNames;

      %%% behaviors/labels
      allbehnames = cell(1,0);
      alllabels = cell(Nobj,1); % contains modified .labels struct for each MacGuffin
      nbehs = nan(Nobj,1);
      for i = 1:Nobj                
        % check no duped behaviors
        % compile combined behavior list        
        bhvrs = m(i).behaviors;
        if ~strcmp(bhvrs.type,m(1).behaviors.type)
          error('Macguffin:incompatibleObjs','Input Macguffins differ unexpectedly in field ''%s''.','behaviors');
        end
        [behs,nobehs] = Labels.verifyBehaviorNames(bhvrs.names);
        if any(ismember(lower(behs),lower(allbehnames)))
          error('Macguffin:dupBehavior','Input Macguffins contain duplicate behaviors.');
        end 
        allbehnames = [allbehnames behs]; %#ok<AGROW>        
        nbehs(i) = numel(behs); % record for comparison to classifierStuff
        
        % convert all 'None's to No_<beh> in labels
        lbls = m(i).labels;
        if numel(behs)==1 && strcmpi(nobehs{1},'none')
          lbls = Labels.renameBehavior(lbls,nobehs{1},Labels.noBehaviorName(behs{1}),behs{1},behs{1});
        end
        alllabels{i} = lbls;
      end
      % set combined .behaviors 
      self.behaviors.type = m(1).behaviors.type;
      self.behaviors.names = [allbehnames cellfun(@Labels.noBehaviorName,allbehnames,'uni',0)];
      self.behaviors.labelcolors = Labels.cropOrAugmentLabelColors(zeros(1,0),numel(self.behaviors.names));
      self.behaviors.labelcolors = Labels.addNoBehColors(self.behaviors.labelcolors);
      self.behaviors.labelcolors = reshape(self.behaviors.labelcolors,3,[])';
      self.behaviors.unknowncolor = m(1).behaviors.unknowncolor;
      self.behaviors.nbeh = numel(allbehnames);      
      
      %%% expdirnames
      % - compile combined expdir list
      allexpdirnames = cat(2,m.expDirNames);
      allexpdirnames = unique(allexpdirnames);      
      self.expDirNames = allexpdirnames;
      warnNoTrace('Macguffin:discardingExpTags','Discarding existing experiment tags in merged jab.');
      % TODO: could be smart about combining expdirtags
      self.expDirTags = ExperimentTags.expTags(self.expDirNames); 
      % compile/combine labels
      self.labels = Labels.compileLabels(allexpdirnames,alllabels,{m.expDirNames});      
        
      %%% classifierStuff
      % TODO Consider moving into ClassifierStuff
      CSFIELDS_SAME = {'type' 'postProcessParams' 'featureNames' 'windowdata' 'savewindowdata'};
      selfCS = ClassifierStuff();
      mCS = [m.classifierStuff];
      for f = CSFIELDS_SAME,f=f{1}; %#ok<FXSET>
        val = {mCS.(f)};
        if ~isequaln(val{:})
          error('Macguffin:incompatibleObjs','Input Macguffins'' classifierStuff differ unexpectedly in field ''%s''.',f);
        end
        selfCS.(f) = val{1};
      end
      FLDS = {'params' 'trainingParams'};
      for f = FLDS,f=f{1}; %#ok<FXSET>
        for i = 1:numel(mCS)
          if ~iscell(mCS(i).(f))
            % not sure if this ever occurs
            mCS(i).(f) = {mCS(i).(f)};
          end
          assert(isrow(mCS(i).(f)),'Expected row vector for classifier.%s',f);
          assert(numel(mCS(i).(f))==nbehs(i),'Mismatch between behaviors and classifier.%s.',f);
        end
        selfCS.(f) = cat(2,mCS.(f));
        assert(numel(selfCS.(f))==self.behaviors.nbeh);
      end
      
      % AL 20140826: new classifierStuffs have independent timestamps
      % for each classifier; old classifierStuffs have a single timestamp
      % even if multiple behaviors. Rather than sort this out, just set
      % classifier timestamps in new (merged) jab to be zeros.
      selfCS.timeStamp = zeros(1,self.behaviors.nbeh);
      
      allConfThresh = nan(0,2);
      for i = 1:Nobj
        ct = m(i).classifierStuff.confThresholds;
        nbeh = numel(Labels.verifyBehaviorNames(m(i).behaviors.names));
        assert(isequal(size(ct),[1 2*nbeh])); % [beh1 beh2 ... No_beh1 No_beh2...] (I think)
        ct = reshape(ct,[],2);
        allConfThresh = [allConfThresh;ct]; %#ok<AGROW>
      end
      selfCS.confThresholds = allConfThresh(:)'; % all "real behavior" confThresholds, then all "No behavior" confThresholds
      assert(self.behaviors.nbeh*2==numel(selfCS.confThresholds));
      selfCS.scoreNorm = [mCS.scoreNorm];
      self.classifierStuff = selfCS;
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
      self.behaviors.names={behaviorName};
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
    
    function addversion(self)
        vid = fopen('version.txt','r');
        vv = textscan(vid,'%s');
        fclose(vid);
        self.version = vv{1};
    end
    
    function modernize(self,dowarn)
      if ~exist('dowarn','var')
        dowarn = false;
      end
      
      for i = 1:numel(self)
        [self(i).labels,tfmodlbl] = Labels.modernizeLabels(self(i).labels,Labels.verifyBehaviorNames(self(i).behaviors.names));        
        tfmodcls = self(i).classifierStuff.modernize();
        % tfmodtags = isequal(self(i).expDirTags,[]);
        % if tfmodtags
        %   self(i).expDirTags = ExperimentTags.expTags(self(i).expDirNames);
        % end
        
        if ischar(self(i).file.scorefilename)
          self(i).file.scorefilename = {self(i).file.scorefilename};
        end          
        if isstruct(self(i).windowFeaturesParams)
          self(i).windowFeaturesParams = {self(i).windowFeaturesParams};
        end
        nCls = numel(self(i).classifierStuff);
        assert(isequal(nCls,numel(self(i).file.scorefilename),...
                       numel(self(i).windowFeaturesParams)));
        
        tfmod = tfmodlbl || tfmodcls; % || tfmodtags;
        if dowarn && tfmod
          warningNoTrace('Macguffin:modernized','Jab contents modernized. Opening and resaving jabfiles in JAABA will update them and eliminate this warning.');
        end
      end
    end
  end  % methods
    
  methods (Static) % jabfile convenience methods
    
    function jabClearLabels(jab,realbeh)
      % Clear all labels for realbeh and No-realbeh 
      
      assert(ischar(jab) && exist(jab,'file')==2,...
        'Cannot find jabfile ''%s''.',jab);
            
      Q = loadAnonymous(jab);
      Q.modernize();
      
      jabrealbehs = Q.behaviors.names(1:Q.behaviors.nbeh);
      tf = strcmpi(realbeh,jabrealbehs);
      assert(any(tf),'Specified behavior ''%s'' not present in jabfile.',realbeh);
      assert(nnz(tf)==1);
      realbeh = jabrealbehs{tf}; % case could be different
      
      tfMultiCls = Q.behaviors.nbeh>1;
      if tfMultiCls
        nobeh = Labels.noBehaviorName(realbeh);
      else
        nobeh = 'None';
      end
      
      Q.labels = Labels.clearLabels(Q.labels,realbeh,realbeh);
      Q.labels = Labels.clearLabels(Q.labels,nobeh,realbeh);
      
      saveAnonymous(jab,Q);
    end
     
    function jabMerge(jabfiles,jabout)
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
    
    function tf = jabfileIsMultiClassifier(jabname)
      % tf = jabfileIsMultiClassifier(filename)
      x = loadAnonymous(jabname);
      tf = x.isMultiClassifier();      
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
  
end  % classdef
