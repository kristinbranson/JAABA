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
    version='0.5.0'
  end  % properties
  
  
  % -----------------------------------------------------------------------
  methods (Access=private)
    % ---------------------------------------------------------------------
    function initFromFeatureLexiconName(self,featureLexiconName)
      self.featureLexiconName=featureLexiconName;
      self.behaviors.names = {};  % no behaviors defined yet
      self.file.perframedir = 'perframe';
      self.file.clipsdir = 'clips';
      %self.windowfeatures = struct;
      self.behaviors.labelcolors = [0.7,0,0,0,0,0.7];
      self.behaviors.unknowncolor = [0,0,0];
      self.trxGraphicParams.colormap = 'jet';
      self.trxGraphicParams.colormap_multiplier = 0.7;
      self.trxGraphicParams.extra_linestyle = '-';
      self.trxGraphicParams.extra_marker = '.';
      self.trxGraphicParams.extra_markersize = 12;
      self.labelGraphicParams.colormap = 'line';
      self.labelGraphicParams.linewidth = 3;
      self.file.scorefilename = '';
      self.file.trxfilename = '';
      self.file.moviefilename = '';
      self.scoreFeatures = struct('classifierfile',{},'ts',{},'scorefilename',{});
      featureLexicon=featureLexiconFromFeatureLexiconName(featureLexiconName);
      self.featureLexicon=featureLexicon;
      featureLexiconPFNames = fieldnames(featureLexicon.perframe);
      self.sublexiconPFNames = featureLexiconPFNames;
      self.windowFeaturesParams=struct();  % scalar sctruct with no fields
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
    end  % method
    
    
    % ---------------------------------------------------------------------
    function initFromJLabelData(self,jld)
      self.featureLexiconName=jld.featureLexiconName;
      self.featureLexicon=jld.featureLexicon;
      self.scoreFeatures=jld.scoreFeatures;
      subdialectPFNames=jld.allperframefns;
      nScoreFeaturess=length(jld.scoreFeatures);
      sublexiconPFNames=subdialectPFNames(1:end-nScoreFeaturess);
      self.sublexiconPFNames=sublexiconPFNames;
      self.behaviors.type=jld.targettype;
      self.behaviors.names=jld.labelnames;
      self.behaviors.labelcolors=jld.labelcolors;
      self.behaviors.unknowncolor=jld.unknowncolor;
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
      [self.labels,self.gtLabels]=jld.storeAndGetLabelsAndGTLabels();
      % Get the window feature params, put in self
      self.windowFeaturesParams=jld.windowfeaturesparams;
      % Put the classifier in self
      self.classifierStuff=jld.getClassifierStuff();
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
    end  % method

    
    % ---------------------------------------------------------------------
    function initFromOldStyleProjectAndClassifier(self, ...
                                                  projectParams, ...
                                                  classifierParams, ...
                                                  gtExpDirNames, ...
                                                  gtLabels)
    
      % get the featureLexicon from the relevant file
      featureLexiconFileNameRel=projectParams.file.featureconfigfile;
      pathToMisc=fileparts(mfilename());
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
        self.trxGraphicParams=projectParams.trx;
      else
        self.trxGraphicParams=projectParams.plot.trx;
      end
      if isfield(projectParams,'labels')
        self.labelGraphicParams=projectParams.labels;
      else
        self.labelGraphicParams=projectParams.plot.labels;
      end
      self.featureLexiconName=featureLexiconName;
      self.featureLexicon=featureLexicon;
      self.sublexiconPFNames=sublexiconPFNames;
      self.scoreFeatures=projectParams.scoresinput;
%       if isempty(projectParams.behaviors.names) ,
%         self.behaviorName='';
%       else
%         self.behaviorName=projectParams.behaviors.names{1};
%       end
      if isfield(projectParams,'landmark_params')
        self.extra.perframe.landmarkParams=projectParams.landmark_params;
      else
        self.extra.perframe.landmarkParams=[];
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
      self.gtLabels=gtLabels;
    end  % method
    
    
    % ---------------------------------------------------------------------
    function appendClassifierAndLabels(self,projectParams,classifierParams)
      % this is used when converting old-style files to .jab files
      self.expDirNames=classifierParams.expdirs;
      self.labels=classifierParams.labels;
      if isfield(projectParams.windowfeatures,'windowfeaturesparams')
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
      if length(varargin)==1 && ischar(varargin{1})
        self.initFromFeatureLexiconName(varargin{1});
      elseif length(varargin)==1 && isequal(class(varargin{1}),'JLabelData')
        jld=varargin{1};
        self.initFromJLabelData(jld);
      elseif length(varargin)==1 && isstruct(varargin{1})
        basicDataStruct=varargin{1};
        self.initFromBasicDataStruct(basicDataStruct);
      elseif length(varargin)==4
        projectParams=varargin{1};
        classifierParams=varargin{2};
        gtExpDirNames=varargin{3};
        gtLabels=varargin{4};
        self.initFromOldStyleProjectAndClassifier(projectParams, ...
                                                  classifierParams, ...
                                                  gtExpDirNames, ...
                                                  gtLabels);
      else
        error('Macguffin:badArgumentsToConstructor', ...
              'The arguments to Macguffin() are no good');
      end
    end  % constructor method
  end  % methods
end  % classdef
