classdef JLabelData < handle
  
  properties (Access=public)

    % type of target (mainly used for plotting
    targettype = 'fly';
    
    % current selection
    
    % currently selected  experiment
    expi = 0;
    % currently selected flies
    flies = [];
    
    % last-used trajectories (one experiment, all flies)
    trx = {};

    % last-used per-frame data (one fly)
    perframedata = {};
    
    % computed and cached window features
    windowdata = struct('X',single([]),'exp',[],'flies',[],'t',[],...
      'labelidx_cur',[],'labelidx_new',[],'labelidx_old',[],...
      'labelidx_imp',[],'featurenames',{{}},...
      'predicted',[],'predicted_probs',[],'isvalidprediction',[],...
      'distNdx',[],'scores',[],'scoreNorm',[],'binVals',[],...
      'scores_old',[],'scores_validated',[],'postprocessed',[]);

    % Score loaded from files. Scores are loaded scores,
    % cur_scores are scores of the current classifier and old_scores are
    % scores of the old classifier.
    % 
    predictdata = struct('exp',[],'flies',[],'t',[],...
          'cur',[],'cur_valid',logical([]),'cur_pp',[],...
          'old',[],'old_valid',logical([]),'old_pp',[],...
          'loaded',[],'loaded_valid',logical([]),'loaded_pp',[],...
          'timestamp',[]);
        
    predictblocks = struct('t0',[],'t1',[],'expi',[],'flies',[]);
   
    fastPredict = struct('classifier',[],...
          'windowfeaturescellparams',[],...
          'wfs',{{}},...
          'pffs',[],'wfidx',[],'ts',[],...
          'wfidx_valid',false);
    
    % constant: radius of window data to compute at a time
    windowdatachunk_radius = 100;
    predictwindowdatachunk_radius = 10000;
    % total number of experiments
    nexps = 0;
    
    % labels struct array
    % labels(expi) is the labeled data for experiment expi
    % labels(expi).t0s are the start frames of all labeled sequences for
    % experiment expi
    % labels(expi).t1s are the corresponding end frames of all labeled
    % sequences for experiment expi
    % labels(expi).names is the cell array of the corresponding behavior
    % names for all labeled sequences for experiment expi
    % labels(expi).flies is the nseq x nflies_labeled matrix of the
    % corresponding flies for all labeled sequences for experiment expi
    % t0s{j}, t1s{j}, names{j}, and flies(j,:) correspond to each other. 
    % labels(expi).off is the offset so that labels(expi).t0s(j) +
    % labels(expi).off corresponds to the frame of the movie (since the
    % first frame for the trajectory(s) may not be 1.
    % labels(expi).timestamp is the Matlab timestamp at which labels(expi)
    % was last set
    labels = struct('t0s',{},'t1s',{},'names',{},'flies',{},'off',{},'timestamp',{},'imp_t0s',{},'imp_t1s',{});
    gt_labels = struct('t0s',{},'t1s',{},'names',{},'flies',{},'off',{},'timestamp',{},'imp_t0s',{},'imp_t1s',{});
    
    % labels for the current experiment and flies, represented as an array
    % such that labelidx(t+labelidx_off) is the index of the behavior for
    % frame t of the movie. labelidx(i) == 0 corresponds to
    % unlabeled/unknown, otherwise labelidx(i) corresponds to behavior
    % labelnames{labelidx{i})
    labelidx = struct('val',[],'imp',[],'timestamp',[]);
    labelidx_off = 0;
    
    labelsLoadedFromClassifier = false;
    % first frame that all flies currently selected are tracked
    t0_curr = 0;
    % last frame that all flies currently selected are tracked
    t1_curr = 0;
    
    % predicted label for current experiment and flies, with the same type
    % of representation as labelidx
    predictedidx = [];
    scoresidx = [];
    scoresidx_old = [];
    scoreTS = [];
    
    % whether the predicted label matches the true label. 0 stands for
    % either not predicted or not labeled, 1 for matching, 2 for not
    % matching. this has the same representation as labelidx.
    erroridx = [];
    
    % TODO: remove this
    % predictedidx for unlabeled data, same representation as labelidx
    suggestedidx = [];
    
    % names of behaviors, corresponding to labelidx
    labelnames = {};
    
%     % colors for plotting each behavior
%     labelcolors = [.7,0,0;0,0,.7];
%     unknowncolor = [0,0,0];
    
    % number of behaviors, including 'none'
    nbehaviors = 0;

    % statistics of labeled data per experiment
    % labelstats(expi).nflies_labeled is the total number of flies labeled,
    % labelstats(expi).nbouts_labeled is the total number of bouts of
    % behaviors labeled, labelstats(expi).
    labelstats = struct('nflies_labeled',{},'nbouts_labeled',{});
    gt_labelstats = struct('nflies_labeled',{},'nbouts_labeled',{});
    
    % computing per-frame properties
    perframe_params = {};
    landmark_params = struct;
    
    % classifier
    
    % type of classifier to use
%    classifiertype = 'ferns';
    classifiertype = 'boosting';
    
    % currently learned classifier. structure depends on the type of
    % classifier. if empty, then no classifier has been trained yet. 
    % ferns:
    % classifier is a struct with the following fields (M is the number of
    % ferns, S is the fern depth, N is the number of training examples). 
    %   .fids     - [MxS] feature ids for each fern for each depth
    %   .thrs     - [MxS] threshold corresponding to each fid
    %   .pFern    - [2^SxHxM] learned log probs at fern leaves
    %   .bayes    - if true combine probs using bayes assumption
    %   .inds     - [NxM] cached indices for original training data
    %   .H        - number classes
    classifier = [];
    classifier_old = [];
    lastFullClassifierTrainingSize = 0;
    
    % Classifiers Time Stamp
    classifierTS = 0;
    
    % parameters to learning the classifier. struct fields depend on type
    % of classifier.
    % TODO
    classifier_params = struct('iter',100,'iter_updates',10,...
      'numSample',2500,'numBins',30,'CVfolds',7,...
      'baseClassifierTypes',{'Decision Stumps'},'baseClassifierSelected',1);
    
    % stuff cached during prediction
    predict_cache = struct;
    
    % name of file containing config parameters
    configfilename = '';
    
    % constant: files per experiment directory
    filetypes = {'movie','trx','label','gt_label','perframedir','clipsdir','scores'};
    
    % config parameters
    
    % locations of files within experiment directories
    moviefilename = 0;
    trxfilename = 0;
    labelfilename = 0;
    gt_labelfilename = 0;
    scorefilename = 0;
    perframedir = 0;
    clipsdir = 0;
    scores = 0;
    
    % whether there is a movie to show
    ismovie = true;
    openmovie = true;
    
    % file containing feature parameters
    featureparamsfilename = 0;
    
    % feature configuration file
    featureConfigFile = '';
    
    % in case we don't want to write to the experiment directory, we will
    % mirror the experiment directory structure in the rootoutput dir
    % this can be the same as the input root directory
    rootoutputdir = 0;

    % name of classifier file to save/load classifier from
    classifierfilename = '';
    
    % experiment info: expi indexes the following
    
    % cell array of input experiment directory paths
    expdirs = {};
    
    % cell array of corresponding experiment names (last part of path)
    expnames = {};
    
    % cell array of corresponding output experiment directory paths
    outexpdirs = {};
    
    % array of number of flies in each experiment
    nflies_per_exp = [];
    
    % cell array of arrays of first frame of each trajectory for each
    % experiment: firstframes_per_exp{expi}(fly) is the first frame of the
    % trajectory of fly for experiment expi. 
    firstframes_per_exp = {};
    % cell array of arrays of end frame of each trajectory for each
    % experiment: endframes_per_exp{expi}(fly) is the last frame of the
    % trajectory of fly for experiment expi. 
    endframes_per_exp = {};
    
    % sex per experiment, fly
    frac_sex_per_exp = {};
    sex_per_exp = {};
    
    % whether sex is computed
    hassex = false;
    % whether sex is computed on a per-frame basis
    hasperframesex = false;
    
    % constant: stuff stored in classifier mat file
    classifiervars = {'expdirs','outexpdirs','expnames',...
      'nflies_per_exp','sex_per_exp','frac_sex_per_exp',...
      'firstframes_per_exp','endframes_per_exp',...
      'moviefilename','trxfilename','labelfilename','perframedir','clipsdir',...,'featureparamsfilename',...
      'configfilename','rootoutputdir','classifiertype','classifier','trainingdata','classifier_params',...
      'classifierTS','confThresholds','scoreNorm','windowfeaturesparams','windowfeaturescellparams',...
      'basicFeatureTable','featureWindowSize','postprocessparams',...
      'featurenames','scorefilename','labels'};
    
    % last used path for loading experiment
    defaultpath = '';
    
    % parameters of window features, represented as a struct
    windowfeaturesparams = struct;
    
    % parameters of window features, represented as a cell array of
    % parameter name, parameter value, so that it can be input to
    % ComputeWindowFeatures
    windowfeaturescellparams = {};
    
    savewindowfeatures = false;
    
    % State of the basic/compact feature table.
    basicFeatureTable = {};
    featureWindowSize = [];
    
    % per-frame features that are used
    allperframefns = {};
    curperframefns = {};
    perframeunits = {};
    scoresasinput = [];
    
    % experiment/file management

    % matrix of size numel(file_types) x nexps, where
    % fileexists(filei,expi) indicates whether file filetypes{filei} exists
    % for experiment expi
    fileexists = [];
    
    % timestamps indicating time the files were last edited, same structure
    % as fileexists
    filetimestamps = [];
    
    % whether all necessary files for all experiments exist
    allfilesexist = true;

    % whether we can generate any missing files
    filesfixable = true;
    
    % whether user has given permission to generate the perframe files
    perframeGenerate = [];
    
    % to overwrite or keep the perframe files.
    perframeOverwrite = [];
    
    % warn about removing arena features.
    arenawarn = true;
    hasarenaparams = [];
    
    % functions for writing text to a status bar
    setstatusfn = '';
    clearstatusfn = '';
    
    % data for show similar frames.
    frameFig = [];
    distMat = [];
    bagModels = {};

    % Confidence Thresholds
    confThresholds = zeros(1,0);
    
    % Retrain properly
    doUpdate = true;
    
    % Ground truthing or not
    gtMode = false;
    advancedMode = false;
    modeSet = false;
    
    % Ground truthing suggestion
    randomGTSuggestions = {};
    thresholdGTSuggestions = [];
    loadedGTSuggestions = {};
    balancedGTSuggestions = {};
    GTSuggestionMode = '';
    
    cacheSize = 4000;
    
    postprocessparams = []
end
  
  methods (Access=private)
    
    
  end
    
  methods (Access=public,Static=true)

    % movie, trx, and perframedir are required for each experiment
    function res = IsRequiredFile(file)
      res = ismember(file,{'movie','trx','perframedir'});
    end    
    
    % perframedir can be generated
    function res = CanGenerateFile(file)
      res = ismember(file,{'perframedir'});
    end
    
    % which files should go in the output directory
    function res = IsOutputFile(file)
      res = ismember(file,{'label','clipsdir','scores','gt_label'});
    end
    
    
    function valid = CheckExp(expi)
      if numel(expi) ~= 1,
        error('Usage: expi must be a scalar');
        valid = false;
      else
        valid = true;
      end
    end
    
    function valid = CheckFlies(flies)
      if size(flies,1) ~= 1,
        error('Usage: one set of flies must be selected');
        valid = false;
      else
        valid = true;
      end
    end      

    function [X,feature_names] = ...
        ComputeWindowDataChunkStatic(curperframefns,allperframefns,perframefile,flies,windowfeaturescellparams,t0,t1)
      
    %function [X,feature_names] = ...
    %    ComputeWindowDataChunkStatic(perframefns,perframedir,flies,windowfeaturecellparams,t0,t1)
    %
    % Computes a chunk of windowdata between frames t0 and t1 for flies
    % flies. 
    %
    % X is the nframes x nfeatures window data, feature_names is a cell array
    % of length nfeatures containing the names of each feature. 
    %
    % It then loops through all the per-frame features, and calls
    % ComputeWindowFeatures to compute all the window data for that
    % per-frame feature.
    
    X = [];
    feature_names = {};
    
    for j = 1:numel(curperframefns),
      fn = curperframefns{j};
      ndx = find(strcmp(fn,allperframefns));

      perframedata = load(perframefile{ndx});
      perframedata = perframedata.data{flies(1)};
      
      t11 = min(t1,numel(perframedata));
      [x_curr,feature_names_curr] = ...
        ComputeWindowFeatures(perframedata,windowfeaturescellparams.(fn){:},'t0',t0,'t1',t11);
      
      if t11 < t1,
        x_curr(:,end+1:end+t1-t11) = nan;
      end
      
      % add the window data for this per-frame feature to X
      nold = size(X,1);
      nnew = size(x_curr,2);
      if nold > nnew,
        warning('Number of examples for per-frame feature %s does not match number of examples for previous features',fn);
        x_curr(:,end+1:end+nold-nnew) = nan;
      elseif nnew > nold && ~isempty(X),
        warning('Number of examples for per-frame feature %s does not match number of examples for previous features',fn);
        X(end+1:end+nnew-nold,:) = nan;
      end
      X = [X,x_curr']; %#ok<AGROW>
      % add the feature names
      feature_names = [feature_names,cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false)]; %#ok<AGROW>
    end
    X = single(X);
    
    end
    
    function params = convertTransTypes2Cell(params)
      % Convert the trans_types field into cell type
      if ~isstruct(params),  return; end
      fnames = fieldnames(params);
      for ndx = 1:numel(fnames)
        if isstruct(params.(fnames{ndx}))
          params.(fnames{ndx}) = JLabelData.convertTransTypes2Cell(params.(fnames{ndx}));
        end
      end
      if isfield(params,'trans_types')&& ~iscell(params.trans_types)
        params.trans_types = {params.trans_types};
      end
    end
    

    
  end
  
  methods (Access=public)


% Constructor


    function obj = JLabelData(varargin)
    % obj = JLabelData(configfilename,...)
    %
    % constructor: first input should be the config file name. All other
    % inputs are optional. if configfilename is not input, user will be
    % prompted for it. 
    % 
    % optional inputs: 
    %
    % TODO: debug this
    % override stuff set in the config file: 
    %
    % moviefilename, trxfilename, labelfilename, perframedir, clipsdir: names of
    % files within experiment directories: 
    % featureparamsfilename: file containing feature parameters
    % rootoutputdir: in case we don't want to write to the experiment
    % directory, we will mirror the experiment directory structure in the
    % rootoutputdir this can be the same as the input root directory
    %
    % defaultpath: default location to look for experiments
    % setstatusfn: handle to function that inputs sprintf-like input and
    % outputs the corresponding string to a status bar.
    % clearstatusfn: handle to function that returns the status to the
    % default string
    % classifierfilename: name of classifier file to save/load classifier from
 
      if nargin == 0 || isempty(varargin{1}),
        [filename,pathname] = uigetfile('*.xml','Choose config XML file');
        if ~ischar(filename),
          return;
        end
        configfilename = fullfile(pathname,filename);
        if ~isempty(varargin),
          varargin = varargin(2:end);
        end
      else
        configfilename = varargin{1};
        varargin = varargin(2:end);
      end
      
      if mod(numel(varargin),2) ~= 0,
        error('Number of inputs to JLabelData constructor must be even');
      end
      
      % config file
      [success,msg] = obj.SetConfigFileName(configfilename);
      if ~success,
        error(msg);
      end
      
      % parse optional arguments in order
      s = varargin(1:2:end);
      v = varargin(2:2:end);
      
      i = find(strcmpi(s,'openmovie'),1);
      if ~isempty(i),
        obj.openmovie = v{i};
      end
      
      % movie
      i = find(strcmpi(s,'moviefilename'),1);
      if ~isempty(i),
        [success,msg] = obj.SetMovieFileName(v{i});
        if ~success,
          error(msg);
        end

      end
      
      % trx
      i = find(strcmpi(s,'trxfilename'),1);
      if ~isempty(i),
        [success,msg] = obj.SetTrxFileName(v{i});
        if ~success,
          error(msg);
        end
      end
      
      % label
      i = find(strcmpi(s,'labelfilename'),1);
      if ~isempty(i),
        [success,msg] = obj.SetLabelFileName(v{i});
        if ~success,
          error(msg);
        end
      end

      % perframedir
      i = find(strcmpi(s,'perframedir'),1);
      if ~isempty(i),
        [success,msg] = obj.SetPerFrameDir(v{i});
        if ~success,
          error(msg);
        end
      end

      % clipsdir
      i = find(strcmpi(s,'clipsdir'),1);
      if ~isempty(i),
        [success,msg] = obj.SetClipsDir(v{i});
        if ~success,
          error(msg);
        end
      end

      % featureparamsfilename
      i = find(strcmpi(s,'featureparamsfilename'),1);
      if ~isempty(i),
        [success,msg] = obj.SetFeatureParamsFileName(v{i});
        if ~success,
          error(msg);
        end
      end
      
      % rootoutputdir
      i = find(strcmpi(s,'rootoutputdir'),1);
      if ~isempty(i),
        [success,msg] = obj.SetRootOutputDir(v{i});
        if ~success,
          error(msg);
        end
      end
      
      % classifier
      i = find(strcmpi(s,'classifierfilename'),1);
      if ~isempty(i),
        [success,msg] = obj.SetClassifierFileName(v{i});
        if ~success,
          error(msg);
        end
      end
      
      % default path
      i = find(strcmpi(s,'defaultpath'),1);
      if ~isempty(i),
        [success,msg] = obj.SetDefaultPath(v{i});
        if ~success,
          warning(msg);
        end
      end
      
      i = find(strcmpi(s,'setstatusfn'),1);
      if ~isempty(i),
        obj.setstatusfn = v{i};
      end

      i = find(strcmpi(s,'clearstatusfn'),1);
      if ~isempty(i),
        obj.clearstatusfn = v{i};
      end
      
      i = find(strcmpi(s,'cacheSize'),1);
      if ~isempty(i),
        obj.cacheSize = v{i};
      end
      
      % make sure everything gets set, one way or another
      requiredfns = {'moviefilename','trxfilename','labelfilename'}; % 'windowfilename',
      for i = 1:numel(requiredfns),
        fn = requiredfns{i};
        if isnumeric(obj.(fn)),
          error('%s did not get initialized',fn);
        end
      end
      
      % initialize the status table describing what required files exist
      [success,msg] = obj.UpdateStatusTable();
      if ~success,
        error(msg);
      end
      
      obj.InitPostprocessparams();
      
    end

    
% Some helper functions.


    function idx = FlyNdx(obj,expi,flies)
      if isempty(obj.windowdata.exp),
        idx = []; return;
      end
      idx = obj.windowdata.exp == expi & all(bsxfun(@eq,obj.windowdata.flies,flies),2);
    end
    
    function idx = FlyNdxPredict(obj,expi,flies)
      if isempty(obj.predictdata.exp),
        idx = []; return;
      end
      idx = obj.predictdata.exp == expi;
      idx(obj.predictdata.flies~=flies) = false;
    end
    
    function val = IsCurFly(obj,expi,flies)
      val = all(flies == obj.flies) && (expi==obj.expi);
    end
    
    function expi = GetExp(obj)
      expi = obj.expi;
    end
    
    function flies = GetFlies(obj)
      flies = obj.flies;
    end
     
    function nflies = GetNumFlies(obj,expi)
      nflies = obj.nflies_per_exp(expi);
    end
    
    
% Configuration settings.


    function [success,msg] = SetConfigFileName(obj,configfilename)
      % [success,msg] = SetConfigFileName(obj,configfilename)
      % Set and read config file.
      % Reads the XML config file, then sets all the file names and paths.
      % I think this currently needs to be called before experiments, labels
      % are loaded in, as locations of files, behaviors can be modified by
      % this step.
      % labelnames, nbehaviors are also set by the config file. If not
      % included explicitly, the 'None' behavior is added. 'None' is put at
      % the end of the behavior list.
      
      success = false;
      msg = '';
      if ~ischar(configfilename),
        return;
      end
      %       try
      [~,~,ext] = fileparts(configfilename);
      if strcmp(ext,'.xml'),
        configparams = ReadXMLParams(configfilename);
      else
        configparams = load(configfilename);
      end
      %       catch ME,
      %         msg = sprintf('Error reading config file %s: %s',configfilename,getReport(ME));
      %         return;
      %       end
      obj.configfilename = configfilename;
      if isfield(configparams,'behaviors'),
        
        % read in behavior names
        if isfield(configparams.behaviors,'names'),
          obj.labelnames = configparams.behaviors.names;
          if ~iscell(obj.labelnames),
            obj.labelnames = {obj.labelnames};
          end
          % add none label
          if ~ismember('none',lower(obj.labelnames)),
            obj.labelnames{end+1} = 'None';
          end
        else
          obj.labelnames = {'Behavior','None'};
        end
        
        obj.nbehaviors = numel(obj.labelnames);
        % allocate configdence thresholds
        obj.confThresholds = zeros(1,obj.nbehaviors);

%         % colors
%         if isfield(configparams.behaviors,'labelcolors'),
%           if numel(configparams.behaviors.labelcolors) == obj.nbehaviors*3,
%             obj.labelcolors = configparams.behaviors.labelcolors;
%           end
%         end
%         if isfield(configparams.behaviors,'unknowncolor'),
%           if numel(configparams.behaviors.unknowncolor) == 3,
%             obj.unknowncolor = configparams.behaviors.unknowncolor;
%           end
%         end
        
        % rearrange so that None is the last label
        nonei = find(strcmpi('None',obj.labelnames),1);
        obj.labelnames = obj.labelnames([1:nonei-1,nonei+1:obj.nbehaviors,nonei]);
        
      end

      if isfield(configparams,'file'),
        if isfield(configparams.file,'moviefilename'),
          [success1,msg] = obj.SetMovieFileName(configparams.file.moviefilename);
          if ~success1,
            return;
          end
        end
        if isfield(configparams.file,'trxfilename'),
          [success1,msg] = obj.SetTrxFileName(configparams.file.trxfilename);
          if ~success1,
            return;
          end
        end
        if isfield(configparams.file,'labelfilename'),
          [success1,msg] = obj.SetLabelFileName(configparams.file.labelfilename);
          if ~success1,
            return;
          end
        end
        if isfield(configparams.file,'gt_labelfilename'),
          [success1,msg] = obj.SetGTLabelFileName(configparams.file.gt_labelfilename);
          if ~success1,
            return;
          end
        end
        if isfield(configparams.file,'scorefilename'),
          scorefilename = configparams.file.scorefilename;
        else
          scorefilename = sprintf('scores_%s.mat',obj.labelnames{1});
        end
        [success1,msg] = obj.SetScoresFileName(scorefilename);
        if ~success1,
          return;
        end
        
        if isfield(configparams.file,'perframedir'),
          [success1,msg] = obj.SetPerFrameDir(configparams.file.perframedir);
          if ~success1,
            return;
          end
        end
        if isfield(configparams.file,'clipsdir') && ~isempty(configparams.file.clipsdir),
          [success1,msg] = obj.SetClipsDir(configparams.file.clipsdir);
          if ~success1,
            return;
          end
        end
        if isfield(configparams.file,'rootoutputdir') && ~isempty(configparams.file.rootoutputdir),
          [success1,msg1] = obj.SetRootOutputDir(configparams.file.rootoutputdir);
          if ~success1,
            uiwait(warndlg(msg1));
          end
        end
        if isfield(configparams.file,'featureconfigfile'),
          [success1,msg] = obj.SetFeatureConfigFile(configparams.file.featureconfigfile);
          if ~success1,
            return;
          end
        end
        if isfield(configparams,'featureparamlist'),
          % read allperframefns from config file
          obj.allperframefns = intersect(obj.allperframefns,...
                              fieldnames(configparams.featureparamlist));
          msg = '';
        end

        if isfield(configparams.file,'featureparamfilename') && ~isempty(configparams.file.featureparamfilename),
          [success1,msg] = obj.SetFeatureParamsFileName(configparams.file.featureparamfilename);
          if ~success1,
            return;
          end
        end
        
        if isfield(configparams,'windowfeatures') && isfield(configparams.windowfeatures,'basicFeatureTable')
          obj.basicFeatureTable = configparams.windowfeatures.basicFeatureTable;
          obj.featureWindowSize = configparams.windowfeatures.featureWindowSize;
          obj.SetPerframeParams(configparams.windowfeatures.windowfeaturesparams,...
            configparams.windowfeatures.windowfeaturescellparams);
        end
        
        if isfield(configparams,'perframe'),
          if isfield(configparams.perframe,'params'),
            pf_fields = fieldnames(configparams.perframe.params);
            for ndx = 1:numel(pf_fields),
              obj.perframe_params.(pf_fields{ndx}) = configparams.perframe.params.(pf_fields{ndx});
            end
          end
          if isfield(configparams.perframe,'landmark_params'),
            obj.landmark_params = configparams.perframe.landmark_params;
          end
        end
        if isfield(configparams,'targets'),
          if isfield(configparams.targets,'type'),
            obj.targettype = configparams.targets.type;
          end
        end
      end
      
      
      if isfield(configparams,'learning'),
        if isfield(configparams.learning,'classifiertype'),
          obj.SetClassifierType(configparams.learning.classifiertype);
        end
      end
      
      if isfield(configparams,'scoresinput'),
        obj.scoresasinput = configparams.scoresinput;
        for ndx = 1:numel(obj.scoresasinput)
          [~,name,~] = fileparts(obj.scoresasinput(ndx).scorefilename);
          obj.allperframefns{end+1} = name;
        end
        
        if ~isempty(obj.basicFeatureTable),
          scoresbasicndx = find(strcmpi(obj.basicFeatureTable(:,1),'scores'));
          if isempty(scoresbasicndx),
            obj.basicFeatureTable(end+1,:) = {'scores','Custom','normal'};
          else
            obj.basicFeatureTable{scoresbasicndx,2} = 'Custom';
          end
        end
      end
      
     
    end
    
    function [success,msg] = SetMovieFileName(obj,moviefilename)
    % change/set the name of the movie within the experiment directory
    % will fail if movie files don't exist for any of the current
    % experiment directories (checked by CheckMovies)

      success = false; msg = '';

      if ischar(moviefilename),
        if ischar(obj.moviefilename) && strcmp(moviefilename,obj.moviefilename),
          success = true;
          return;
        end
        oldmoviefilename = obj.moviefilename;
        obj.moviefilename = moviefilename;
        obj.ismovie = ~isempty(moviefilename) && obj.openmovie;
        [success1,msg] = obj.CheckMovies();
        if ~success1,
          obj.moviefilename = oldmoviefilename;
          return;
        end
        [success,msg] = obj.UpdateStatusTable('movie');
      end
      
    end
    
    function [successes,msg] = CheckMovies(obj,expis)
    % [successes,msg] = CheckMovies(obj,expis)
    % check that the movie files exist and can be read for the input
    % experiments.
      
      successes = []; msg = '';
      
      if nargin < 2,
        expis = 1:obj.nexps;
      end
      
      if isempty(expis),
        return;
      end
      
      successes = true(1,numel(expis));
      
      if ~obj.ismovie,
        return;
      end
      
      for i = 1:numel(expis),
        moviefilename = obj.GetFile('movie',expis(i));
        obj.SetStatus('Checking movie %s...',moviefilename);
        
        % check for file existence
        if ~exist(moviefilename,'file'),
          successes(i) = false;
          msg1 = sprintf('File %s missing',moviefilename);
          if isempty(msg),
            msg = msg1;
          else
            msg = sprintf('%s\n%s',msg,msg1);
          end
        else
          
          % try reading a frame
%           try
            [readframe,~,movie_fid] = ...
              get_readframe_fcn(moviefilename);
            if movie_fid <= 0,
              error('Could not open movie %s for reading',moviefilename);
            end
            readframe(1);
            fclose(movie_fid);
%           catch ME,
%             successes(i) = false;
%             msg1 = sprintf('Could not parse movie %s: %s',moviefilename,getReport(ME));
%             if isempty(msg),
%               msg = msg1;
%             else
%               msg = sprintf('%s\n%s',msg,msg1);
%             end
%           end
          
        end
      end
      
      obj.ClearStatus();
      
    end
    
    function [success,msg] = SetTrxFileName(obj,trxfilename)
      % [success,msg] = SetTrxFileName(obj,trxfilename)
    % set the name of the trx file within the experiment directory. this
    % does not currently check for missing/bad trx files, or replace
    % preloaded trx data, so you really shouldn't call it if expdirs are
    % loaded. (TODO)

      success = false;
      msg = '';
      if ischar(trxfilename),
        if ischar(obj.trxfilename) && strcmp(trxfilename,obj.trxfilename),
          success = true;
          return;
        end
        obj.trxfilename = trxfilename;
        [success,msg] = obj.UpdateStatusTable('trx');        
        % TODO: check that trx are parsable, remove bad experiments, update
        % preloaded trx
      end
      
    end
    
    function [success,msg] = SetLabelFileName(obj,labelfilename)
    % [success,msg] = SetLabelFileName(obj,labelfilename)
    % set the name of the label file within the experiment directory. this
    % does not currently update labelidx, and probably should not be called
    % once an experiment is open. 
      
      success = false;
      msg = '';

      if ischar(labelfilename),
        if ischar(obj.labelfilename) && strcmp(labelfilename,obj.labelfilename),
          success = true;
          return;
        end

        % reload labels from file
        for expi = 1:obj.nexps,
          [success1,msg] = obj.LoadLabelsFromFile(expi);
          if ~success1,
            return;
          end
        end
        
        obj.labelfilename = labelfilename;
        [success,msg] = obj.UpdateStatusTable('label');   
      end
      
    end

    function [success,msg] = SetGTLabelFileName(obj,gt_labelfilename)
    % [success,msg] = SetGTLabelFileName(obj,labelfilename)
    % set the name of the *ground truth* label file within the experiment directory. this
    % does not currently update labelidx, and probably should not be called
    % once an experiment is open. 

      success = false;
      msg = '';

      if ischar(gt_labelfilename),
        if ischar(obj.gt_labelfilename) && strcmp(gt_labelfilename,obj.gt_labelfilename),
          success = true;
          return;
        end

        % reload labels from file
        for expi = 1:obj.nexps,
          [success1,msg] = obj.LoadLabelsFromFile(expi);
          if ~success1,
            return;
          end
        end
        
        obj.gt_labelfilename = gt_labelfilename;
        [success,msg] = obj.UpdateStatusTable('gt_label');   
      end
      
    end
    
    function [success,msg] = SetScoresFileName(obj,scorefilename)
    % [success,msg] = SetGTLabelFileName(obj,labelfilename)
    % set the name of the *ground truth* label file within the experiment directory. this
    % does not currently update labelidx, and probably should not be called
    % once an experiment is open. 

      success = false;
      msg = '';

      obj.scorefilename = scorefilename;
      [success,msg] = obj.UpdateStatusTable('scores');
      
    end

    function [success,msg] = SetClassifierType(obj,classifiertype)

      success = true;
      msg = '';
      
      % TODO: retrain classifier if necessary
      if strcmpi(classifiertype,obj.classifiertype),
        return;
      end
      
      obj.classifiertype = classifiertype;
      
    end
    
    function [success,msg] = LoadLabelsFromExternalFile(obj,expdir,labelfilename)
      
      success = false; msg = '';
      if ~exist(labelfilename,'file'),
        msg = sprintf('Label file %s does not exist',labelfilename);
        return;
      end
      expi = find(strcmp(expdir,obj.expdirs),1);
      if isempty(expi),
        msg = sprintf('Experiment %s not loaded in',expdir);
        return;
      end

      obj.SetStatus('Loading labels for %s from %s',obj.expdirs{expi},labelfilename);
      
      %         try
      loadedlabels = load(labelfilename,'t0s','t1s','names','flies','off','timestamp');
      
      if obj.IsGTMode(),
        obj.gt_labels(expi).t0s = loadedlabels.t0s;
        obj.gt_labels(expi).t1s = loadedlabels.t1s;
        obj.gt_labels(expi).names = loadedlabels.names;
        obj.gt_labels(expi).flies = loadedlabels.flies;
        obj.gt_labels(expi).off = loadedlabels.off;
        obj.gt_labelstats(expi).nflies_labeled = size(loadedlabels.flies,1);
        obj.gt_labelstats(expi).nbouts_labeled = numel([loadedlabels.t0s{:}]);
        
        if iscell(loadedlabels.timestamp)
          obj.gt_labels(expi).timestamp = loadedlabels.timestamp;
        else
          for ndx = 1:numel(loadedlabels.flies)
            nBouts = numel(loadedlabels.t0s{ndx});
            if isempty(loadedlabels.timestamp)
              obj.gt_labels(expi).timestamp{ndx}(1:nBouts) = now;
            else
              obj.gt_labels(expi).timestamp{ndx}(1:nBouts) = loadedlabels.timestamp;
            end
          end
        end
        
        if ~isempty(whos('-file',labelfilename,'imp_t0s'))
          loadedimp = load(labelfilename,'imp_t0s','imp_t1s');
          obj.gt_labels(expi).imp_t0s = loadedimp.imp_t0s;
          obj.gt_labels(expi).imp_t1s = loadedimp.imp_t1s;
        else
          obj.gt_labels(expi).imp_t0s = cell(1,numel(loadedlabels.flies));
          obj.gt_labels(expi).imp_t1s = cell(1,numel(loadedlabels.flies));
        end
        
      else
        
        obj.labels(expi).t0s = loadedlabels.t0s;
        obj.labels(expi).t1s = loadedlabels.t1s;
        obj.labels(expi).names = loadedlabels.names;
        obj.labels(expi).flies = loadedlabels.flies;
        obj.labels(expi).off = loadedlabels.off;
        obj.labelstats(expi).nflies_labeled = size(loadedlabels.flies,1);
        obj.labelstats(expi).nbouts_labeled = numel([loadedlabels.t0s{:}]);
        if iscell(loadedlabels.timestamp)
          obj.labels(expi).timestamp = loadedlabels.timestamp;
        else
          for ndx = 1:numel(loadedlabels.flies)
            nBouts = numel(loadedlabels.t0s{ndx});
            if isempty(loadedlabels.timestamp)
              obj.labels(expi).timestamp{ndx}(1:nBouts) = now;
            else
              obj.labels(expi).timestamp{ndx}(1:nBouts) = loadedlabels.timestamp;
            end
          end
        end
        if ~isempty(whos('-file',labelfilename,'imp_t0s'))
          loadedimp = load(labelfilename,'imp_t0s','imp_t1s');
          obj.labels(expi).imp_t0s = loadedimp.imp_t0s;
          obj.labels(expi).imp_t1s = loadedimp.imp_t1s;
        else
          obj.labels(expi).imp_t0s = cell(1,numel(loadedlabels.flies));
          obj.labels(expi).imp_t1s = cell(1,numel(loadedlabels.flies));
        end
        
      end
        %         catch ME,
        %           msg = getReport(ME);
        %           obj.ClearStatus();
        %           return;
        %         end
      
      obj.ClearStatus();
      success = true;
      
    end
    
    function [success,msg] = LoadLabelsFromFile(obj,expi)
    % [success,msg] = LoadLabelsFromFile(obj,expi)
    % If the label file exists, this function loads labels for experiment
    % expi into obj.labels. Otherwise, it sets the labels to be empty. This
    % does not currently update the windowdata and labelidx (TODO). 
      
      success = false; msg = '';
      labelfilename = obj.GetFile('label',expi);
      
      if exist(labelfilename,'file'),

        obj.SetStatus('Loading labels for %s',obj.expdirs{expi});
        
        %         try
        loadedlabels = load(labelfilename,'t0s','t1s','names','flies','off','timestamp');
        obj.labels(expi).t0s = loadedlabels.t0s;
        obj.labels(expi).t1s = loadedlabels.t1s;
        obj.labels(expi).names = loadedlabels.names;
        obj.labels(expi).flies = loadedlabels.flies;
        obj.labels(expi).off = loadedlabels.off;
        obj.labelstats(expi).nflies_labeled = size(loadedlabels.flies,1);
        obj.labelstats(expi).nbouts_labeled = numel([loadedlabels.t0s{:}]);
        if iscell(loadedlabels.timestamp)
          obj.labels(expi).timestamp = loadedlabels.timestamp;
        else
          for ndx = 1:numel(loadedlabels.flies)
            nBouts = numel(loadedlabels.t0s{ndx});
            if isempty(loadedlabels.timestamp)
              obj.labels(expi).timestamp{ndx}(1:nBouts) = now;
            else
              obj.labels(expi).timestamp{ndx}(1:nBouts) = loadedlabels.timestamp;
            end
          end
        end
        if ~isempty(whos('-file',labelfilename,'imp_t0s'))
          loadedimp = load(labelfilename,'imp_t0s','imp_t1s');
          obj.labels(expi).imp_t0s = loadedimp.imp_t0s;
          obj.labels(expi).imp_t1s = loadedimp.imp_t1s;
        else
          obj.labels(expi).imp_t0s = cell(1,numel(loadedlabels.flies));
          obj.labels(expi).imp_t1s = cell(1,numel(loadedlabels.flies));
        end
        %         catch ME,
        %           msg = getReport(ME);
        %           obj.ClearStatus();
        %           return;
        %         end
        
        obj.ClearStatus();
        
      else
        
        obj.labels(expi).t0s = {};
        obj.labels(expi).t1s = {};
        obj.labels(expi).names = {};
        obj.labels(expi).flies = [];
        obj.labels(expi).off = [];
        obj.labels(expi).timestamp = {};
        obj.labels(expi).imp_t0s = {};
        obj.labels(expi).imp_t1s = {};
        obj.labelstats(expi).nflies_labeled = 0;
        obj.labelstats(expi).nbouts_labeled = 0;

      end
      
      % TODO: update windowdata
      
      success = true;
      
    end
 
    function [success,msg] = LoadGTLabelsFromFile(obj,expi)
    % [success,msg] = LoadGTLabelsFromFile(obj,expi)
    % If the label file exists, this function loads labels for experiment
    % expi into obj.gt_labels. Otherwise, it sets the gt_labels to be empty. This
    % does not currently update the windowdata and labelidx (TODO). 
      
      success = false; msg = '';
      
      labelfilename = obj.GetFile('gt_label',expi);

      if exist(labelfilename,'file'),
        
        obj.SetStatus('Loading labels for %s',obj.expdirs{expi});
        
        %         try
        loadedlabels = load(labelfilename,'t0s','t1s','names','flies','off','timestamp');
        obj.gt_labels(expi).t0s = loadedlabels.t0s;
        obj.gt_labels(expi).t1s = loadedlabels.t1s;
        obj.gt_labels(expi).names = loadedlabels.names;
        obj.gt_labels(expi).flies = loadedlabels.flies;
        obj.gt_labels(expi).off = loadedlabels.off;
        obj.gt_labelstats(expi).nflies_labeled = size(loadedlabels.flies,1);
        obj.gt_labelstats(expi).nbouts_labeled = numel([loadedlabels.t0s{:}]);

        if iscell(loadedlabels.timestamp)
          obj.gt_labels(expi).timestamp = loadedlabels.timestamp;
        else
          for ndx = 1:numel(loadedlabels.flies)
            nBouts = numel(loadedlabels.t0s{ndx});
            if isempty(loadedlabels.timestamp)
              obj.gt_labels(expi).timestamp{ndx}(1:nBouts) = now;
            else
              obj.gt_labels(expi).timestamp{ndx}(1:nBouts) = loadedlabels.timestamp;
            end
          end
        end
        
        if ~isempty(whos('-file',labelfilename,'imp_t0s'))
          loadedimp = load(labelfilename,'imp_t0s','imp_t1s');
          obj.gt_labels(expi).imp_t0s = loadedimp.imp_t0s;
          obj.gt_labels(expi).imp_t1s = loadedimp.imp_t1s;
        else
          obj.gt_labels(expi).imp_t0s = cell(1,numel(loadedlabels.flies));
          obj.gt_labels(expi).imp_t1s = cell(1,numel(loadedlabels.flies));
        end
        %         catch ME,
        %           msg = getReport(ME);
        %           obj.ClearStatus();
        %           return;
        %         end
        
        obj.ClearStatus();
        
      else
        
        obj.gt_labels(expi).t0s = {};
        obj.gt_labels(expi).t1s = {};
        obj.gt_labels(expi).names = {};
        obj.gt_labels(expi).flies = [];
        obj.gt_labels(expi).off = [];
        obj.gt_labels(expi).timestamp = {};
        obj.gt_labels(expi).imp_t0s = {};
        obj.gt_labels(expi).imp_t1s = {};
        obj.gt_labelstats(expi).nflies_labeled = 0;
        obj.gt_labelstats(expi).nbouts_labeled = 0;
      end
      
      % TODO: update windowdata
      
      success = true;
      
    end
    
    function [success,msg] = SetPerFrameDir(obj,perframedir)
      % [success,msg] = SetPerFrameDir(obj,perframedir)
      % Sets the per-frame directory name within the experiment directory.
      % Currently, this does not change the cached per-frame data or check
      % that all the per-frame files necessary are within the directory
      % (TODO).
      
      success = false; msg = '';
      
      if ischar(perframedir),
        if ischar(obj.perframedir) && strcmp(perframedir,obj.perframedir),
          success = true;
          return;
        end
        
        obj.perframedir = perframedir;
        
        % TODO: check per-frame directories are okay, remove bad
        % experiments
        
        [success,msg] = obj.UpdateStatusTable('perframedir');
      end
      
    end
    
    function [success,msg] = SetClipsDir(obj,clipsdir)
    % [success,msg] = SetClipsDir(obj,clipsdir)
    % Sets the clips directory name within the experiment directory.
      
      success = false;
      msg = '';

      if ischar(clipsdir),
        for i = 1:numel(obj.expdirs),
          clipsdircurr = fullfile(obj.expdirs{i},clipsdir);
%           if exist(obj.expdirs{i},'dir') && ~exist(clipsdircurr,'dir'),
%             mkdir(clipsdircurr);
%           end
        end
        if ischar(obj.clipsdir) && strcmp(clipsdir,obj.clipsdir),
          success = true;
          return;
        end

        obj.clipsdir = clipsdir;        
        [success,msg] = obj.UpdateStatusTable('clipsdir');
      end
      
    end

    function [success,msg] = SetDefaultPath(obj,defaultpath)
    % [success,msg] = SetDefaultPath(obj,defaultpath)
    % sets the default path to load experiments from. only checks for
    % existence of the directory.
      
      success = false;
      msg = '';
      
      if ischar(defaultpath),
        
        if ~isempty(defaultpath) && ~exist(defaultpath,'file'),
          msg = sprintf('defaultpath directory %s does not exist',defaultpath);
          return;
        end
          
        obj.defaultpath = defaultpath;
        success = true;
      end

    end
    
    function [success,msg] = SetRootOutputDir(obj,rootoutputdir)
    % [success,msg] = SetRootOutputDir(obj,rootoutputdir)
    % sets the root directory for outputing files. currently, it does not
    % update labels, etc. or recheck for the existence of all the required
    % files. (TODO)
      
      success = true;
      msg = '';
      if ischar(rootoutputdir),
        if ischar(obj.rootoutputdir) && strcmp(obj.rootoutputdir,rootoutputdir),
          success = true;
          return;
        end
        if ~exist(rootoutputdir,'file'),
          msg = sprintf('root output directory %s does not exist, outputs will be stored in the experiment directories',...
            rootoutputdir);
          success = false;
          return;
        end
        obj.rootoutputdir = rootoutputdir;
        for i = 1:obj.nexps,
          obj.outexpdirs{i} = fullfile(rootoutputdir,obj.expnames{i});
        end
        % TODO: check all files are okay, remove bad experiments
        
        [success,msg] = obj.UpdateStatusTable();
      end
      
    end    
    
    function [success,msg] = SetClassifierFileName(obj,classifierfilename,varargin)
    % [success,msg] = SetClassifierFileName(obj,classifierfilename)
    % Sets the name of the classifier file. If the classifier file exists, 
    % it loads the data stored in the file. This involves removing all the
    % experiments and data currently loaded, setting the config file,
    % setting all the file names set in the config file, setting the
    % experiments to be those listed in the classifier file, clearing all
    % the previously computed window data and computing the window data for
    % all the labeled frames. 
      
    [classifierlabels,doreadconfigfile] = myparse(varargin,...
      'classifierlabels',false,...
      'doreadconfigfile',true);
    
      success = false;
      msg = '';
      
      obj.classifierfilename = classifierfilename;
      if ~isempty(classifierfilename) && exist(classifierfilename,'file'),
%         try

          loadeddata = load(obj.classifierfilename); %,obj.classifiervars{:});
          
          if ~strcmp(loadeddata.labelfilename,obj.labelfilename),
            success = false;
            msg = ['Label files specified for the project doesn''t match' ...
              ' the labelfiles used to train the classifier. Not loading the classifier'];
            return;
          end
          
          obj.SetStatus('Loading classifier from %s',obj.classifierfilename);

          % remove all experiments
          obj.RemoveExpDirs(1:obj.nexps);
          
          if doreadconfigfile,
          
          % set config file
          if ~strcmp(obj.configfilename,'configfilename'),
            obj.SetConfigFileName(loadeddata.configfilename);
          end

          % set movie
          [success,msg] = obj.SetMovieFileName(loadeddata.moviefilename);
          if ~success,error(msg);end

          % trx
          [success,msg] = obj.SetTrxFileName(loadeddata.trxfilename);
          if ~success,error(msg);end
      
          % label
          [success,msg] = obj.SetLabelFileName(loadeddata.labelfilename);
          if ~success,error(msg);end
      
          % perframedir
          [success,msg] = obj.SetPerFrameDir(loadeddata.perframedir);
          if ~success,error(msg);end

          % clipsdir
          [success,msg] = obj.SetClipsDir(loadeddata.clipsdir);
          if ~success,error(msg);end
          
          end
          
          % featureparamsfilename
%           [success,msg] = obj.SetFeatureParamsFileName(loadeddata.featureparamsfilename);
%           if ~success,error(msg);end
          % load actual window features params instead of filename.
          if all( isfield(loadeddata,{'windowfeaturesparams','windowfeaturescellparams',...
              'basicFeatureTable','featureWindowSize'}))
            
            loadeddata.windowfeaturesparams = JLabelData.convertTransTypes2Cell(loadeddata.windowfeaturesparams);
            if ~( isequal(obj.windowfeaturesparams,loadeddata.windowfeaturesparams) && ...
                  isequal(obj.featureWindowSize,loadeddata.featureWindowSize)),
                str = sprintf('Window feature parameters in the configuration file');
                str = sprintf('%s\ndo not match the parameters saved in the classifier',str);
                str = sprintf('%s\nUsing parameters stored in the classifier file',str);
                uiwait(warndlg(str));
            end
            obj.UpdatePerframeParams(loadeddata.windowfeaturesparams,...
              loadeddata.windowfeaturescellparams,loadeddata.basicFeatureTable,...
              loadeddata.featureWindowSize);
          end
          
          if ~isfield(loadeddata,'featurenames')
            feature_names = {};
            for j = 1:numel(obj.curperframefns),
              fn = obj.curperframefns{j};
              [~,feature_names_curr] = ComputeWindowFeatures([0,0],...
                obj.windowfeaturescellparams.(fn){:});
              feature_names_curr = cellfun(@(x) [{fn},x],feature_names_curr,'UniformOutput',false);
              feature_names = [feature_names,feature_names_curr]; %#ok<AGROW>
            end
            obj.windowdata.featurenames = feature_names;
          else
            obj.windowdata.featurenames = loadeddata.featurenames;
          end
          
      
          % rootoutputdir
%           [success,msg] = obj.SetRootOutputDir(loadeddata.rootoutputdir);
%           if ~success,error(msg); end
           
          % set experiment directories
          if classifierlabels && isfield(loadeddata,'labels'),
            [success,msg] = obj.SetExpDirs(loadeddata.expdirs,loadeddata.outexpdirs,...
              loadeddata.nflies_per_exp,loadeddata.sex_per_exp,loadeddata.frac_sex_per_exp,...
              loadeddata.firstframes_per_exp,loadeddata.endframes_per_exp);
            if ~success,error(msg); end
            obj.labels = loadeddata.labels;
            [obj.labelidx,obj.t0_curr,obj.t1_curr] = obj.GetLabelIdx(expi,flies);
            obj.labelidx_off = 1 - obj.t0_curr;
            [success,msg] = obj.PreLoadLabeledData();
            if ~success,error(msg); end
            obj.labelsLoadedFromClassifier = true;
            
          else
            if classifierlabels,
              uiwait(warndlg('The classifier file didn''t have any labels. Loading the current labels'));
            end
            [success,msg] = obj.SetExpDirs(loadeddata.expdirs,loadeddata.outexpdirs,...
              loadeddata.nflies_per_exp,loadeddata.sex_per_exp,loadeddata.frac_sex_per_exp,...
              loadeddata.firstframes_per_exp,loadeddata.endframes_per_exp);
            if ~success,error(msg); end
          end
          [success,msg] = obj.UpdateStatusTable();
          if ~success, error(msg); end
          
          % update cached data
%           obj.windowdata = struct('X',[],'exp',[],'flies',[],'t',[],...
%             'labelidx_cur',[],'labelidx_new',[],'featurenames',{{}},...
%             'predicted',[],'predicted_probs',[],'isvalidprediction',[]);
          [success,msg] = obj.PreLoadLabeledData();
          if ~success,error(msg);end
                                       
          obj.classifier = loadeddata.classifier;
          obj.classifiertype = loadeddata.classifiertype;
          obj.classifierTS = loadeddata.classifierTS;
          obj.windowdata.scoreNorm = loadeddata.scoreNorm;
          obj.confThresholds = loadeddata.confThresholds;
          if isfield(loadeddata,'postprocessparams')
            obj.postprocessparams = loadeddata.postprocessparams;
          end
          
          paramFields = fieldnames(loadeddata.classifier_params);
          for ndx = 1:numel(paramFields)
            obj.classifier_params.(paramFields{ndx}) = loadeddata.classifier_params.(paramFields{ndx});
          end
          % predict for all loaded examples
          obj.PredictLoaded();
          
          % set labelidx_cur
          obj.SetTrainingData(loadeddata.trainingdata);

%           if strcmp(obj.classifiertype,'boosting'),
%             [obj.windowdata.binVals, obj.windowdata.bins] = findThresholds(obj.windowdata.X);
%           end
          
          % make sure inds is ordered correctly
          if ~isempty(obj.classifier),
            switch obj.classifiertype,
              
              case 'ferns',
                waslabeled = obj.windowdata.labelidx_cur ~= 0;
                obj.classifier.inds = obj.predict_cache.last_predicted_inds(waslabeled,:);
            
            end
          end
          
          % clear the cached per-frame, trx data
          obj.ClearCachedPerExpData();
          
%         catch ME,
%           errordlg(getReport(ME),'Error loading classifier from file');
%         end
        
        obj.ClearStatus();
        
        obj.classifierfilename = classifierfilename;
 
        obj.FindFastPredictParams();
 
      end

    end
    
    function [success,msg] = SetClassifierFileNameWoExp(obj,classifierfilename)

      success = false;
      msg = '';
      
      obj.classifierfilename = classifierfilename;
      if ~isempty(classifierfilename) && exist(classifierfilename,'file'),
%         try

          loadeddata = load(obj.classifierfilename); %,obj.classifiervars{:});

          if ~strcmp(loadeddata.labelfilename,obj.labelfilename),
            success = false;
            msg = ['Label files specified for the project doesn''t match the labelfiles '...
              'used to train the classifier. Not loading the classifier'];
            return;
          end
          
          obj.SetStatus('Loading classifier from %s',obj.classifierfilename);

          % remove all experiments
          % obj.RemoveExpDirs(1:obj.nexps);

          % set config file
%           if ~strcmp(obj.configfilename,'configfilename'),
%             obj.SetConfigFileName(loadeddata.configfilename);
%           end

          % set movie
          [success,msg] = obj.SetMovieFileName(loadeddata.moviefilename);
          if ~success,error(msg);end

          % trx
          [success,msg] = obj.SetTrxFileName(loadeddata.trxfilename);
          if ~success,error(msg);end
      
          % label
          [success,msg] = obj.SetLabelFileName(loadeddata.labelfilename);
          if ~success,error(msg);end
      
          % perframedir
          [success,msg] = obj.SetPerFrameDir(loadeddata.perframedir);
          if ~success,error(msg);end

          % clipsdir
          [success,msg] = obj.SetClipsDir(loadeddata.clipsdir);
          if ~success,error(msg);end

          % featureparamsfilename
%           [success,msg] = obj.SetFeatureParamsFileName(loadeddata.featureparamsfilename);
%           if ~success,error(msg);end

          % load actual window features params instead of filename.
          if all( isfield(loadeddata,{'windowfeaturesparams','windowfeaturescellparams',...
              'basicFeatureTable','featureWindowSize'}))
            if ~( isequal(obj.windowfeaturesparams,loadeddata.windowfeaturesparams) && ...
                  isequal(obj.windowfeaturescellparams,loadeddata.windowfeaturescellparams) && ...
                  isequal(obj.featureWindowSize,loadeddata.featureWindowSize)),
                str = sprintf('Window feature parameters in the configuration file');
                str = sprintf('%s\ndo not match the parameters saved in the classifier',str);
                str = sprintf('%s\nUsing parameters stored in the classifier file',str);
                uiwait(warndlg(str));
                obj.UpdatePerframeParams(loadeddata.windowfeaturesparams,...
                  loadeddata.windowfeaturescellparams,loadeddata.basicFeatureTable,...
                  loadeddata.featureWindowSize);
            end
          end
          
          if ~isfield(loadeddata,'featurenames')
            feature_names = {};
            for j = 1:numel(obj.curperframefns),
              fn = obj.curperframefns{j};
              [~,feature_names_curr] = ComputeWindowFeatures([0,0],...
                obj.windowfeaturescellparams.(fn){:});
              feature_names_curr = cellfun(@(x) [{fn},x],feature_names_curr,'UniformOutput',false);
              feature_names = [feature_names,feature_names_curr]; %#ok<AGROW>
            end
            obj.windowdata.featurenames = feature_names;
          else
            obj.windowdata.featurenames = loadeddata.featurenames;
          end
          
          % rootoutputdir
%           [success,msg] = obj.SetRootOutputDir(loadeddata.rootoutputdir);
%           if ~success,error(msg); end

          [success,msg] = obj.UpdateStatusTable();
          if ~success, error(msg); end
                         

          obj.classifier = loadeddata.classifier;
          obj.classifiertype = loadeddata.classifiertype;
          obj.classifierTS = loadeddata.classifierTS;
          obj.classifier_params = loadeddata.classifier_params;
          obj.windowdata.scoreNorm = loadeddata.scoreNorm;
          obj.confThresholds = loadeddata.confThresholds;
          obj.classifierfilename = classifierfilename;
          paramFields = fieldnames(loadeddata.classifier_params);
          for ndx = 1:numel(paramFields)
            obj.classifier_params.(paramFields{ndx}) = loadeddata.classifier_params.(paramFields{ndx});
          end
          % obj.ClearCachedPerExpData();
          obj.ClearStatus();
          obj.FindFastPredictParams();
 
      end
    end
    
    function [success,msg] = SetExpDirs(obj,expdirs,outexpdirs,nflies_per_exp,...
        sex_per_exp,frac_sex_per_exp,firstframes_per_exp,endframes_per_exp)
    % [success,msg] = SetExpDirs(obj,[expdirs,outexpdirs,nflies_per_exp,firstframes_per_exp,endframes_per_exp])
    % Changes what experiments are currently being used for this
    % classifier. This function calls RemoveExpDirs to remove all current
    % experiments not in expdirs, then calls AddExpDirs to add the new
    % experiment directories. 

      success = false;
      msg = '';
      
      if isnumeric(expdirs),
        return;
      end
      
      if nargin < 2,
        error('Usage: obj.SetExpDirs(expdirs,[outexpdirs],[nflies_per_exp])');
      end
      
      isoutexpdirs = nargin > 2 && ~isnumeric(outexpdirs);
      isnflies = nargin > 3 && ~isempty(nflies_per_exp);
      issex = nargin > 4 && ~isempty(sex_per_exp);
      isfracsex = nargin > 5 && ~isempty(frac_sex_per_exp);
      isfirstframes = nargin > 6 && ~isempty(firstframes_per_exp);
      isendframes = nargin > 7 && ~isempty(endframes_per_exp);
      
      % check inputs
      
      % sizes must match
      if isoutexpdirs && numel(expdirs) ~= numel(outexpdirs),
        error('expdirs and outexpdirs do not match size');
      end
      if isnflies && numel(expdirs) ~= numel(nflies_per_exp),
        error('expdirs and nflies_per_exp do not match size');
      end
      
      oldexpdirs = obj.expdirs;
      
      % remove oldexpdirs
      
      [success1,msg] = obj.RemoveExpDirs(find(~ismember(oldexpdirs,expdirs))); %#ok<FNDSB>
      if ~success1,
        return;
      end

      % add new expdirs
      idx = find(~ismember(expdirs,oldexpdirs));
      success = true;
      for i = idx,
        params = cell(1,nargin-1);
        params{1} = expdirs{i};
        if isoutexpdirs,
          params{2} = outexpdirs{i};
        end
        if isnflies,
          params{3} = nflies_per_exp(i);
        end
        if issex,
          params{4} = sex_per_exp{i};
        end
        if isfracsex,
          params{5} = frac_sex_per_exp{i};
        end
        if isfirstframes,
          params{6} = firstframes_per_exp{i};
        end
        if isendframes,
          params{7} = endframes_per_exp{i};
        end
        [success1,msg1] = obj.AddExpDir(params{:});
        success = success && success1;
        if isempty(msg),
          msg = msg1;
        else
          msg = sprintf('%s\n%s',msg,msg1);
        end        
      end

    end

    
% Saving and loading    


    function SaveScores(obj,allScores,expi,sfn)
    % Save prediction scores for the whole experiment.
    % The scores are stored as a cell array.
     if nargin< 4
      sfn = obj.GetFile('scores',expi,true);
     end
      obj.SetStatus('Saving scores for experiment %s to %s',obj.expnames{expi},sfn);

      didbak = false;
      if exist(sfn,'file'),
        [didbak,msg] = copyfile(sfn,[sfn,'~']);
        if ~didbak,
          warning('Could not create backup of %s: %s',sfn,msg);
        end
      end
      timestamp = obj.classifierTS;
      save(sfn,'allScores','timestamp');
      obj.ClearStatus();
    end
    
    function AddScores(obj,expi,allScores,timestamp,classifierfilename,updateCurrent)
%       obj.predictdata.classifierfilenames{expi} = classifierfilename;
      obj.SetStatus('Updating Predictions ...');
      for ndx = 1:numel(allScores.scores)
        idxcurr = obj.FlyNdxPredict(expi,ndx);
        
        tStart = allScores.tStart(ndx);
        tEnd = allScores.tEnd(ndx);
        sz = tEnd-tStart+1;
        if (nnz(idxcurr) ~= sz),
          uiwait(warndlg(['Cannot load scores for experiment:%s. Number of frames '...
            'in the score file dont match the number of frames in the trax file for '...
            'fly:%d'],obj.expnames{expi},ndx));
          return;
        end
        curScores = allScores.scores{ndx}(tStart:tEnd);
        if updateCurrent,
          obj.predictdata.cur(idxcurr) = curScores;
          obj.predictdata.cur_valid(idxcurr) = true;
        else
          obj.predictdata.loaded(idxcurr) = curScores;
          obj.predictdata.loaded_valid(idxcurr) = true;
        end
      end
      

      if isempty(obj.windowdata.scoreNorm) || isnan(obj.windowdata.scoreNorm)
        if ~isempty(obj.predictdata.loaded)
          ss = obj.predictdata.loaded;
          ss = ss(~isnan(ss));
          scoreNorm = prctile(abs(ss),80);
          obj.windowdata.scoreNorm = scoreNorm;
        end
      end
      
      if ~isempty(obj.postprocessparams),
        [success,msg] = obj.ApplyPostprocessing();
        if ~success,
          uiwait(warndlg(['Couldn''t apply postprocessing to the scores: ' msg]));
        end
      elseif ~updateCurrent,
        for ndx = 1:numel(allScores.loaded)
          idxcurr = obj.FlyNdxPredict(expi,ndx);
          tStart = allScores.tStart(ndx);
          tEnd = allScores.tEnd(ndx);
          if isfield(allScores,'postprocessedscores');
            obj.postprocessparams = allScores.postprocessparams;
            curpostprocessedscores = allScores.postprocessedscores{ndx}(tStart:tEnd);
            obj.predictdata.loaded_pp(idxcurr) = curpostprocessedscores;
          else
            obj.predictdata.loaded_pp(idxcurr) = 0;
          end
        end
      end
      
      obj.UpdatePredictedIdx();
      obj.ClearStatus();

    end
    
    function SaveCurScores(obj,expi,sfn)
    % Saves the current scores to a file.
      if nargin < 3
        sfn = obj.GetFile('scores',expi,true);
      end
    
      if isempty(obj.predictdata.exp),
        uiwait(warndlg('No scores to save'));
        return;
      end
      
      allScores = struct('scores',{{}},'tStart',[],'tEnd',[],...
        'postprocessed',{{}},'postprocessedparams',[]);
      scores_valid = true;
      for fly = 1:obj.nflies_per_exp(expi)
        idxcurr = obj.FlyNdxPredict(expi,fly);
        
        curt = obj.predictdata.t(idxcurr);
        if any(curt(2:end)-curt(1:end-1) ~= 1)
          uiwait(warndlg('Scores are out of order. This shouldn''t happen. Not saving them'));
          return;
        end
        
        if ~all(obj.predictdata.cur(idxcurr)), 
          scores_valid = false; 
          break; 
        end
        
        tStart = obj.firstframes_per_exp{expi}(fly);
        tEnd = obj.endframes_per_exp{expi}(fly);
        
        allScores.scores{fly}(tStart:tEnd) = obj.predictdata.cur(idxcurr);
        allScores.tStart(fly) = tStart;
        allScores.tEnd(fly) = tEnd;
        allScores.postprocessed{fly}(tStart:tEnd) = obj.predictdata.cur_pp(idxcurr);
      end
      
      if ~scores_valid,
        uiwait(warndlg(['Scores have not been computed for all the frames for experiment ' ...
         '%s. Cannot save the scores.'],obj.expnames{expi}));
        return;
      end
      allScores.postprocessedparams = obj.postprocessparams;
      obj.SaveScores(allScores,expi,sfn);
      
    end
    
    function LoadScores(obj,expi,sfn)
      
      obj.SetStatus('Loading scores for experiment %s from %s',obj.expnames{expi},sfn);
      if ~exist(sfn,'file')
        warndlg('Score file %s does not exist. Not loading scores',sfn);
        return;
      end
      load(sfn,'allScores','timestamp');
      if ~isempty(whos('-file',sfn,'classifierfilename'))
        S = load(sfn,'classifierfilename');
        classifierfilename = S.classifierfilename;
      else
        classifierfilename = '';
      end
      
      obj.AddScores(expi,allScores,timestamp,classifierfilename,false);
      
      obj.ClearStatus();

    end
    
    function LoadScoresDefault(obj,expi)
      sfn = obj.GetFile('scores',expi);
      if ~exist(sfn,'file')
        warndlg(sprintf('No scores file %s at the default location',...
          sfn));
        return;
      end
      obj.LoadScores(expi,sfn);
    end
    
    function SaveClassifier(obj)
    % SaveClassifier(obj)
    % This function saves the current classifier to the file
    % ons.classifierfilename. It first constructs a struct representing the
    % training data last used to train the classifier, then adds all the
    % data described in obj.classifiervars.       
      
      s = struct;
      s.classifierTS = obj.classifierTS;
      s.trainingdata = obj.SummarizeTrainingData();
%       try
        for i = 1:numel(obj.classifiervars),
          fn = obj.classifiervars{i};
          if isfield(s,fn),
          % elseif isprop(obj,fn),
          % isprop doesn't work right on 2010b
          elseif ismember(fn,properties(obj))
            s.(fn) = obj.(fn);
          elseif isstruct(obj.windowdata) && isfield(obj.windowdata,fn),
            s.(fn) = obj.windowdata.(fn);
          else
            error('Unknown field %s',fn);
          end
            
        end
        save(obj.classifierfilename,'-struct','s');
%       catch ME,
%         errordlg(getReport(ME),'Error saving classifier to file');
%       end      
      
    end

    function SaveLabels(obj,expis)
    % SaveLabels(obj,expis)
    % For each experiment in expis, save the current set of labels to file.
    % A backup of old labels is made if they exist and stored in
    % <labelfilename>~
    
      if isempty(obj.labels), return; end
    
      if nargin<2
        expis = 1:obj.nexps;
      end
      
      if obj.labelsLoadedFromClassifier,
        res = questdlg(['Labels were loaded from the classifier. Saving the'...
          ' labels will overwrite the current labels. Overwrite?'],...
          'Overwrite Current Labels?','Yes','No','Cancel','No');
        if ~strcmpi(res,'Yes'), return, end
      end
      
      % store labels in labelidx
      obj.StoreLabels();
      
      for i = expis,
        
        
        lfn = GetFile(obj,'label',i,true);
        obj.SetStatus('Saving labels for experiment %s to %s',obj.expnames{i},lfn);

        didbak = false;
        if exist(lfn,'file'),
          [didbak,msg] = copyfile(lfn,[lfn,'~']);
          if ~didbak,
            warning('Could not create backup of %s: %s',lfn,msg);
          end
        end

        t0s = obj.labels(i).t0s; %#ok<NASGU>
        t1s = obj.labels(i).t1s; %#ok<NASGU>
        names = obj.labels(i).names; %#ok<NASGU>
        flies = obj.labels(i).flies; %#ok<NASGU>
        off = obj.labels(i).off; %#ok<NASGU>
        timestamp = obj.labels(i).timestamp; %#ok<NASGU>
        imp_t0s = obj.labels(i).imp_t0s; %#ok<NASGU>
        imp_t1s = obj.labels(i).imp_t1s; %#ok<NASGU>
        
%         try
          save(lfn,'t0s','t1s','names','flies','off','timestamp','imp_t0s','imp_t1s');
%         catch ME,
%           if didbak,
%             [didundo,msg] = copyfile([lfn,'~'],lfn);
%             if ~didundo, warning('Error copying backup file for %s: %s',lfn,msg); end
%           end
%           errordlg(sprintf('Error saving label file %s: %s.',lfn,getReport(ME)),'Error saving labels');
%         end
      end

      
      [success,msg] = obj.UpdateStatusTable('label');
      if ~success,
        error(msg);
      end

      obj.ClearStatus();
    end

    function SaveGTLabels(obj,expis)
    % SaveGTLabels(obj,expis)
    % For each experiment in expis, save the current set of labels to file.
    % A backup of old labels is made if they exist and stored in
    % <labelfilename>~

      if isempty(obj.gt_labels), return; end
      
      if nargin<2
        expis = 1:obj.nexps;
      end
      
      % store labels in labelidx
      obj.StoreLabels();
      
      for i = expis,

        lfn = GetFile(obj,'gt_label',i,true);
        obj.SetStatus('Saving labels for experiment %s to %s',obj.expnames{i},lfn);

        didbak = false;
        if exist(lfn,'file'),
          [didbak,msg] = copyfile(lfn,[lfn,'~']);
          if ~didbak,
            warning('Could not create backup of %s: %s',lfn,msg);
          end
        end

        t0s = obj.gt_labels(i).t0s; %#ok<NASGU>
        t1s = obj.gt_labels(i).t1s; %#ok<NASGU>
        names = obj.gt_labels(i).names; %#ok<NASGU>
        flies = obj.gt_labels(i).flies; %#ok<NASGU>
        off = obj.gt_labels(i).off; %#ok<NASGU>
        timestamp = obj.gt_labels(i).timestamp; %#ok<NASGU>
        imp_t0s = obj.gt_labels(i).imp_t0s; %#ok<NASGU>
        imp_t1s = obj.gt_labels(i).imp_t1s; %#ok<NASGU>
        
%         try
          save(lfn,'t0s','t1s','names','flies','off','timestamp','imp_t0s','imp_t1s');
%         catch ME,
%           if didbak,
%             [didundo,msg] = copyfile([lfn,'~'],lfn);
%             if ~didundo, warning('Error copying backup file for %s: %s',lfn,msg); end
%           end
%           errordlg(sprintf('Error saving label file %s: %s.',lfn,getReport(ME)),'Error saving labels');
%         end
      end

      
      [success,msg] = obj.UpdateStatusTable('gt_label');
      if ~success,
        error(msg);
      end

      obj.ClearStatus();
    end

    
% Experiment handling


    function [success,msg] = AddExpDir(obj,expdir,outexpdir,nflies_per_exp,sex_per_exp,...
        frac_sex_per_exp,firstframes_per_exp,endframes_per_exp)
    % [success,msg] = AddExpDir(obj,expdir,outexpdir,nflies_per_exp,firstframes_per_exp,endframes_per_exp)
    % Add a new experiment to the GUI. If this is the first experiment,
    % then it will be preloaded. 

      success = false; msg = '';
      
      if isnumeric(expdir), return; end
      
      if nargin < 2,
        error('Usage: obj.AddExpDirs(expdir,[outexpdir],[nflies_per_exp])');
      end

      % make sure directory exists
      if ~exist(expdir,'file'),
        msg = sprintf('expdir %s does not exist',expdir);
        return;
      end
      
      isoutexpdir = nargin > 2 && ~isnumeric(outexpdir);
      istrxinfo = nargin > 7 && ~isempty(nflies_per_exp);

      % base name
      [~,expname] = myfileparts(expdir);
      
      % expnames and rootoutputdir must match
      if isoutexpdir,
        [rootoutputdir,outname] = myfileparts(outexpdir); %#ok<*PROP>
        if ~strcmp(expname,outname),
          msg = sprintf('expdir and outexpdir do not match base names: %s ~= %s',expname,outname);
          return;
        end
%         if ischar(obj.rootoutputdir) && ~strcmp(rootoutputdir,obj.rootoutputdir),
%           msg = sprintf('Inconsistent root output directory: %s ~= %s',rootoutputdir,obj.rootoutputdir);
%           return;
%         end
      elseif ~ischar(obj.rootoutputdir),
        outexpdir = expdir;
        rootoutputdir = 0;
      else
        rootoutputdir = obj.rootoutputdir;        
      end
      
      if ischar(obj.rootoutputdir) && ~isoutexpdir,
        outexpdir = fullfile(rootoutputdir,expname);
      end
      
      % create missing outexpdirs
      if ~exist(outexpdir,'dir'),
        [success1,msg1] = mkdir(rootoutputdir,expname);
        if ~success1,
          msg = (sprintf('Could not create output directory %s, failed to set expdirs: %s',outexpdir,msg1));
          return;
        end
      end

      % create clips dir
      clipsdir = obj.GetFileName('clipsdir');
      outclipsdir = fullfile(outexpdir,clipsdir);
%       if ~exist(outclipsdir,'dir'),
%         [success1,msg1] = mkdir(outexpdir,clipsdir);
%         if ~success1,
%           msg = (sprintf('Could not create output clip directory %s, failed to set expdirs: %s',outclipsdir,msg1));
%           return;
%         end
%       end

      % okay, checks succeeded, start storing stuff
      obj.nexps = obj.nexps + 1;
      obj.expdirs{end+1} = expdir;
      obj.expnames{end+1} = expname;
      %obj.rootoutputdir = rootoutputdir;
      obj.outexpdirs{end+1} = outexpdir;
      
      % load labels for this experiment
      [success1,msg] = obj.LoadLabelsFromFile(obj.nexps);
      if ~success1,
        obj.RemoveExpDirs(obj.nexps);
        return;
      end
      [success1,msg] = obj.LoadGTLabelsFromFile(obj.nexps);
      if ~success1,
        obj.RemoveExpDirs(obj.nexps);
        return;
      end

      [success1,msg1,missingfiles] = obj.UpdateStatusTable('',obj.nexps);
      missingfiles = missingfiles{obj.nexps};
      if ~success1,
        msg = msg1;
        obj.RemoveExpDirs(obj.nexps);
        return;
      end
      
      % check for existence of necessary files in this directory
      if ~obj.filesfixable,
        msg = sprintf(['Experiment %s is missing required files that cannot '...
          'be generated within this interface. Removing...'],expdir);
        success = false;
        % undo
        obj.RemoveExpDirs(obj.nexps);
        return;
      end
      
      if obj.filesfixable && ~obj.allfilesexist,
        if ~isdeployed 
          if isempty(obj.GetGenerateMissingFiles) || ~obj.GetGenerateMissingFiles()
            res = questdlg(sprintf(['Experiment %s is missing required files:%s. '...
              'Generate now?'],expdir,sprintf(' %s',missingfiles{:})),...
              'Generate missing files?','Yes','Cancel','Yes');
            if strcmpi(res,'Yes')
              obj.SetGenerateMissingFiles();
            end
          else obj.GetGenerateMissingFiles()
            res = 'Yes';
          end
        else
          res = 'Yes';
        end
        
        if strcmpi(res,'Yes'),
          [success,msg] = obj.GenerateMissingFiles(obj.nexps);
          if ~success,
            msg = sprintf(['Error generating missing required files %s '...
              'for experiment %s: %s. Removing...'],...
              sprintf(' %s',missingfiles{:}),expdir,msg);
            obj.RemoveExpDirs(obj.nexps);
            return;
          end
          
        else
          obj.RemoveExpDirs(obj.nexps);
          return;
        end
      end
      
      % Convert the scores file into perframe files.
      
      for i = 1:numel(obj.scoresasinput)
        [success,msg] = obj.ScoresToPerframe(obj.nexps,obj.scoresasinput(i).scorefilename,...
          obj.scoresasinput(i).ts);
          if ~success,
            obj.RemoveExpDirs(obj.nexps);
            return;
          end
      end
      
%       for i = 1:numel(obj.allperframefns),
%         fn = obj.allperframefns{i};
%         if numel(fn)>7 && strcmpi('score',fn(1:5))
%           [success,msg] = obj.ScoresToPerframe(obj.nexps,fn);
%           if ~success,
%             obj.RemoveExpDirs(obj.nexps);
%             return;
%           end
%         end
%       end
%       
      % preload this experiment if this is the first experiment added
      if obj.nexps == 1,
        % TODO: make this work with multiple flies
        [success1,msg1] = obj.PreLoad(1,1);
        if ~success1,
          msg = sprintf('Error getting basic trx info: %s',msg1);
          uiwait(warndlg(msg));
          obj.RemoveExpDirs(obj.nexps);
          obj.ClearStatus();
          return;
        end
      elseif istrxinfo,
        obj.nflies_per_exp(end+1) = nflies_per_exp;
        obj.sex_per_exp{end+1} = sex_per_exp;
        obj.frac_sex_per_exp{end+1} = frac_sex_per_exp;
        obj.firstframes_per_exp{end+1} = firstframes_per_exp;
        obj.endframes_per_exp{end+1} = endframes_per_exp;
        
%         if obj.nexps == 1 % This will set hassex and hasperframesex.
%           [success1,msg1] = obj.GetTrxInfo(obj.nexps,true,obj.trx);
%           if ~success1,
%             msg = sprintf('Error getting basic trx info: %s',msg1);
%             obj.RemoveExpDirs(obj.nexps);
%             return;
%           end
%         end
        
      else
        obj.nflies_per_exp(end+1) = nan;
        obj.sex_per_exp{end+1} = {};
        obj.frac_sex_per_exp{end+1} = struct('M',{},'F',{});
        obj.firstframes_per_exp{end+1} = [];
        obj.endframes_per_exp{end+1} = [];
        [success1,msg1] = obj.GetTrxInfo(obj.nexps);
        if ~success1,
          msg = sprintf('Error getting basic trx info: %s',msg1);
          obj.RemoveExpDirs(obj.nexps);
          return;
        end
      end
      
      obj.InitPredictiondata(obj.nexps);
      
      [success1,msg1] = obj.PreLoadLabeledData();
      if ~success1,
        msg = msg1;
        obj.RemoveExpDirs(obj.nexps);
        return;
      end
      
      
      % save default path
      obj.defaultpath = expdir;
      
      success = true;
      
    end
   
    function [success,msg] = AddExpDirNoPreload(obj,expdir,outexpdir,nflies_per_exp,...
        sex_per_exp,frac_sex_per_exp,firstframes_per_exp,endframes_per_exp)
    % [success,msg] = AddExpDir(obj,expdir,outexpdir,nflies_per_exp,firstframes_per_exp,endframes_per_exp)
    % Add a new experiment to the GUI. If this is the first experiment,
    % then it will be preloaded. 

      success = false; msg = '';
      
      if isnumeric(expdir), return; end
      
      if nargin < 2,
        error('Usage: obj.AddExpDirs(expdir,[outexpdir],[nflies_per_exp])');
      end

      % make sure directory exists
      if ~exist(expdir,'file'),
        msg = sprintf('expdir %s does not exist',expdir);
        return;
      end
      
      isoutexpdir = nargin > 2 && ~isnumeric(outexpdir);
      istrxinfo = nargin > 7 && ~isempty(nflies_per_exp);

      % base name
      [~,expname] = myfileparts(expdir);
      
      % expnames and rootoutputdir must match
      if isoutexpdir,
        [rootoutputdir,outname] = myfileparts(outexpdir); %#ok<*PROP>
        if ~strcmp(expname,outname),
          msg = sprintf('expdir and outexpdir do not match base names: %s ~= %s',expname,outname);
          return;
        end
%         if ischar(obj.rootoutputdir) && ~strcmp(rootoutputdir,obj.rootoutputdir),
%           msg = sprintf('Inconsistent root output directory: %s ~= %s',rootoutputdir,obj.rootoutputdir);
%           return;
%         end
      elseif ~ischar(obj.rootoutputdir),
        outexpdir = expdir;
        rootoutputdir = 0;
      else
        rootoutputdir = obj.rootoutputdir;        
      end
      
      if ischar(obj.rootoutputdir) && ~isoutexpdir,
        outexpdir = fullfile(rootoutputdir,expname);
      end
      
      % create missing outexpdirs
      if ~exist(outexpdir,'dir'),
        [success1,msg1] = mkdir(rootoutputdir,expname);
        if ~success1,
          msg = (sprintf(['Could not create output directory %s, failed to '...
            'set expdirs: %s'],outexpdir,msg1));
          return;
        end
      end

      % create clips dir
      clipsdir = obj.GetFileName('clipsdir');
      outclipsdir = fullfile(outexpdir,clipsdir);

      % okay, checks succeeded, start storing stuff
      obj.nexps = obj.nexps + 1;
      obj.expdirs{end+1} = expdir;
      obj.expnames{end+1} = expname;
      %obj.rootoutputdir = rootoutputdir;
      obj.outexpdirs{end+1} = outexpdir;
      
      % load labels for this experiment
      [success1,msg] = obj.LoadLabelsFromFile(obj.nexps);
      if ~success1,
        obj.RemoveExpDirs(obj.nexps);
        return;
      end
      [success1,msg] = obj.LoadGTLabelsFromFile(obj.nexps);
      if ~success1,
        obj.RemoveExpDirs(obj.nexps);
        return;
      end

      [success1,msg1] = obj.UpdateStatusTable('',obj.nexps);
      if ~success1,
        msg = msg1;
        obj.RemoveExpDirs(obj.nexps);
        return;
      end
      
      % check for existence of necessary files in this directory
      if ~obj.filesfixable,
        msg = sprintf(['Experiment %s is missing required files that cannot '...
          'be generated within this interface. Removing...'],expdir);
        success = false;
        % undo
        obj.RemoveExpDirs(obj.nexps);
        return;
      end
      
      if obj.filesfixable && ~obj.allfilesexist,
%        if ~isdeployed
%          res = questdlg(sprintf('Experiment %s is missing required files. ...
%           Generate now?',expdir),'Generate missing files?','Yes','Cancel','Yes');
%        else
          res = 'Yes';
%        end
        if strcmpi(res,'Yes'),
          [success,msg] = obj.GenerateMissingFiles(obj.nexps,false);
          if ~success,
            msg = sprintf('Error generating missing required files for experiment %s: %s. Removing...',expdir,msg);
            obj.RemoveExpDirs(obj.nexps);
            return;
          end
          
        else
          obj.RemoveExpDirs(obj.nexps);
        end
      end
      
      for i = 1:numel(obj.scoresasinput)
        [success,msg] = obj.ScoresToPerframe(obj.nexps,obj.scoresasinput(i).scorefilename,...
          obj.scoresasinput(i).ts);
        if ~success,
          obj.RemoveExpDirs(obj.nexps);
          return;
        end
      end
      % Convert the scores file into perframe files.
%       for i = 1:numel(obj.allperframefns),
%         fn = obj.allperframefns{i};
%         if numel(fn)>7 && strcmpi('score',fn(1:5))
%           obj.ScoresToPerframe(obj.nexps,fn);
%         end
%       end
%       

      
      % preload this experiment if this is the first experiment added
      if obj.nexps == 1,
        % TODO: make this work with multiple flies
        [success1,msg1] = obj.PreLoad(1,1);
        if ~success1,
          msg = sprintf('Error getting basic trx info: %s',msg1);
          obj.RemoveExpDirs(obj.nexps);
          return;
        end
      elseif istrxinfo,
        obj.nflies_per_exp(end+1) = nflies_per_exp;
        obj.sex_per_exp{end+1} = sex_per_exp;
        obj.frac_sex_per_exp{end+1} = frac_sex_per_exp;
        obj.firstframes_per_exp{end+1} = firstframes_per_exp;
        obj.endframes_per_exp{end+1} = endframes_per_exp;
        
%         if obj.nexps == 1 % This will set hassex and hasperframesex.
%           [success1,msg1] = obj.GetTrxInfo(obj.nexps,true,obj.trx);
%           if ~success1,
%             msg = sprintf('Error getting basic trx info: %s',msg1);
%             obj.RemoveExpDirs(obj.nexps);
%             return;
%           end
%         end
        
      else
        obj.nflies_per_exp(end+1) = nan;
        obj.sex_per_exp{end+1} = {};
        obj.frac_sex_per_exp{end+1} = struct('M',{},'F',{});
        obj.firstframes_per_exp{end+1} = [];
        obj.endframes_per_exp{end+1} = [];
        [success1,msg1] = obj.GetTrxInfo(obj.nexps);
        if ~success1,
          msg = sprintf('Error getting basic trx info: %s',msg1);
          obj.RemoveExpDirs(obj.nexps);
          return;
        end
      end      
      
      obj.InitPredictiondata(obj.nexps);
   
    % save default path
      obj.defaultpath = expdir;
      
      success = true;
      
    end
    
    function [success,msg] = RemoveExpDirs(obj,expi)
      % [success,msg] = RemoveExpDirs(obj,expi)
    % Removes experiments in expi from the GUI. If the currently loaded
    % experiment is removed, then a different experiment may be preloaded. 

      success = false;
      msg = '';
      
      if any(obj.nexps < expi) || any(expi < 1),
        msg = sprintf('expi = %s must be in the range 1 < expi < nexps = %d',mat2str(expi),obj.nexps);
        return;
      end
      
      newExpNumbers = [];
      for ndx = 1:obj.nexps
        if ismember(ndx,expi);
          newExpNumbers(1,ndx) = 0;
        else
          newExpNumbers(1,ndx) = ndx-nnz(expi<ndx);
        end
      end
      
      if ~(numel(obj.expdirs)<expi); obj.expdirs(expi) = []; end
      if ~(numel(obj.expnames)<expi); obj.expnames(expi) = []; end
      if ~(numel(obj.outexpdirs)<expi); obj.outexpdirs(expi) = []; end
      if ~(numel(obj.nflies_per_exp)<expi); obj.nflies_per_exp(expi) = []; end
      if ~(numel(obj.sex_per_exp)<expi); obj.sex_per_exp(expi) = []; end
      if ~(numel(obj.frac_sex_per_exp)<expi); obj.frac_sex_per_exp(expi) = []; end
      if ~(numel(obj.firstframes_per_exp)<expi); obj.firstframes_per_exp(expi) = []; end
      if ~(numel(obj.endframes_per_exp)<expi); obj.endframes_per_exp(expi) = []; end
      if ~(numel(obj.labels)<expi); obj.labels(expi) = []; end
      if ~(numel(obj.labelstats)<expi); obj.labelstats(expi) = []; end
      if ~(numel(obj.gt_labels)<expi); obj.gt_labels(expi) = []; end
      if ~(numel(obj.gt_labelstats)<expi); obj.gt_labelstats(expi) = []; end
      obj.nexps = obj.nexps - numel(expi);
      % TODO: exp2labeloff

      % Clean window features.
      idxcurr = ismember(obj.windowdata.exp, expi);
      obj.windowdata.X(idxcurr,:) = [];
      obj.windowdata.exp(idxcurr) = [];
      obj.windowdata.exp = newExpNumbers(obj.windowdata.exp);
      obj.windowdata.flies(idxcurr) =[];
      obj.windowdata.t(idxcurr) =[];
      obj.windowdata.labelidx_new(idxcurr) = [];
      obj.windowdata.labelidx_imp(idxcurr) = [];
      obj.windowdata.isvalidprediction(...
        idxcurr(1:numel(obj.windowdata.isvalidprediction))) = [];
      obj.windowdata.labelidx_cur(...
        idxcurr(1:numel(obj.windowdata.labelidx_cur))) = [];
      obj.windowdata.predicted(...
        idxcurr(1:numel(obj.windowdata.predicted))) = [];
      obj.windowdata.predicted_probs(...
        idxcurr(1:numel(obj.windowdata.predicted_probs))) = [];
      obj.windowdata.scores(...
        idxcurr(1:numel(obj.windowdata.scores))) = [];
      obj.windowdata.scores_old(...
        idxcurr(1:numel(obj.windowdata.scores_old))) = [];
      obj.windowdata.scores_validated(...
        idxcurr(1:numel(obj.windowdata.scores_validated)),:) = [];
      obj.windowdata.postprocessed(...
        idxcurr(1:numel(obj.windowdata.postprocessed)),:) = [];
      obj.windowdata.distNdx = [];
      obj.windowdata.binVals=[];
      obj.windowdata.exp = newExpNumbers(obj.windowdata.exp);

      idxcurr = ismember(obj.predictdata.exp, expi);
      fnames = fieldnames(obj.predictdata);
      for fndx = 1:numel(fnames)
        obj.predictdata.(fnames{fndx})(idxcurr) = [];
      end
      obj.predictdata.exp = newExpNumbers(obj.predictdata.exp);
      
      if ~isempty(obj.predictblocks.expi)
        idxcurr = ismember(obj.predictblocks.expi,expi);
        fnames = fieldnames(obj.predictblocks);
        for fndx = 1:numel(fnames)
          obj.predictblocks.(fnames{fndx})(idxcurr) = [];
        end
        obj.predictblocks.expi = newExpNumbers(obj.predictblocks.expi);
      end
      
      if ~isempty(obj.randomGTSuggestions)
        obj.randomGTSuggestions(expi) = [];
      end
      
      if ~isempty(obj.loadedGTSuggestions) && numel(obj.loadedGTSuggestions)>=expi
        obj.loadedGTSuggestions(expi) = [];
      end

      % update current exp, flies
      if ~isempty(obj.expi) && obj.expi > 0 && ismember(obj.expi,expi),
        
        % change to different experiment, by default choose fly 1
        % TODO: allow for more than one fly to be selected at once
        obj.expi = 0;
        obj.flies = nan(size(obj.flies));

        if obj.nexps > 0,
          obj.PreLoad(obj.nexps,1);
        end

      elseif ~isempty(obj.expi) && obj.expi > 0 && ~ismember(obj.expi,expi)
        obj.expi = obj.expi - nnz(ismember(1:obj.expi,expi));
      end
      
      success = true;
      
    end

    
% File Handling 
    

    function res = GetFileName(obj,file)
    % res = GetFileName(obj,file)
    % Get base name of file of the input type file.
      switch file,
        case 'movie',
          res = obj.moviefilename;
        case 'trx',
          res = obj.trxfilename;
        case 'label',
          res = obj.labelfilename;
        case 'gt_label',
            res = obj.gt_labelfilename;            
        case {'perframedir','perframe'},
          res = obj.perframedir;
        case {'clipsdir','clips'},
          if ischar(obj.clipsdir)
            res = strtrim(obj.clipsdir);
          else
            res = 'clips';
          end            
        case 'scores',
          res = obj.scorefilename;
        otherwise
          error('Unknown file type %s',file);
      end
    end
    
    function [filename,timestamp] = GetFile(obj,file,expi,dowrite)
        % [filename,timestamp] = GetFile(obj,file,expi)
    % Get the full path to the file of type file for experiment expi. 
  
      if nargin < 4,
        dowrite = false;
      end
      
      % base name
      fn = obj.GetFileName(file);
      
      % if this is an output file, only look in output experiment directory
      if dowrite && JLabelData.IsOutputFile(file),
        expdirs_try = obj.outexpdirs(expi);
      else
        % otherwise, first look in output directory, then look in input
        % directory
        expdirs_try = {obj.outexpdirs{expi},obj.expdirs{expi}};
      end
      
      % initialize timestamp = -inf if we never set
      timestamp = -inf;

      % loop through directories to look in
      for j = 1:numel(expdirs_try),
        expdir = expdirs_try{j}; 
        
        % just one file to look for
        filename = fullfile(expdir,fn);
        if exist(filename,'file'),
          tmp = dir(filename);
          timestamp = tmp.datenum;
          break;
        end
        % check for lnk files
        if ispc && exist([filename,'.lnk'],'file'),
          isseq = ~isempty(regexp(fn,'\.seq$','once'));
          % for seq file, just keep the soft link, get_readframe_fcn will
          % deal with it
          [actualfilename,didfind] = GetPCShortcutFileActualPath(filename);
          if didfind,
            tmp = dir(actualfilename);
            timestamp = tmp.datenum;
            if ~isseq || ~strcmpi(file,'movie'),
              filename = actualfilename;
            end              
            break;
          end
        end
      end
      
    end
    
    function SetGenerateMissingFiles(obj)
      obj.perframeGenerate = true;
    end
    
    function perframeGenerate = GetGenerateMissingFiles(obj)
      perframeGenerate = obj.perframeGenerate;
    end
    
    function [success,msg] = GenerateMissingFiles(obj,expi,isInteractive)
    % [success,msg] = GenerateMissingFiles(obj,expi)
    % Generate required, missing files for experiments expi. 
    % TODO: implement this!
      
      success = true;
      msg = '';
      
      if nargin< 3
        if isdeployed 
          isInteractive = false;
        else
          isInteractive = true;
        end
      end
      
      for i = 1:numel(obj.filetypes),
        file = obj.filetypes{i};
        if obj.IsRequiredFile(file) && obj.CanGenerateFile(file) && ...
            ~obj.FileExists(file,expi),
          fprintf('Generating %s for %s...\n',file,obj.expnames{expi});
          switch file,
%             case 'window',
%               [success1,msg1] = obj.GenerateWindowFeaturesFiles(expi);
%               success = success && success1;
%               if ~success1,
%                 msg = [msg,'\n',msg1]; %#ok<AGROW>
%               end
            case 'perframedir',
              [success1,msg1] = obj.GeneratePerFrameFiles(expi,isInteractive);
              success = success && success1;
              if ~success1,
                msg = [msg,'\n',msg1]; %#ok<AGROW>
              end
          end
        end
      end
      [success1,msg1] = obj.UpdateStatusTable();
      success = success && success1;
      if isempty(msg),
        msg = msg1;
      else
        msg = sprintf('%s\n%s',msg,msg1);
      end
      
    end
    
    function RemoveArenaPFs(obj)
      settings = ReadXMLParams(obj.featureConfigFile);
      toRemove = [];
      for i = 1:numel(obj.allperframefns)
        curpf = obj.allperframefns{i};
        if any(strcmp(curpf,{obj.scoresasinput(:).scorefilename})), continue; end
        curtypes = settings.perframe.(curpf).type;
        if any(strcmpi(curtypes,'arena')) || any(strcmpi(curtypes,'position'))
          toRemove(end+1) = i;
          if isfield(obj.windowfeaturesparams,curpf)
            obj.windowfeaturesparams = rmfield(obj.windowfeaturesparams,curpf);
            obj.windowfeaturescellparams = rmfield(obj.windowfeaturescellparams,curpf);
          end
          curndx = strcmp(obj.curperframefns,curpf);
          obj.curperframefns(curndx) = [];
        end
        
      end
      obj.allperframefns(toRemove) = [];
    end
    
    function [success,msg] = GeneratePerFrameFiles(obj,expi,isInteractive)
      success = false; %#ok<NASGU>
      msg = '';

      perframedir = obj.GetFile('perframedir',expi);
      
      if ~isInteractive,
        dooverwrite = false;
      elseif ~isempty(obj.perframeOverwrite) 
        if obj.perframeOverwrite
          dooverwrite = true;
        else
          dooverwrite = false;
        end
      elseif exist(perframedir,'dir'),
        res = questdlg('Do you want to overwrite existing files or keep them?',...
          'Regenerate files?','Overwrite','Keep','Keep');
        dooverwrite = strcmpi(res,'Overwrite');
        obj.perframeOverwrite = dooverwrite;
      else
        dooverwrite = true;
      end
      
      expdir = obj.expdirs{expi};
      
      if isInteractive
        hwait = mywaitbar(0,sprintf('Initializing perframe directory for %s',expdir),'interpreter','none');
      else
        fprintf('Initializing perframe directory for %s\n',expdir);
      end
      
      perframetrx = Trx('trxfilestr',obj.GetFileName('trx'),...
        'moviefilestr',obj.GetFileName('movie'),...
        'perframedir',obj.GetFileName('perframedir'),...
        'default_landmark_params',obj.landmark_params,...
        'perframe_params',obj.perframe_params,...
        'rootwritedir',obj.rootoutputdir);
      
      perframetrx.AddExpDir(expdir,'dooverwrite',dooverwrite,'openmovie',false);
      
      
      if isempty(fieldnames(obj.landmark_params)) && ~perframetrx.HasLandmarkParams && obj.arenawarn,
        if isInteractive,
          uiwait(warndlg(['Landmark params were not defined in the configuration file'...
            ' or in the trx file. Not computing arena features and removing them from the perframe list']));
        else
          fprintf('Landmark params were not defined in the configuration file. Not computing arena features and removing them from the perframe list');
        end
        obj.RemoveArenaPFs();
        obj.arenawarn = false;
      end
      
      perframefiles = obj.GetPerframeFiles(expi);
      for i = 1:numel(obj.allperframefns),
        fn = obj.allperframefns{i};
        %ndx = find(strcmp(fn,obj.allperframefns));
        file = perframefiles{i};
        if ~dooverwrite && exist(file,'file'),
          continue;
        end
        if isInteractive
          hwait = mywaitbar(i/numel(obj.allperframefns),hwait,...
            sprintf('Computing %s and saving to file %s',fn,file));
        else
          fprintf('Computing %s and saving to file %s\n',fn,file);
        end
        
        % Don't generate the per-frame files from scores here anymore..
        if ~any(strcmp(fn,{obj.scoresasinput(:).scorefilename}))
          perframetrx.(fn);
        end        
      end
      
      if isInteractive && ishandle(hwait),
        delete(hwait);
      end
      
      success = true;
      
    end
    
    
    
    function [success, msg] = ScoresToPerframe(obj,expi,fn,ts)
      success = true; msg = '';
      outdir = obj.outexpdirs{expi};
      scoresFileIn = [fullfile(outdir,fn) '.mat'];
      scoresFileOut = [fullfile(outdir,obj.GetFileName('perframe'),fn) '.mat'];
      if ~exist(scoresFileIn,'file'),
        success = false; 
        msg = sprintf('Scores file %s does not exist to be used as perframe feature',scoresFileIn);
      end
      Q = load(scoresFileIn);
      if Q.timestamp ~= ts, % check the timestamps match the classifier's timestamp.
        success = false; 
        msg = sprintf(['The scores file %s were generated using a classifier' ...
          'that was saved on %s while the classifier chosen was saved on %s'],...
          scoresFileIn,datestr(Q.timestamp),datestr(ts));
      end
      OUT = struct();
      OUT.units = struct(); OUT.units.num = {'scores'};
      OUT.units.den = {''};
      for ndx = 1:numel(Q.allScores.scores)
        t0 = Q.allScores.tStart(ndx);
        t1 = Q.allScores.tEnd(ndx);
        OUT.data{ndx} = Q.allScores.scores{ndx}(t0:t1);
      end
      try
        save(scoresFileOut,'-struct','OUT');
      catch ME,
        success = false;
        msg = ME.message;
      end
    end
    
    function [success,msg] = SetFeatureConfigFile(obj,configfile)
      success = false;
      msg = '';
      
      obj.featureConfigFile = configfile;
      settings = ReadXMLParams(configfile);
      
      if isfield(settings,'perframe_params'),
        pf_fields = fieldnames(settings.perframe_params);
        for ndx = 1:numel(pf_fields),
          obj.perframe_params.(pf_fields{ndx}) = settings.perframe_params.(pf_fields{ndx});
        end
      end

      obj.allperframefns =  fieldnames(settings.perframe);
      
      if isempty(obj.allperframefns)
        msg = 'No perframefns defined';
        return;
      end
      success = true;
      
    end
    
    function [success,msg] = SetFeatureParamsFileName(obj,featureparamsfilename)
    % [success,msg] = SetFeatureParamsFileName(obj,featureparamsfilename)
    % Sets the name of the file describing the features to use to
    % featureparamsfilename. These parameters are read in. Currently, the
    % window data and classifier, predictions are not changed. (TODO)

    success = false;
      msg = '';
      
      if ischar(obj.featureparamsfilename) && strcmp(featureparamsfilename,obj.featureparamsfilename),
        success = true;
        return;
      end
      
      if obj.nexps > 0,
        msg = 'Currently, feature params file can only be changed when no experiments are loaded';
        return;
      end
      
      if ~exist(featureparamsfilename,'file'),
        success = true; 
        msg = '';
        return;
      end

      
%       try
        [windowfeaturesparams,windowfeaturescellparams,basicFeatureTable,featureWindowSize] = ...
          ReadPerFrameParams(featureparamsfilename,obj.featureConfigFile); %#ok<PROP>
%       catch ME,
%         msg = sprintf('Error reading feature parameters file %s: %s',...
%           params.featureparamsfilename,getReport(ME));
%         return;
%       end
      obj.SetPerframeParams(windowfeaturesparams,windowfeaturescellparams); %#ok<PROP>
      obj.featureparamsfilename = featureparamsfilename;
      obj.basicFeatureTable = basicFeatureTable;
      obj.featureWindowSize = featureWindowSize;
      success = true;
    end
    
    function SetPerframeParams(obj,windowfeaturesparams,windowfeaturescellparams)
      obj.windowfeaturesparams = windowfeaturesparams; %#ok<PROP>
      obj.windowfeaturescellparams = windowfeaturescellparams; %#ok<PROP>
      obj.curperframefns = fieldnames(windowfeaturesparams);
    end  
    
    
    function ret = NeedSaveProject(obj)
      ret = obj.savewindowfeatures;
    end
    
    function ResetSaveProject(obj)
      obj.savewindowfeatures = false;
    end
    
    function SaveProject(obj)
      configfilename = obj.configfilename;
      [~,~,ext] = fileparts(configfilename);
      if strcmp(ext,'.xml'),
        uiwait(warndlg('Project file is saved in the old format. Cannot save the window features to the project file'));
        return;
      end
      windowfeatures = struct('windowfeaturesparams',obj.windowfeaturesparams,...
        'windowfeaturescellparams',obj.windowfeaturescellparams,...
        'basicFeatureTable',obj.basicFeatureTable,...
        'featureWindowSize',obj.featureWindowSize); %#ok<NASGU>
      if exist(configfilename,'file')
        [didbak,msg] = copyfile(configfilename,[configfilename '~']);
        if ~didbak,
          warning('Could not create backup of %s: %s',configfilename,msg);
        end
      end
      
      save(configfilename,'windowfeatures','-append');
      obj.ResetSaveProject(obj);
    end
    
    function [windowfeaturesparams,windowfeaturescellparams] = GetPerframeParams(obj)
      windowfeaturesparams = obj.windowfeaturesparams; %#ok<PROP>
      windowfeaturescellparams = obj.windowfeaturescellparams; %#ok<PROP>
    end  
    
    function [filenames,timestamps] = GetPerframeFiles(obj,expi,dowrite)
    % [filenames,timestamps] = GetPerFrameFiles(obj,file,expi)
    % Get the full path to the per-frame mat files for experiment expi
      
      if nargin < 3,
        dowrite = false;
      end
      
      fn = obj.GetFileName('perframedir');
      
      % if this is an output file, only look in output experiment directory
      if dowrite && JLabelData.IsOutputFile('perframedir'),
        expdirs_try = obj.outexpdirs(expi);
      else
        % otherwise, first look in output directory, then look in input
        % directory
        expdirs_try = {obj.outexpdirs{expi},obj.expdirs{expi}};
      end
      
      filenames = cell(1,numel(obj.allperframefns));
      timestamps = -inf(1,numel(obj.allperframefns));
      
      for i = 1:numel(obj.allperframefns),

        % loop through directories to look in
        for j = 1:numel(expdirs_try),
          expdir = expdirs_try{j};
          perframedir = fullfile(expdir,fn);
          if ispc && ~exist(perframedir,'dir'),
            [actualperframedir,didfind] = GetPCShortcutFileActualPath(perframedir);
            if didfind,
              perframedir = actualperframedir;
            end
          end
          filename = fullfile(perframedir,[obj.allperframefns{i},'.mat']);
          if ispc && ~exist(filename,'file'),
            [actualfilename,didfind] = GetPCShortcutFileActualPath(filename);
            if didfind,
              filename = actualfilename;
            end
          end
          
          if exist(filename,'file'),
            filenames{i} = filename;
            tmp = dir(filename);
            timestamps(i) = tmp.datenum;
          elseif j == 1,
            filenames{i} = filename;
          end
          
        end
      end
    end
    
    function [success,msg,missingfiles] = UpdateStatusTable(obj,filetypes,expis)
    % [success,msg] = UpdateStatusTable(obj,filetypes,expis)
    % Update the tables of what files exist for what experiments. This
    % returns false if all files were in existence or could be generated
    % and now they are/can not. 

      msg = '';
      success = false;

      if nargin > 1 && ~isempty(filetypes),
        [ism,fileis] = ismember(filetypes,obj.filetypes);
        if any(~ism),
          msg = 'Unknown filetypes';
          return;
        end
      else
        fileis = 1:numel(obj.filetypes);
      end
      if nargin <= 2 || isempty(expis),
        expis = 1:obj.nexps;
      end

      missingfiles = cell(1,obj.nexps);
      for i = 1:obj.nexps,
        missingfiles{i} = {};
      end
      
      % initialize fileexists table
      obj.fileexists(expis,fileis) = false;
      obj.filetimestamps(expis,fileis) = nan;
      
      
      % loop through all file types
      for filei = fileis,
        file = obj.filetypes{filei};
        % loop through experiments
        for expi = expis,
          
          if strcmpi(file,'perframedir'),
            [fn,timestamps] = obj.GetPerframeFiles(expi);
            if isempty(fn),
              obj.fileexists(expi,filei) = false;
              obj.filetimestamps(expi,filei) = -inf;
            else
              pfexists = cellfun(@(s) exist(s,'file'),fn);
              obj.fileexists(expi,filei) = all(pfexists);
              if ~obj.fileexists(expi,filei) && JLabelData.IsRequiredFile(file),
                for tmpi = find(~pfexists(:)'),
                  [~,missingfiles{expi}{end+1}] = myfileparts(fn{tmpi});
                  missingfiles{expi}{end} = ['perframe_',missingfiles{expi}{end}];
                end
              end
              obj.filetimestamps(expi,filei) = max(timestamps);
            end
          else
          
            % check for existence of current file(s)
            [fn,obj.filetimestamps(expi,filei)] = obj.GetFile(file,expi);
            if iscell(fn),
              obj.fileexists(expi,filei) = ~isinf(obj.filetimestamps(expi,filei)) || ...
                all(cellfun(@(s) exist(s,'file'),fn));
            else
              obj.fileexists(expi,filei) = ~isinf(obj.filetimestamps(expi,filei)) || exist(fn,'file');
            end
            if ~obj.fileexists(expi,filei) && JLabelData.IsRequiredFile(file)
              missingfiles{expi}{end+1} = file;
            end
              
          end
          
        end
      end

      % store old values to see if latest change broke something
      old_filesfixable = obj.filesfixable;
      old_allfilesexist = obj.allfilesexist;

      % initialize summaries to true
      obj.filesfixable = true;
      obj.allfilesexist = true;

      for filei = 1:numel(obj.filetypes),
        file = obj.filetypes{filei};
        % loop through experiments
        for expi = 1:obj.nexps,
          
          % if file doesn't exist and is required, then not all files exist
          if ~obj.fileexists(expi,filei),
            if JLabelData.IsRequiredFile(file),
              obj.allfilesexist = false;
              % if furthermore file can't be generated, then not fixable
              if ~JLabelData.CanGenerateFile(file),
                obj.filesfixable = false;                
                msg1 = sprintf('%s missing and cannot be generated.',file);
                if isempty(msg),
                  msg = msg1;
                else
                  msg = sprintf('%s\n%s',msg,msg1);
                end
              end
            end
            if ~isempty(missingfiles{expi}),
              msg = [msg,'\n','Missing',sprintf(' %s',missingfiles{expi}{:})];
            end
          end
        end
      end

      % fail if was ok and now not ok
      success = ~(old_allfilesexist || old_filesfixable) || ...
        (obj.allfilesexist || obj.filesfixable);
      
    end

    function [fe,ft] = FileExists(obj,file,expi)
    % [fe,ft] = FileExists(obj,file,expi)
    % Returns whether the input file exists for the input experiment. 
      filei = find(strcmpi(file,obj.filetypes),1);
      if isempty(filei),
        error('file type %s does not match any known file type',file);
      end
      if nargin < 3,
        expi = 1:obj.nexps;
      end
      fe = obj.fileexists(expi,filei);
      ft = obj.filetimestamps(expi,filei);
    end
    
    
% Tracking information   


    function [success,msg] = GetTrxInfo(obj,expi,canusecache,trx)
    % [success,msg] = GetTrxInfo(obj,expi)
    % Fills in nflies_per_exp, firstframes_per_exp, and endframes_per_exp
    % for experiment expi. This may require loading in trajectories. 
      success = true;
      msg = '';
      if nargin < 3,
        canusecache = true;
      end
%       canusecache = false;
      istrxinput = nargin >= 4;
      
      obj.SetStatus('Reading trx info for experiment %s',obj.expdirs{expi});
      if numel(obj.nflies_per_exp) < expi || ...
          numel(obj.sex_per_exp) < expi || ...
          numel(obj.frac_sex_per_exp) < expi || ...
          numel(obj.firstframes_per_exp) < expi || ...
          numel(obj.endframes_per_exp) < expi || ...
          isnan(obj.nflies_per_exp(expi)),
        if ~istrxinput,

          trxfile = fullfile(obj.expdirs{expi},obj.GetFileName('trx'));
          if ~exist(trxfile,'file'),
            msg = sprintf('Trx file %s does not exist, cannot count flies',trxfile);
            success = false;
            return;
          else
          
            if isempty(obj.expi) || obj.expi == 0,
              % TODO: make this work for multiple flies
              obj.PreLoad(expi,1);
              trx = obj.trx;
            elseif canusecache && expi == obj.expi,
              trx = obj.trx;
            else
%               try
                % REMOVE THIS
                global CACHED_TRX; %#ok<TLEV>
                global CACHED_TRX_EXPNAME; %#ok<TLEV>
                if isempty(CACHED_TRX) || isempty(CACHED_TRX_EXPNAME) || ...
                    ~strcmp(obj.expnames{expi},CACHED_TRX_EXPNAME),
                  hwait = mywaitbar(0,sprintf('Loading trx to determine number of flies for %s',...
                    obj.expnames{expi}),'interpreter','none');
                  trx = load_tracks(trxfile);
                  if ishandle(hwait), delete(hwait); end
                  CACHED_TRX = trx;
                  CACHED_TRX_EXPNAME = obj.expnames{expi};
                else
                  fprintf('DEBUG: Using CACHED_TRX. REMOVE THIS\n');
                  trx = CACHED_TRX;
                end
%               catch ME,
%                 msg = sprintf(['Could not load trx file for experiment %s '...
%                     'to count flies: %s'],obj.expdirs{expi},getReport(ME));
%               end
            end
          end
        end
        obj.nflies_per_exp(expi) = numel(trx);
        obj.firstframes_per_exp{expi} = [trx.firstframe];
        obj.endframes_per_exp{expi} = [trx.endframe];

        obj.hassex = obj.hassex || isfield(trx,'sex');
        
        % store sex info
        tmp = repmat({nan},[1,numel(trx)]);
        obj.frac_sex_per_exp{expi} = struct('M',tmp,'F',tmp);
        obj.sex_per_exp{expi} = repmat({'?'},[1,numel(trx)]);
        if isfield(trx,'sex'),
          obj.hasperframesex = iscell(trx(1).sex);
          if obj.hasperframesex,
            for fly = 1:numel(trx),
              n = numel(trx(fly).sex);
              nmale = nnz(strcmpi(trx(fly).sex,'M'));
              nfemale = nnz(strcmpi(trx(fly).sex,'F'));
              obj.frac_sex_per_exp{expi}(fly).M = nmale/n;
              obj.frac_sex_per_exp{expi}(fly).F = nfemale/n;
              if nmale > nfemale,
                obj.sex_per_exp{expi}{fly} = 'M';
              elseif nfemale > nmale,
                obj.sex_per_exp{expi}{fly} = 'F';
              else
                obj.sex_per_exp{expi}{fly} = '?';
              end
            end
          else
            for fly = 1:numel(trx),
              obj.sex_per_exp{expi}{fly} = trx(fly).sex;
              if strcmpi(trx(fly).sex,'M'),
                obj.frac_sex_per_exp{expi}(fly).M = 1;
                obj.frac_sex_per_exp{expi}(fly).F = 0;
              elseif strcmpi(trx(fly).sex,'F'),
                obj.frac_sex_per_exp{expi}(fly).M = 0;
                obj.frac_sex_per_exp{expi}(fly).F = 1;
              end
            end
          end
        end
      end
      if isfield(trx,'arena')
        obj.hasarenaparams(expi) = true;
      else
        obj.hasarenaparams(expi) = false;
      end
      
      obj.ClearStatus();
      
    end

    function out = GetTrxValues(obj,infoType,expi,flies,ts)
    % A generic function that return track info.

      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,1);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 4,     % No flies given
        switch infoType
          case 'Trx'
            out = obj.trx;
          case 'X'
            out = {obj.trx.x};
          case 'Y'
            out = {obj.trx.y};
          case 'A'
            out = {obj.trx.a};
          case 'B'
            out = {obj.trx.b};
          case 'Theta'
            out = {obj.trx.theta};
          otherwise
            error('Incorrect infotype requested from GetTrxValues with less than 4 arguments');
        end
        return;

      
      elseif nargin < 5, % No ts given
        switch infoType
          case 'Trx'
            out = obj.trx(flies);
          case 'X'
            out = {obj.trx(flies).x};
          case 'Y'
            out = {obj.trx(flies).y};
          case 'A'
            out = {obj.trx(flies).a};
          case 'B'
            out = {obj.trx(flies).b};
          case 'Theta'
            out = {obj.trx(flies).theta};
          case 'X1'
            out = [obj.trx(flies).x];
          case 'Y1'
            out = [obj.trx(flies).y];
          case 'A1'
            out = [obj.trx(flies).a];
          case 'B1'
            out = [obj.trx(flies).b];
          case 'Theta1'
            out = [obj.trx(flies).theta];
          otherwise
            error('Incorrect infotype requested from GetTrxValues');
        end
        return
      else               % Everything is given
        nflies = numel(flies);
        fly = flies(1);
        switch infoType
          case 'Trx'
            c = cell(1,nflies);
            trx = struct('x',c,'y',c,'a',c,'b',c,'theta',c,'ts',c,'firstframe',c,'endframe',c);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              trx(i).x = obj.trx(fly).x(js);
              trx(i).y = obj.trx(fly).y(js);
              trx(i).a = obj.trx(fly).a(js);
              trx(i).b = obj.trx(fly).b(js);
              trx(i).theta = obj.trx(fly).theta(js);
              trx(i).ts = js-obj.trx(fly).off;
              trx(i).firstframe = trx(i).ts(1);
              trx(i).endframe = trx(i).ts(end);
            end
            out = trx;
          case 'X'
            x = cell(1,nflies);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              x{i} = obj.trx(fly).x(js);
            end
            out = x;
          case 'Y'
            x = cell(1,nflies);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              x{i} = obj.trx(fly).y(js);
            end
            out = x;
          case 'A'
            x = cell(1,nflies);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              x{i} = obj.trx(fly).a(js);
            end
            out = x;
          case 'B'
            x = cell(1,nflies);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              x{i} = obj.trx(fly).b(js);
            end
            out = x;
          case 'Theta'
            x = cell(1,nflies);
            for i = 1:numel(flies),
              fly = flies(i);
              js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
              x{i} = obj.trx(fly).theta(js);
            end
            out = x;
          case 'X1'
            out = obj.trx(fly).x(ts + obj.trx(fly).off);
          case 'Y1'
            out = obj.trx(fly).y(ts + obj.trx(fly).off);
          case 'A1'
            out = obj.trx(fly).a(ts + obj.trx(fly).off);
          case 'B1'
            out = obj.trx(fly).b(ts + obj.trx(fly).off);
          case 'Theta1'
            out = obj.trx(fly).theta(ts + obj.trx(fly).off);
          otherwise
            error('Incorrect infotype requested from GetTrxValues');
         end
      end
      
    end
    
    function pos = GetTrxPos1(varargin)
    % [x,y,theta,a,b] = GetTrxPos1(obj,expi,fly,ts)
    % Returns the position for the input experiment, SINGLE fly, and
    % frames. If ts is not input, then all frames are returned. 

      % moved to separate file so that function could be easily modified
      pos = JLabelData_GetTrxPos(varargin{:});

    end

    function sex = GetSex(obj,expi,fly,ts,fast)
    % x = GetSex(obj,expi,fly,ts)
    % Returns the sex for the input experiment, SINGLE fly, and
    % frames. If ts is not input, then all frames are returned. 

      if ~obj.hassex,
        sex = '?';
        return;
      end
      
      if nargin < 5,
        fast = false;
      end
      
      if ~obj.hasperframesex || fast,
        sex = obj.sex_per_exp{expi}(fly);
        return;
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,fly);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end
      
      if nargin < 4,
        sex = obj.trx(fly).sex;
        return;
      end
      
      sex = obj.trx(fly).sex(ts + obj.trx(fly).off);

    end

    function sex = GetSex1(obj,expi,fly,t)
    % x = GetSex1(obj,expi,fly,t)
    % Returns the sex for the input experiment, SINGLE fly, and
    % SINGLE frame. 

      if ~obj.hassex,
        sex = '?';
        return;
      end
            
      if ~obj.hasperframesex,
        sex = obj.sex_per_exp{expi}(fly);
        if iscell(sex),
          sex = sex{1};
        end
        return;
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,fly);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end
            
      sex = obj.trx(fly).sex{t + obj.trx(fly).off};

    end
    
    function sexfrac = GetSexFrac(obj,expi,fly)
    % x = GetSexFrac(obj,expi,fly)
    % Returns a struct indicating the fraction of frames for which the sex
    % of the fly is M, F

      sexfrac = obj.frac_sex_per_exp{expi}(fly);

    end
    
    function t0 = GetTrxFirstFrame(obj,expi,flies)
    % t0 = GetTrxFirstFrame(obj,expi,flies)
    % Returns the firstframes for the input experiment and flies. If flies
    % is not input, then all flies are returned. 
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if nargin < 3,
        t0 = obj.firstframes_per_exp{expi};
        return;
      end

      t0 = obj.firstframes_per_exp{expi}(flies);
      
    end

    function t1 = GetTrxEndFrame(obj,expi,flies)
    % t1 = GetTrxEndFrame(obj,expi,flies)
    % Returns the endframes for the input experiment and flies. If flies
    % is not input, then all flies are returned. 

      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if nargin < 3,
        t1 = obj.endframes_per_exp{expi};
        return;
      end

      t1 = obj.endframes_per_exp{expi}(flies);
      
    end

    function SetConfidenceThreshold(obj,thresholds,ndx)
      obj.confThresholds(ndx) = thresholds;
    end
    
    function thresholds = GetConfidenceThreshold(obj,ndx)
      thresholds =obj.confThresholds(ndx) ;
    end
    
    
% Labels and predictions    
    

    function [success,msg] = PreLoad(obj,expi,flies)
    % [success,msg] = PreLoad(obj,expi,flies)
    % Preloads data associated with the input experiment and flies. If
    % neither the experiment nor flies are changing, then we do nothing. If
    % there is currently a preloaded experiment, then we store the labels
    % in labelidx into labels using StoreLabels. We then load from labels
    % into labelidx for the new experiment and flies. We load the per-frame
    % data for this experiment and flies. If this is a different
    % experiment, then we load in the trajectories for this experiment.  
      
      success = false;
      msg = '';
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end

      if numel(unique(flies)) ~= numel(flies),
        msg = 'flies must all be unique';
        return;
      end
      
      diffexpi = isempty(obj.expi) || expi ~= obj.expi;
      diffflies = diffexpi || numel(flies) ~= numel(obj.flies) || ~all(flies == obj.flies);
      % nothing to do
      if ~diffflies,
        success = true;
        return;
      end

      if ~isempty(obj.expi) && obj.expi > 0,
        % store labels currently in labelidx to labels
        obj.StoreLabels();
      end
      
      if diffexpi,
        
        % load trx
%         try
          trxfilename = obj.GetFile('trx',expi);
          if ~exist(trxfilename,'file')
            msg = sprintf('Trx file %s does not exist',trxfilename);
            success = false;
            return;
          end
          
          obj.SetStatus('Loading trx for experiment %s',obj.expnames{expi});
                    
          % TODO: remove this
          global CACHED_TRX; %#ok<TLEV>
          global CACHED_TRX_EXPNAME; %#ok<TLEV>
          if isempty(CACHED_TRX) || isempty(CACHED_TRX_EXPNAME) || ...
              ~strcmp(obj.expnames{expi},CACHED_TRX_EXPNAME),
            obj.trx = load_tracks(trxfilename);
            CACHED_TRX = obj.trx;
            CACHED_TRX_EXPNAME = obj.expnames{expi};
          else
            fprintf('DEBUG: Using CACHED_TRX. REMOVE THIS\n');
            obj.trx = CACHED_TRX;
          end
          % store trx_info, in case this is the first time these trx have
          % been loaded
          [success,msg] = obj.GetTrxInfo(expi,true,obj.trx);
          if ~success,
            return;
          end
          
%         catch ME,
%           msg = sprintf('Error loading trx from file %s: %s',trxfilename,getReport(ME));
%           if ishandle(hwait),
%             delete(hwait);
%             drawnow;
%           end
%           return;
%         end
 
      end

      % set labelidx from labels
      obj.SetStatus('Caching labels for experiment %s, flies%s',obj.expnames{expi},sprintf(' %d',flies));
      [obj.labelidx,obj.t0_curr,obj.t1_curr] = obj.GetLabelIdx(expi,flies);
      obj.labelidx_off = 1 - obj.t0_curr;
      
      % load perframedata
      obj.SetStatus('Loading per-frame data for %s, flies %s',obj.expdirs{expi},mat2str(flies));
      file = obj.GetPerframeFiles(expi);
      for j = 1:numel(obj.allperframefns),
        if ~exist(file{j},'file'),
          success = false;
          msg = sprintf('Per-frame data file %s does not exist',file{j});
          return;
        end
%         try
          tmp = load(file{j});
          obj.perframedata{j} = tmp.data{flies(1)};
          obj.perframeunits{j} = tmp.units;
%         catch ME,
%           msg = getReport(ME);
%         end
      end
      
      obj.expi = expi;
      obj.flies = flies;

      obj.UpdatePredictedIdx();
      obj.ClearStatus();
           
      success = true;
      
    end
    
    function ClearCachedPerExpData(obj)
    % ClearCachedPerExpData(obj)
    % Clears all cached data for the currently loaded experiment
      obj.trx = {};
      obj.expi = 0;
      obj.flies = nan(size(obj.flies));
      obj.perframedata = {};
      obj.labelidx = struct('vals',[],'imp',[],'timestamp',[]);
      obj.labelidx_off = 0;
      obj.t0_curr = 0;
      obj.t1_curr = 0;
      obj.predictedidx = [];
      obj.scoresidx = [];
      obj.scoresidx_old = [];
      obj.erroridx = [];
      obj.suggestedidx = [];
    end

    function timestamp = GetLabelTimestamps(obj,expis,flies,ts)
      
      timestamp = nan(size(ts));      
      for expi = 1:obj.nexps,
        expidx = expis == expi;
        if ~any(expidx),
          continue;
        end
        % TODO: extend to multiple flies
        for fly = 1:obj.nflies_per_exp(expi),
          flyidx = expidx & flies == fly;
          if ~any(flyidx),
            continue;
          end
          [labelidx,T0] = obj.GetLabelIdx(expi,fly);
          timestamp(flyidx) = labelidx.timestamp(ts(flyidx)-T0+1);
        end
        
        
      end
      
    end
    
    function [labelidx,T0,T1] = GetLabelIdx(obj,expi,flies,T0,T1)
    % [labelidx,T0,T1] = GetLabelIdx(obj,expi,flies)
    % Returns the labelidx for the input experiment and flies read from
    % labels. 

      if ~isempty(obj.expi) && numel(flies) == numel(obj.flies) && obj.IsCurFly(expi,flies),
        if nargin < 4,
          labelidx = obj.labelidx;
          T0 = obj.t0_curr;
          T1 = obj.t1_curr;
        else
          labelidx.vals = obj.labelidx.vals(T0+obj.labelidx_off:T1+obj.labelidx_off);
          labelidx.imp = obj.labelidx.imp(T0+obj.labelidx_off:T1+obj.labelidx_off);
          labelidx.timestamp = obj.labelidx.timestamp(T0+obj.labelidx_off:T1+obj.labelidx_off);
        end
        return;
      end
      
      if nargin < 4,
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      n = T1-T0+1;
      off = 1 - T0;
      labels_curr = obj.GetLabels(expi,flies);
      labelidx.vals = zeros(1,n);
      labelidx.imp = zeros(1,n);
      labelidx.timestamp = zeros(1,n);
      
      for i = 1:obj.nbehaviors,
        for j = find(strcmp(labels_curr.names,obj.labelnames{i})),
          t0 = labels_curr.t0s(j);
          t1 = labels_curr.t1s(j);
          if t0>T1 || t1<T0; continue;end
          t0 = max(T0,t0);
          t1 = min(T1+1,t1);
          labelidx.vals(t0+off:t1-1+off) = i;
          labelidx.timestamp(t0+off:t1-1+off) = labels_curr.timestamp(j); 
        end
      end
      for j = 1:numel(labels_curr.imp_t0s)
        t0 = labels_curr.imp_t0s(j); t1 = labels_curr.imp_t1s(j);
        if t0>T1 || t1<T0; continue;end
        t0 = max(T0,t0);
        t1 = min(T1+1,t1);
        labelidx.imp(t0+off:t1-1+off) = 1;
      end
      
    end

    function [perframedata,T0,T1] = GetPerFrameData(obj,expi,flies,prop,T0,T1)
    % [perframedata,T0,T1] = GetPerFrameData(obj,expi,flies,prop,T0,T1)
    % Returns the per-frame data for the input experiment, flies, and
    % property. 

      if ischar(prop),
        prop = find(strcmp(prop,obj.allperframefns),1);
        if isempty(prop),
          error('Property %s is not a per-frame property');
        end
      end
      
      if obj.IsCurFly(expi,flies) 
        if nargin < 5,
          perframedata = obj.perframedata{prop};
          T0 = obj.t0_curr;
          T1 = obj.t0_curr + numel(perframedata) - 1;
        else
          T0 = max(T0,obj.t0_curr);
          T1 = min(T1,obj.t0_curr+numel(obj.perframedata{prop})-1);
          i0 = T0 - obj.t0_curr + 1;
          i1 = T1 - obj.t0_curr + 1;
          perframedata = obj.perframedata{prop}(i0:i1);
        end
        return;
      end
      
      perframedir = obj.GetFile('perframedir',expi);
      tmp = load(fullfile(perframedir,[obj.allperframefns{prop},'.mat']));
      if nargin < 5,
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        % TODO: generalize to multi-fly
        perframedata = tmp.data{flies(1)};
        T1 = T0 + numel(perframedata) - 1;
        return;
      end
      off = 1 - obj.GetTrxFirstFrame(expi,flies);
      i0 = T0 + off;
      i1 = T1 + off;
      perframedata = tmp.data{flies(1)}(i0:i1);

    end

    function perframedata = GetPerFrameData1(obj,expi,flies,prop,t)
    % perframedata = GetPerFrameData1(obj,expi,flies,prop,t)
    % Returns the per-frame data for the input experiment, flies, and
    % property. 

%       if ischar(prop),
%         prop = find(strcmp(prop,handles.perframefn),1);
%         if isempty(prop),
%           error('Property %s is not a per-frame property');
%         end
%       end
      
      if ~isempty(obj.expi) && expi == obj.expi && numel(flies) == numel(obj.flies) && all(flies == obj.flies),
        is = t-obj.t0_curr+1;
        badidx = is > numel(obj.perframedata{prop});
        if any(badidx),
          perframedata = nan(size(is));
          perframedata(~badidx) = obj.perframedata{prop}(is(~badidx));
        else
          perframedata = obj.perframedata{prop}(is);
        end
        return;
      end
      
      perframedir = obj.GetFile('perframedir',expi);
      tmp = load(fullfile(perframedir,[obj.allperframefns{prop},'.mat']));
      off = 1 - obj.GetTrxFirstFrame(expi,flies);
      perframedata = tmp.data{flies(1)}(t+off);

    end
    
    function [prediction,T0,T1] = GetPredictedIdx(obj,expi,flies,T0,T1)

      if ~isempty(obj.expi) && numel(flies) == numel(obj.flies) && obj.IsCurFly(expi,flies),
        if nargin < 4,
          prediction = struct('predictedidx',obj.predictedidx,...
                              'scoresidx', obj.scoresidx);
          T0 = obj.t0_curr;
          T1 = obj.t1_curr;
        else
          prediction = struct(...
            'predictedidx', obj.predictedidx(T0+obj.labelidx_off:T1+obj.labelidx_off),...
            'scoresidx',  obj.scoresidx(T0+obj.labelidx_off:T1+obj.labelidx_off));
        end
        return;
      end
      
      if nargin < 4,
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      
      n = T1-T0+1;
      off = 1 - T0;
      prediction = struct('predictedidx', zeros(1,n),...
                         'scoresidx', zeros(1,n));
      
      
      if ~isempty(obj.predictdata.exp)
        idxcurr = obj.FlyNdx(expi,flies) & ...
          obj.predictdata.t >= T0 & obj.predictdata.t <= T1 & ...
          obj.predictdata.cur_valid;
        prediction.predictedidx(obj.predictdata.t(idxcurr)+off) = ...
          -sign(obj.predictdata.cur(idxcurr))*0.5+1.5;
        prediction.scoresidx(obj.predictdata.t(idxcurr)+off) = ...
          sign(obj.predictdata.cur(idxcurr));
      end
    end
    
    function scores = GetValidatedScores(obj,expi,flies,T0,T1)
      if nargin<4
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      
      n = T1-T0+1;
      off = 1 - T0;
      scores = zeros(1,n);
      
      if ~isempty(obj.windowdata.scores_validated) 
        idxcurr = obj.FlyNdx(expi,flies) & ...
          obj.windowdata.t >= T0 & obj.windowdata.t <= T1;
        scores(obj.windowdata.t(idxcurr)+off) = ...
          obj.windowdata.scores_validated(idxcurr);
      end
      
    end
        
    function scores = GetLoadedScores(obj,expi,flies,T0,T1)
      if nargin<4
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      
      n = T1-T0+1;
      off = 1 - T0;
      scores = zeros(1,n);

      if ~isempty(obj.predictdata.exp)                 
        idxcurr = obj.FlyNdxPredict(expi,flies) &...
          obj.predictdata.t >= T0 & obj.predictdata.t <= T1;
        scores(obj.predictdata.t(idxcurr)+off) = ...
          obj.predictdata.loaded(idxcurr);      
      end
      
    end
    
    function [scores,predictions] = GetPostprocessedScores(obj,expi,flies,T0,T1)
      if nargin<4
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      
      n = T1-T0+1;
      off = 1 - T0;
      scores = zeros(1,n);
      predictions = zeros(1,n);
      
      if any(obj.predictdata.cur_valid)
        idxcurr = obj.FlyNdxPredict(expi,flies) & obj.predictdata.cur_valid;
        scores(obj.predictdata.t(idxcurr)+off) = obj.predictdata.cur(idxcurr);
        predictions(obj.predictdata.t(idxcurr)+off) = 2 - obj.predictdata.cur_pp(idxcurr);
      else
        idxcurr = obj.FlyNdxPredict(expi,flies) & obj.predictdata.loaded_valid;
        scores(obj.predictdata.t(idxcurr)+off) = obj.predictdata.loaded(idxcurr);      
        predictions(obj.predictdata.t(idxcurr)+off) = 2 - obj.predictdata.loaded_pp(idxcurr);      
      end
      
    end
    
    
    function scores = GetOldScores(obj,expi,flies)
      T0 = max(obj.GetTrxFirstFrame(expi,flies));
      T1 = min(obj.GetTrxEndFrame(expi,flies));
      
      n = T1-T0+1;
      off = 1 - T0;
      scores = zeros(1,n);
      
      if ~isempty(obj.predictdata.exp) 
        idxcurr = obj.FlyNdxPredict(expi,flies) & ...
          obj.predictdata.old_valid;
        scores(obj.predictdata.t(idxcurr)+off) = ...
          obj.predictdata.old(idxcurr);
      end
      
    end
    
    function [idx,T0,T1] = IsBehavior(obj,behaviori,expi,flies,T0,T1)
    % [idx,T0,T1] = IsBehavior(obj,behaviori,expi,flies,T0,T1)
    % Returns whether the behavior is labeled as behaviori for experiment
    % expi, flies from frames T0 to T1. If T0 and T1 are not input, then
    % firstframe to endframe are used. 

      if ~isempty(obj.expi) && expi == obj.expi && numel(flies) == numel(obj.flies) && all(flies == obj.flies),
        if nargin < 4,
          idx = obj.labelidx.vals == behaviori;
          T0 = obj.t0_curr;
          T1 = obj.t1_curr;
        else
          idx = obj.labelidx.vals(T0+obj.labelidx_off:T1+obj.labelidx_off) == behaviori;
        end
        return;
      end
      
      if nargin < 4,
        T0 = max(obj.GetTrxFirstFrame(expi,flies));
        T1 = min(obj.GetTrxEndFrame(expi,flies));
      end
      n = T1-T0+1;
      off = 1 - T0;
      labels_curr = obj.GetLabels(expi,flies);
      idx = false(1,n);
      for j = find(strcmp(labels_curr.names,obj.labelnames{behaviori})),
        t0 = labels_curr.t0s(j);
        t1 = labels_curr.t1s(j);
        idx(t0+off:t1-1+off) = true;
      end
      
    end

    function labels_curr = GetLabels(obj,expi,flies)
    % labels_curr = GetLabels(obj,expi,flies)
    % Returns the labels for the input 

      labels_curr = struct('t0s',[],'t1s',[],'names',{{}},'timestamp',[],'off',0,'imp_t0s',[],'imp_t1s',[]);
      
      if nargin < 2 || isempty(expi),
        expi = obj.expi;
      end
      
      if nargin < 3 || isempty(flies),
        flies = obj.flies;
      end

      % cache these labels if current experiment and flies selected
      if expi == obj.expi && all(flies == obj.flies),
        obj.StoreLabels();
      end
      
      if obj.IsGTMode()
        labelsToUse = obj.gt_labels;
      else
        labelsToUse = obj.labels;
      end

      [ism,fliesi] = ismember(flies,labelsToUse(expi).flies,'rows');
      if ism,
        labels_curr.t0s = labelsToUse(expi).t0s{fliesi};
        labels_curr.t1s = labelsToUse(expi).t1s{fliesi};
        labels_curr.names = labelsToUse(expi).names{fliesi};
        labels_curr.off = labelsToUse(expi).off(fliesi);
        if isfield(labelsToUse(expi),'imp_t0s')
          labels_curr.imp_t0s = labelsToUse(expi).imp_t0s{fliesi};
          labels_curr.imp_t1s = labelsToUse(expi).imp_t1s{fliesi};
        end
        labels_curr.timestamp = labelsToUse(expi).timestamp{fliesi};
      else
%         if expi ~= obj.expi,
%           error('This should never happen -- only should get new labels for current experiment');
%         end
        t0_curr = max(obj.GetTrxFirstFrame(expi,flies));
        labels_curr.off = 1-t0_curr;
      end

      
    end

    function StoreLabels(obj)
    % Store labels cached in labelidx for the current experiment and flies
    % to labels structure. This is when the timestamp on labels gets
    % updated. 
      
      % flies not yet initialized
      if isempty(obj.flies) || all(isnan(obj.flies)) || isempty(obj.labelidx.vals),
        return;
      end
      
      obj.StoreLabels1(obj.expi,obj.flies,obj.labelidx,obj.labelidx_off);
            
      % preload labeled window data while we have the per-frame data loaded
      ts = find(obj.labelidx.vals~=0) - obj.labelidx_off;
      if ~obj.IsGTMode() || isempty(obj.predictdata.exp),
        [success,msg] = obj.PreLoadWindowData(obj.expi,obj.flies,ts);
        if ~success,
          warning(msg);
        end
      end
      % update windowdata's labelidx_new
      if ~isempty(obj.windowdata.exp),
        idxcurr = obj.windowdata.exp == obj.expi & ...
          all(bsxfun(@eq,obj.windowdata.flies,obj.flies),2);
        obj.windowdata.labelidx_new(idxcurr) = obj.labelidx.vals(obj.windowdata.t(idxcurr)+obj.labelidx_off);
        obj.windowdata.labelidx_imp(idxcurr) = obj.labelidx.imp(obj.windowdata.t(idxcurr)+obj.labelidx_off);
      end
      
      %obj.UpdateWindowDataLabeled(obj.expi,obj.flies);
      
    end

    function StoreLabels1(obj,expi,flies,labelidx,labelidx_off)
      
      % update labels
      newlabels = struct('t0s',[],'t1s',[],'names',{{}},'flies',[],'timestamp',[],'imp_t0s',[],'imp_t1s',[]);
      for j = 1:obj.nbehaviors,
        [i0s,i1s] = get_interval_ends(labelidx.vals==j);
        
        if ~isempty(i0s),
          n = numel(i0s);
          newlabels.t0s(end+1:end+n) = i0s - labelidx_off;
          newlabels.t1s(end+1:end+n) = i1s - labelidx_off;
          newlabels.names(end+1:end+n) = repmat(obj.labelnames(j),[1,n]);
          newlabels.timestamp(end+1:end+n) = labelidx.timestamp(i0s);
        end
      end
      
      [i0s,i1s] = get_interval_ends(labelidx.imp);
      if ~isempty(i0s),
        newlabels.imp_t0s = i0s - labelidx_off;
        newlabels.imp_t1s = i1s - labelidx_off;
      end
      
      % Store labels according to the mode
      if obj.IsGTMode(),
        labelsToUse = 'gt_labels';
        labelstatsToUse = 'gt_labelstats';
      else
        labelsToUse = 'labels';
        labelstatsToUse = 'labelstats';
      end
      
      [ism,j] = ismember(flies,obj.(labelsToUse)(expi).flies,'rows');
      if ~ism,
        j = size(obj.(labelsToUse)(expi).flies,1)+1;
      end

      obj.(labelsToUse)(expi).t0s{j} = newlabels.t0s;
      obj.(labelsToUse)(expi).t1s{j} = newlabels.t1s;
      obj.(labelsToUse)(expi).names{j} = newlabels.names;
      obj.(labelsToUse)(expi).flies(j,:) = flies;
      obj.(labelsToUse)(expi).off(j) = labelidx_off;
      obj.(labelsToUse)(expi).timestamp{j} = newlabels.timestamp;
      obj.(labelsToUse)(expi).imp_t0s{j} = newlabels.imp_t0s;
      obj.(labelsToUse)(expi).imp_t1s{j} = newlabels.imp_t1s;

      % store labelstats
      obj.(labelstatsToUse)(expi).nflies_labeled = numel(unique(obj.(labelsToUse)(expi).flies));
      obj.(labelstatsToUse)(expi).nbouts_labeled = numel(newlabels.t1s);
            
    end

    function isstart = IsLabelStart(obj,expi,flies,ts)
      
      if obj.expi == expi && all(flies == obj.flies),
        isstart = obj.labelidx.vals(ts+obj.labelidx_off) ~= 0 & ...
          obj.labelidx.vals(ts+obj.labelidx_off-1) ~= obj.labelidx.vals(ts+obj.labelidx_off);
      else
        
        if obj.IsGTMode(),
          labelsToUse = 'gt_labels';
        else
          labelsToUse = 'labels';
        end

        [ism,fliesi] = ismember(flies,obj.(labelsToUse)(expi).flies,'rows');
        if ism,
          isstart = ismember(ts,obj.(labelsToUse)(expi).t0s{fliesi});
        else
          isstart = false(size(ts));
        end
      end
      
    end

    function ClearLabels(obj,expi,flies)
      
      if obj.nexps == 0,
        return;
      end
      
      if obj.IsGTMode()
        labelsToUse = 'gt_labels';
        labelstatsToUse = 'gt_labelstats';
      else
        labelsToUse = 'labels';
        labelstatsToUse = 'labelstats';
      end
      
      timestamp = now;
      
      % use all experiments by default
      if nargin < 2,
        expi = 1:obj.nexps;
      end
      
      % delete all flies by default
      if nargin < 3,
        for i = expi(:)',
          obj.(labelsToUse)(expi).t0s = {};
          obj.(labelsToUse)(expi).t1s = {};
          obj.(labelsToUse)(expi).names = {};
          obj.(labelsToUse)(expi).flies = [];
          obj.(labelsToUse)(expi).off = [];
          obj.(labelsToUse)(expi).timestamp = {};
          obj.(labelstatsToUse)(expi).nflies_labeled = 0;
          obj.(labelstatsToUse)(expi).nbouts_labeled = 0;
          obj.(labelsToUse)(expi).imp_t0s = {};
          obj.(labelsToUse)(expi).imp_t1s = {};
        end
      else
        if numel(expi) > 1,
          error('If flies input to ClearLabels, expi must be a single experiment');
        end
        % no labels
        if numel(obj.(labelsToUse)) < expi,
          return;
        end
        % which index of labels
        [~,flyis] = ismember(obj.(labelsToUse)(expi).flies,flies,'rows');
        for flyi = flyis(:)',
          % keep track of number of bouts so that we can update stats
          ncurr = numel(obj.(labelsToUse)(expi).t0s{flyi});
          obj.(labelsToUse)(expi).t0s{flyi} = [];
          obj.(labelsToUse)(expi).t1s{flyi} = [];
          obj.(labelsToUse)(expi).names{flyi} = {};
          obj.(labelsToUse)(expi).timestamp{flyi} = [];
          obj.(labelsToUse)(expi).imp_t0s{flyi} = [];
          obj.(labelsToUse)(expi).imp_t1s{flyi} = [];
          % update stats
          obj.(labelstatsToUse)(expi).nflies_labeled = obj.(labelstatsToUse)(expi).nflies_labeled - 1;
          obj.(labelstatsToUse)(expi).nbouts_labeled = obj.(labelstatsToUse)(expi).nbouts_labeled - ncurr;
        end
      end
      
      % clear labelidx if nec
      if ismember(obj.expi,expi) && ((nargin < 3) || ismember(obj.flies,flies,'rows')),
        obj.labelidx.vals(:) = 0;
        obj.labelidx.imp(:) = 0;
        obj.labelidx.timestamp(:) = 0;
      end
      
      % clear windowdata labelidx_new
      for i = expi(:)',
        if nargin < 3,
          idx = obj.windowdata.exp == i;
        else
          idx = obj.windowdata.exp == i & ismember(obj.windowdata.flies,flies,'rows');
        end
        obj.windowdata.labelidx_new(idx) = 0;
        obj.windowdata.labelidx_imp(idx) = 0;
        obj.UpdateErrorIdx();data
      end
      
    end
    
    function SetLabel(obj,expi,flies,ts,behaviori,important)
    % SetLabel(obj,expi,flies,ts,behaviori)
    % Set label for experiment expi, flies, and frames ts to behaviori. If
    % expi, flies match current expi, flies, then we only set labelidx.
    % Otherwise, we set labels. 
      
      if obj.IsCurFly(expi,flies),
        obj.labelidx.vals(ts+obj.labelidx_off) = behaviori;
        obj.labelidx.imp(ts+obj.labelidx_off) = important;
        obj.labelidx.timestamp(ts+obj.labelidx_off) = now;
      else
        [labelidx,T0] = obj.GetLabelIdx(expi,flies);
        labelidx.vals(ts+1-T0) = behaviori;
        labelidx.imp(ts+1-T0) = important;
        labelidx.timestamp(ts+1-T0) = now;
        obj.StoreLabels1(expi,flies,labelidx,1-T0);        
      end
      
    end
    
    
% Window data computation.

    function InitPredictiondata(obj,expi)
      
      % First remove the old scores if any for that experiment.
      idx = obj.predictdata.exp == expi;
      fnames = fieldnames(obj.predictdata);
      set_to_nan = {'cur','cur_pp','old','old_pp',...
        'loaded','loaded_pp','timestamp'};
      for ndx = 1:numel(fnames)
        obj.predictdata.(fnames{ndx})(idx) = [];
      end
      
      for flies = 1:obj.nflies_per_exp(expi)
        firstframe = obj.firstframes_per_exp{expi}(flies);
        endframe = obj.endframes_per_exp{expi}(flies);
        nframes = endframe-firstframe+1;
        for sndx = 1:numel(set_to_nan)
          obj.predictdata.(set_to_nan{sndx})(end+1:end+nframes) = nan;
        end
        obj.predictdata.exp(end+1:end+nframes) = expi;
        obj.predictdata.flies(end+1:end+nframes) = flies;
        obj.predictdata.t(end+1:end+nframes) = firstframe:endframe;
        obj.predictdata.cur_valid(end+1:end+nframes) = false;
        obj.predictdata.old_valid(end+1:end+nframes) = false;
        obj.predictdata.loaded_valid(end+1:end+nframes) = false;
      end

    end

    function [success,msg] = PreLoadWindowData(obj,expi,flies,ts)
    % [success,msg] = PreLoadWindowData(obj,expi,flies,ts)
    % Compute and store the window data for experiment expi, flies flies,
    % and all frames ts. 
    % This function finds all frames that currently do not have window data
    % cached. In a loop, it finds the first frame that is missing window
    % data, and computes window data for all frames in a chunk of size
    % 2*obj.windowdatachunk_radius + 1 after this frame using the function
    % ComputeWindowDataChunk. Then, it updates the frames that are missing
    % window data. It proceeds in this loop until there are not frames
    % in the input ts missing window data. 
      
      success = false; msg = '';
      obj.CheckExp(expi); obj.CheckFlies(flies);
      
      % which frames don't have window data yet
      if isempty(obj.windowdata.exp),
        missingts = ts;
        tscurr = [];
      else      
        idxcurr = obj.FlyNdx(expi,flies);
        tscurr = obj.windowdata.t(idxcurr);
        missingts = setdiff(ts,tscurr);
      end
        
      % no frames missing data?
      if isempty(missingts),
        success = true;
        return;
      end

      % get labels for current flies -- will be used when filling in
      % windowdata
      [labelidxStruct,t0_labelidx] = obj.GetLabelIdx(expi,flies);

      % total number of frames to compute window data for -- used for
      % showing prctage complete. 
      nts0 = numel(missingts);
      
      while true,

        % choose a frame missing window data
        %t = missingts(1);
        t = median(missingts);
        if ~ismember(t,missingts),
          t = missingts(argmin(abs(t-missingts)));
        end
        
        % update the status
        obj.SetStatus('Computing window data for exp %s, fly%s: %d%% done...',...
          obj.expnames{expi},sprintf(' %d',flies),round(100*(nts0-numel(missingts))/nts0));
        
        % compute window data for a chunk starting at t
        [success1,msg,t0,t1,X,feature_names] = obj.ComputeWindowDataChunk(expi,flies,t,'center');
        if ~success1, warning(msg); return; end

        % only store window data that isn't already cached
        tsnew = t0:t1;
        idxnew = (~ismember(tsnew,tscurr)) &ismember(tsnew,missingts);
        m = nnz(idxnew);
        if m==0; return; end

        % Add this to predict blocks.
        obj.predictblocks.expi(end+1) = expi;
        obj.predictblocks.flies(end+1) = flies;
        obj.predictblocks.t0(end+1) = t0;
        obj.predictblocks.t1(end+1) = t1;
        

        % add to windowdata
        obj.windowdata.X(end+1:end+m,:) = X(idxnew,:);
        obj.windowdata.exp(end+1:end+m,1) = expi;
        obj.windowdata.flies(end+1:end+m,:) = repmat(flies,[m,1]);
        obj.windowdata.t(end+1:end+m,1) = tsnew(idxnew);
        obj.windowdata.labelidx_cur(end+1:end+m,1) = 0;
        tempLabelsNew = labelidxStruct.vals(t0-t0_labelidx+1:t1-t0_labelidx+1);
        obj.windowdata.labelidx_new(end+1:end+m,1) = tempLabelsNew(idxnew);
        tempLabelsImp = labelidxStruct.imp(t0-t0_labelidx+1:t1-t0_labelidx+1);        
        obj.windowdata.labelidx_imp(end+1:end+m,1) = tempLabelsImp(idxnew);        
        obj.windowdata.labelidx_old(end+1:end+m,1) = 0;
        obj.windowdata.predicted(end+1:end+m,1) = 0;
        obj.windowdata.scores(end+1:end+m,1) = 0;
        obj.windowdata.scores_old(end+1:end+m,1) = 0;   
        obj.windowdata.scores_validated(end+1:end+m,1) = 0;           
        obj.windowdata.postprocessed(end+1:end+m,1) = 0;           
        obj.windowdata.isvalidprediction(end+1:end+m,1) = false;

        % remove from missingts all ts that were computed in this chunk
        missingts(missingts >= t0 & missingts <= t1) = [];

        % stop if we're done
        if isempty(missingts),
          obj.ClearStatus();
          break;
        end
        
      end
      
      % Clean the window data.
%       obj.CleanWindowData();
      
      % store feature_names -- these shouldn't really change
      obj.windowdata.featurenames = feature_names;
      
      success = true;
%       obj.TrimWindowData();
      
    end

    function [success,msg,t0,t1,X,feature_names] = ComputeWindowDataChunk(obj,expi,flies,t,mode,forceCalc)
    % [success,msg,t0,t1,X,feature_names] = ComputeWindowDataChunk(obj,expi,flies,t)
    % Computes a chunk of windowdata near frame t for experiment expi and
    % flies flies. if mode is 'start', then the chunk will start at t. if
    % it is 'center', the chunk will be centered at t. if mode is 'end',
    % the chunk will end at t. by default, mode is 'center'. 
    % t0 and t1 define the bounds of the chunk of window data computed. X
    % is the nframes x nfeatures window data, feature_names is a cell array
    % of length nfeatures containing the names of each feature. 
    %
    % This function first chooses an interval of frames around t, depending 
    % on the mode. it then chooses a subinterval of this interval that
    % covers all frames in this interval that do not have window data. This
    % defines t0 and t1. 
    % 
    % It then loops through all the per-frame features, and calls
    % ComputeWindowFeatures to compute all the window data for that
    % per-frame feature. 
    %
    % To predict over the whole movie we use forceCalc which
    % forces the function to recalculate all the features even though they
    % were calculated before.
      
      success = false; msg = '';
      
      if ~exist('mode','var'), mode = 'center'; end
      if ~exist('forceCalc','var'), forceCalc = false; end
      
      % Check if the features have been configured.
      if isempty(fieldnames(obj.windowfeaturesparams))
        obj.ShowSelectFeatures();
        if isempty(fieldnames(obj.windowfeaturesparams))
          error('No features selected!');
        end
      end
      
      % choose frames to compute:
      
      % bound at start and end frame of these flies
      T0 = max(obj.GetTrxFirstFrame(expi,flies));
      T1 = min(obj.GetTrxEndFrame(expi,flies));
      
      switch lower(mode),
        case 'center',
          % go forward r to find the end of the chunk
          t1 = min(t+obj.windowdatachunk_radius,T1);
          % go backward 2*r to find the start of the chunk
          t0 = max(t1-2*obj.windowdatachunk_radius,T0);
          % go forward 2*r again to find the end of the chunk
          t1 = min(t0+2*obj.windowdatachunk_radius,T1);
        case 'start',
          t0 = max(t,T0);
          t1 = min(t0+2*obj.windowdatachunk_radius,T1);
        case 'end',
          t1 = min(t,T1);
          t0 = max(t1-2*obj.windowdatachunk_radius,T0);
        otherwise
          error('Unknown mode %s',mode);
      end
      
      % find a continuous interval that covers all uncomputed ts between t0
      % and t1
      off = 1-t0;
      n = t1-t0+1;
      docompute = true(1,n);
      if ~isempty(obj.windowdata.exp) && ~forceCalc,
        tscomputed = obj.windowdata.t(obj.FlyNdx(expi,flies));
        tscomputed = tscomputed(tscomputed >= t0 & tscomputed <= t1);
        docompute(tscomputed+off) = false;
      end
      
      X = [];
      feature_names = {};
      if ~any(docompute),
        t1 = t0-1;
        success = true;
        return;
      end
      
      t0 = find(docompute,1,'first') - off;
      t1 = find(docompute,1,'last') - off;
      i0 = t0 - obj.GetTrxFirstFrame(expi,flies) + 1;
      i1 = t1 - obj.GetTrxFirstFrame(expi,flies) + 1;
      
      %       try
      
      curperframefns = obj.curperframefns;
      allperframefns = obj.allperframefns;
      perframeInMemory = ~isempty(obj.flies) && obj.IsCurFly(expi,flies);
      perframedata_all = obj.perframedata;
      perframefile = obj.GetPerframeFiles(expi);
      x_curr_all = cell(1,numel(curperframefns));
      feature_names_all = cell(1,numel(curperframefns));
      windowfeaturescellparams = obj.windowfeaturescellparams;
      
      % loop through per-frame fields to check that they exist.
      for j = 1:numel(curperframefns),
        fn = curperframefns{j};
        
        % get per-frame data
        ndx = find(strcmp(fn,allperframefns));
        if isempty(ndx),
          success = false;
          msg = 'Window features config file has a perframe feature that is not defined in params file';
          return;
        end
        
        if ~exist(perframefile{ndx},'file'),
          if ~isdeployed 
            if isempty(obj.GetGenerateMissingFiles) || ~obj.GenerateMissingFiles()
              res = questdlg(sprintf(['Experiment %s is missing some perframe files '...
                '(%s, possibly more). Generate now?'],obj.expnames{expi},perframefile{ndx}),...
                'Generate missing files?','Yes','Cancel','Yes');
              if strcmpi(res,'Yes');
                obj.SetGenerateMissingFiles();
              end
            else obj.GetGenerateMissingFiles()
              res = 'Yes';
            end
          else
            res = 'Yes';
          end
          
          if strcmpi(res,'Yes'),
            for ndx = 1:obj.nexps  
              [success1,msg1] = obj.GenerateMissingFiles(ndx);
              if ~success1,
                success = success1; msg = msg1;
                return;
              end
            end
            
          else
            success = false;
            msg = sprintf('Cannot compute window data for %s ',expdir);
            return;
          end

          
        end
        
      end
      
      parfor j = 1:numel(curperframefns),
        fn = curperframefns{j};
        
        % get per-frame data
        ndx = find(strcmp(fn,allperframefns));
        if perframeInMemory,
          perframedata = perframedata_all{ndx};
        else
          perframedata = load(perframefile{ndx});
          perframedata = perframedata.data{flies(1)};
        end
        
        i11 = min(i1,numel(perframedata));
        [x_curr,feature_names_curr] = ...
          ComputeWindowFeatures(perframedata,windowfeaturescellparams.(fn){:},'t0',i0,'t1',i11);
        if any(imag(x_curr(:)))
          fprintf('Feature values are complex, check input\n');
        end
        
        if i11 < i1,
          x_curr(:,end+1:end+i1-i11) = nan;
        end
        
        x_curr_all{j} = x_curr;
        feature_names_all{j} = feature_names_curr;
        
      end
      
      for j = 1:numel(curperframefns),
        fn = curperframefns{j};
        x_curr = x_curr_all{j};
        feature_names_curr = feature_names_all{j};
        % add the window data for this per-frame feature to X
        nold = size(X,1);
        nnew = size(x_curr,2);
        if nold > nnew,
          warning(['Number of examples for per-frame feature %s does not match number '...
            'of examples for previous features'],fn);
          x_curr(:,end+1:end+nold-nnew) = nan;
        elseif nnew > nold && ~isempty(X),
          warning(['Number of examples for per-frame feature %s does not match number '...
            'of examples for previous features'],fn);
          X(end+1:end+nnew-nold,:) = nan;
        end
        X = [X,x_curr']; %#ok<AGROW>
        % add the feature names
        feature_names = [feature_names,cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false)]; %#ok<AGROW>
      end
      %       catch ME,
      %         msg = getReport(ME);
      %         return;
      %       end
      X = single(X);
      success = true;
     
    end
    
    function ClearWindowData(obj)
      % Clears window features and predictions for a clean start when selecting
      % features.
      obj.windowdata.X = [];
      obj.windowdata.exp = [];
      obj.windowdata.flies=[];
      obj.windowdata.t=[];
      obj.windowdata.labelidx_cur=[];
      obj.windowdata.labelidx_new=[];
      obj.windowdata.labelidx_imp=[];
      obj.windowdata.labelidx_old=[];      
      obj.windowdata.featurenames={{}};
      obj.windowdata.predicted=[];
      obj.windowdata.predicted_probs=[];
      obj.windowdata.isvalidprediction=[];
      obj.windowdata.distNdx=[];
      obj.windowdata.scores=[];
      obj.windowdata.scores_old=[];
      obj.windowdata.scores_validated=[];
      obj.windowdata.postprocessed =[];
      obj.windowdata.scoreNorm=[];
      obj.windowdata.binVals=[];
      
      obj.predictblocks.t0 = [];
      obj.predictblocks.t1 = [];
      obj.predictblocks.flies = [];
      obj.predictblocks.expi = [];
      
      obj.UpdatePredictedIdx();

    end
  
    function TrimWindowData(obj,doforce)
      % If the size of windowdata is too large, removes windowdata for
      % unlabeled examples.
%       sizeLimit = obj.cacheSize*1e6;
%       classSize = 4;
%       ratioLimit = 0.2;
%       
%       numUnlabeled = nnz(obj.windowdata.labelidx_new==0);
%       numLabeled = nnz(obj.windowdata.labelidx_new);
%       
%       if (nargin < 2 || ~doforce) && (numel(obj.windowdata.X)*classSize < sizeLimit || numUnlabeled/numLabeled<ratioLimit);
%         return;
%       end
%       
%       idx2remove = obj.windowdata.labelidx_new==0 & ...
%         ~obj.FlyNdx(obj.expi,obj.flies);

      idx2remove = obj.windowdata.labelidx_new==0;
      if ~any(idx2remove); return; end
      obj.windowdata.X(idx2remove,:) = [];
      obj.windowdata.exp(idx2remove,:) = [];
      obj.windowdata.flies(idx2remove,:) = [];
      obj.windowdata.t(idx2remove,:) = [];
      obj.windowdata.labelidx_cur(idx2remove,:) = [];
      obj.windowdata.labelidx_new(idx2remove,:) = [];
      obj.windowdata.labelidx_imp(idx2remove,:) = [];
      obj.windowdata.labelidx_old(idx2remove,:) = [];
      obj.windowdata.binVals = [];
      
    end
    
    function UpdatePerframeParams(obj,params,cellParams,basicFeatureTable,featureWindowSize)
    % Updates the feature params. Called by SelectFeatures
      if ~isempty(obj.classifier),
        hasClassifier = true;
      else
        hasClassifier = false;
      end
      
      obj.SetPerframeParams(params,cellParams);
      if nargin>2
        obj.basicFeatureTable = basicFeatureTable;
        obj.featureWindowSize = featureWindowSize;
      end
      obj.savewindowfeatures = true;

      obj.ClearWindowData();
      obj.classifier = [];
      obj.classifier_old = [];
      obj.PreLoadLabeledData();
      if hasClassifier,
        obj.StoreLabels();
        obj.Train();
      end
      % TODO: remove clearwindow features.
    end
    
    function [success,msg] = PreLoadLabeledData(obj)
    % [success,msg] = PreLoadLabeledData(obj)
    % This function precomputes any missing window data for all labeled
    % training examples by calling PreLoadWindowData on all labeled frames.

      success = false; msg = '';
      
      for expi = 1:obj.nexps,
        if obj.IsGTMode(),
          flies_curr = obj.gt_labels(expi).flies;
        else
          flies_curr = obj.labels(expi).flies;
        end
        for i = 1:size(flies_curr,1),
          flies = flies_curr(i,:);
          labels_curr = obj.GetLabels(expi,flies);
          ts = [];
          
          for j = 1:numel(labels_curr.t0s),
            ts = [ts,labels_curr.t0s(j):labels_curr.t1s(j)-1]; %#ok<AGROW>
          end
          
          [success1,msg] = obj.PreLoadWindowData(expi,flies,ts);
          if ~success1,return;end            
          
        end
      end
      success = true;
      
    end
    
    function ShowSelectFeatures(obj)
      obj.SetStatus('Set the window computation features...');
      selHandle = SelectFeatures(obj,obj.scoresasinput);
      uiwait(selHandle);
      obj.ClearStatus();
    end
    
    function wsize = GetFeatureWindowSize(obj)
      wsize = obj.featureWindowSize;
    end

    function UpdateBoostingBins(obj)
      
      islabeled = obj.windowdata.labelidx_cur ~= 0;
      obj.windowdata.binVals = findThresholds(obj.windowdata.X(islabeled,:),obj.classifier_params);
    end
    
    
% Training and prediction.
    

    function Train(obj,doFastUpdates,timerange)
    % Train(obj)
    % Updates the classifier to reflect the current labels. This involves
    % first loading/precomputing the training features. Then, the clasifier
    % is trained/updated. Finally, predictions for the currently loaded
    % window data are updated. Currently, the only implemented classifier is 
    % random ferns. If the classifier exists, then it is updated instead of
    % retrained from scratch. This involves three steps -- replacing labels
    % for frames which have changed label, removing examples for frames
    % which have been removed the training set, and adding new examples for
    % newly labeled frames. If the classifier has not yet been trained, it
    % is trained from scratch. 
      
      % load all labeled data
      [success,msg] = obj.PreLoadLabeledData();
      if ~success,
        warning(msg);
        return;
      end
      
      islabeled = (obj.windowdata.labelidx_new ~= 0) & (obj.windowdata.labelidx_imp);
      if nargin >= 3 && ~isempty(timerange),
        label_timestamp = obj.GetLabelTimestamps(obj.windowdata.exp(islabeled),...
          obj.windowdata.flies(islabeled,:),obj.windowdata.t(islabeled));
        islabeled(islabeled) = label_timestamp >= timerange(1) & label_timestamp < timerange(2);
      end

      if ~any(islabeled),
        uiwait(warndlg('No frames have been labeled. Not doing any training'));
        return;
      end
      
      switch obj.classifiertype,
      
        case 'ferns',
          if isempty(obj.classifier),
            
            % train classifier
            obj.SetStatus('Training fern classifier from %d examples...',numel(islabeled));

            s = struct2paramscell(obj.classifier_params);
            obj.classifier = fernsClfTrain( obj.windowdata.X(islabeled,:), obj.windowdata.labelidx_new(islabeled), s{:} );
            obj.windowdata.labelidx_cur = obj.windowdata.labelidx_new;
                        
          else
            
            % new data added to windowdata at the end, so classifier.inds still
            % matches windowdata(:,1:Nprev)
            Nprev = numel(obj.windowdata.labelidx_cur);
            Ncurr = numel(obj.windowdata.labelidx_new);
            waslabeled = obj.windowdata.labelidx_cur(1:Nprev) ~= 0;
            islabeled = obj.windowdata.labelidx_new(1:Nprev) ~= 0;
            
            % replace labels for examples that have been relabeled:
            % islabeled & waslabeled will not change
            idx_relabel = islabeled & waslabeled & ...
              (obj.windowdata.labelidx_new(1:Nprev) ~= obj.windowdata.labelidx_cur(1:Nprev));
            if any(idx_relabel),
              obj.SetStatus('Updating fern classifier for %d relabeled examples...',nnz(idx_relabel));
              [obj.classifier] = fernsClfRelabelTrainingData( obj.windowdata.labelidx_cur(waslabeled), ...
                obj.windowdata.labelidx_new(waslabeled), obj.classifier );
              % update labelidx_cur
              obj.windowdata.labelidx_cur(idx_relabel) = obj.windowdata.labelidx_new(idx_relabel);
            end
            
            % remove training examples that were labeled but now aren't
            idx_remove = waslabeled & ~islabeled(1:Nprev);
            if any(idx_remove),
              obj.SetStatus('Removing %d training examples from fern classifier',nnz(idx_remove));
              [obj.classifier] = fernsClfRemoveTrainingData(obj.windowdata.labelidx_cur(waslabeled), ...
                idx_remove(waslabeled), obj.classifier );
              % update labelidx_cur
              obj.windowdata.labelidx_cur(idx_remove) = 0;
            end
            % update islabeled and waslabeled
            islabeled = obj.windowdata.labelidx_new ~= 0;
            waslabeled = [obj.windowdata.labelidx_cur ~= 0;false(Ncurr-Nprev,1)];
            % now only examples with islabeled should be in training set
            
            % add training examples that are labeled now but weren't before
            idx_add = ~waslabeled(islabeled);
            if any(idx_add),
              obj.SetStatus('Adding %d new examples to fern classifier...',nnz(idx_add));
              [obj.classifier] = fernsClfAddTrainingData( obj.windowdata.X(islabeled,:), ...
                obj.windowdata.labelidx_new(islabeled), find(idx_add), obj.classifier );
              % update labelidx_cur
              obj.windowdata.labelidx_cur(~waslabeled&islabeled) = ...
                obj.windowdata.labelidx_new(~waslabeled&islabeled);
            end
            
            % labelidx_cur and new should match
            if ~all(obj.windowdata.labelidx_cur == obj.windowdata.labelidx_new),
              error('Sanity check: labelidx_cur and labelidx_new should match');
            end
          end
          
          obj.classifierTS = now();
          obj.windowdata.isvalidprediction(:) = false;
          obj.windowdata.scoreNorm = [];
          obj.isValidated = false;
          % predict for all window data
          obj.PredictLoaded();
          
        case 'boosting',
          if nargin<2
            doFastUpdates = false;
          end

          if obj.DoFullTraining(doFastUpdates),
            obj.SetStatus('Training boosting classifier from %d examples...',nnz(islabeled));

            obj.classifier_old = obj.classifier;
            [obj.windowdata.binVals] = findThresholds(obj.windowdata.X(islabeled,:),obj.classifier_params);
            bins = findThresholdBins(obj.windowdata.X(islabeled,:),obj.windowdata.binVals);
            [obj.classifier, ~] =...
                boostingWrapper( obj.windowdata.X(islabeled,:), ...
                                 obj.windowdata.labelidx_new(islabeled),obj,...
                                 obj.windowdata.binVals,...
                                 bins,obj.classifier_params);
            obj.lastFullClassifierTrainingSize = nnz(islabeled);
            
          else
            oldNumPts = nnz(obj.windowdata.labelidx_cur ~= 0 & obj.windowdata.labelidx_imp );
            newNumPts = nnz(obj.windowdata.labelidx_new ~= 0 & obj.windowdata.labelidx_imp );
            newData = newNumPts - oldNumPts;
            obj.SetStatus('Updating boosting classifier with %d examples...',newData);
            
            bins = findThresholdBins(obj.windowdata.X(islabeled,:),obj.windowdata.binVals);
            
            obj.classifier_old = obj.classifier;
            [obj.classifier, ~] = boostingUpdate(obj.windowdata.X(islabeled,:),...
                                          obj.windowdata.labelidx_new(islabeled),...
                                          obj.classifier,obj.windowdata.binVals,...
                                          bins,obj.classifier_params);
          end
          obj.classifierTS = now();
          obj.windowdata.labelidx_old = obj.windowdata.labelidx_cur;
          obj.windowdata.labelidx_cur = obj.windowdata.labelidx_new;
          
          obj.predictdata.old = obj.predictdata.cur;
          obj.predictdata.old_valid = obj.predictdata.cur_valid;
          obj.predictdata.old_pp = obj.predictdata.cur_pp;
          obj.predictdata.cur_valid(:) = false;
          
          obj.windowdata.scoreNorm = [];
          % To later find out where each example came from.

%           obj.windowdata.isvalidprediction = false(numel(islabeled),1);
          
          obj.FindFastPredictParams();
          obj.predictdata.cur_valid(:) = false;
          obj.PredictLoaded();
      end

      obj.ClearStatus();
      
    end
    
    function res = DoFullTraining(obj,doFastUpdates)
      % Check if we should do fast updates or not.
      res = true; return; % Always do complete training.
      if ~doFastUpdates, return; end
      if isempty(obj.classifier), return; end
      if isempty(obj.windowdata.binVals), return; end
      
      if (numel(obj.classifier) - obj.classifier_params.iter)/obj.classifier_params.iter_updates > 4
        return;
      end
      
      oldNumPts = obj.lastFullClassifierTrainingSize;
      newNumPts = nnz(obj.windowdata.labelidx_new ~= 0 & obj.windowdata.labelidx_imp );
      newData = newNumPts - oldNumPts;
      if (newData/oldNumPts)>0.25, return; end
      
      res = false;
    end
    
    function PredictLoaded(obj)
    % PredictLoaded(obj)
    % Runs the classifier on all preloaded window data. 
      
      if isempty(obj.classifier),
        return;
      end
      
      % apply classifier
      switch obj.classifiertype,
        
        case 'ferns',
          obj.SetStatus('Applying fern classifier to %d windows',size(obj.windowdata.X,1));
          [obj.windowdata.predicted,...
            obj.windowdata.predicted_probs,...
            obj.predict_cache.last_predicted_inds] = ...
            fernsClfApply(obj.windowdata.X,obj.classifier);
          obj.windowdata.isvalidprediction(:) = true;
          s = exp(obj.windowdata.predicted_probs);
          s = bsxfun(@rdivide,s,sum(s,2));
          scores = max(s,[],2);
          idx0 = obj.windowdata.predicted == 1;
          idx1 = obj.windowdata.predicted > 1;
          obj.windowdata.scores(idx1) = -scores(idx1);
          obj.windowdata.scores(idx0) = scores(idx0);
          obj.ClearStatus();
        case 'boosting',
          if(isempty(obj.predictblocks.t0)), return, end
          
          for ndx = 1:numel(obj.predictblocks.t0)
            curex = obj.predictblocks.expi(ndx);
            flies = obj.predictblocks.flies(ndx);
            numcurex = nnz(obj.predictblocks.expi(:) == curex & ...
              obj.predictblocks.flies(:) == flies);
            numcurexdone = nnz(obj.predictblocks.expi(1:ndx) == curex & ...
              obj.predictblocks.flies(1:ndx) == flies);
            obj.SetStatus('Predicting for exp %s fly %d ... %d%% done',...
            obj.expnames{curex},flies,round(100*numcurexdone/numcurex));
            obj.PredictFast(obj.predictblocks.expi(ndx),...
              obj.predictblocks.flies(ndx),...
              obj.predictblocks.t0(ndx):obj.predictblocks.t1(ndx));
          end
            obj.NormalizeScores([]);
            obj.ApplyPostprocessing();
            obj.ClearStatus();

%{          
%           toPredict = ~obj.windowdata.isvalidprediction;
%           if any(toPredict),
%             obj.SetStatus('Applying boosting classifier to %d windows',sum(toPredict));
%             scores = myBoostClassify(obj.windowdata.X(toPredict,:),obj.classifier);
%             obj.windowdata.predicted(toPredict) = -sign(scores)*0.5+1.5;
%             obj.windowdata.scores(toPredict) = scores;
%             obj.windowdata.isvalidprediction(toPredict) = true;
%             if ~isempty(obj.classifier_old),
%               obj.windowdata.scores_old(toPredict) = ...
%                 myBoostClassify(obj.windowdata.X(toPredict,:),obj.classifier_old);
%             else
%               obj.windowdata.scores_old(toPredict) = 0;
%             end
%             obj.NormalizeScores([]);
%             obj.ApplyPostprocessing();
%             obj.ClearStatus();
%           end
%}
          
      end
            
      % transfer to predictidx for current fly
      if ~isempty(obj.expi) && obj.expi > 0 && ~isempty(obj.flies) && all(obj.flies > 0),
        obj.UpdatePredictedIdx();
      end
      
    end
    
    function SetTrainingData(obj,trainingdata)
    % SetTrainingData(obj,trainingdata)
    % Sets the labelidx_cur of windowdata based on the input training data.
    % This reflects the set of labels the classifier was last trained on. 

      for i = 1:numel(trainingdata),
        [ism,labelidx] = ismember(trainingdata(i).names,obj.labelnames);
        if any(~ism),
          tmp = unique(trainingdata(i).names(~ism));
          error('Unknown labels %s',sprintf('%s ',tmp{:})); %#ok<SPERR>
        end
        isexp = obj.windowdata.exp == i;
        for j = 1:numel(trainingdata(i).t0s),
          t0 = trainingdata(i).t0s(j);
          t1 = trainingdata(i).t1s(j);
          l = labelidx(j);
          flies = trainingdata(i).flies(j,:);
          isflies = isexp & all(bsxfun(@eq,obj.windowdata.flies,flies),2);
          ist = isflies & obj.windowdata.t >= t0 & obj.windowdata.t < t1;
          if nnz(ist) ~= (t1-t0),
            error('Sanity check: number of training examples does not match windowdata');
          end
          obj.windowdata.labelidx_cur(ist) = l;
        end
      end
            
    end

    function trainingdata = SummarizeTrainingData(obj)
    % trainingdata = SummarizeTrainingData(obj)
    % Summarize labelidx_cur into trainingdata, which is similar to the
    % form of labels.
      
      trainingdata = struct('t0s',{},'t1s',{},'names',{},'flies',{});
      waslabeled = obj.windowdata.labelidx_cur;
      for expi = 1:obj.nexps,
        trainingdata(expi) = struct('t0s',[],'t1s',[],'names',{{}},'flies',[]);
        isexp = waslabeled & obj.windowdata.exp == expi;
        if ~any(isexp),
          continue;
        end
        fliess = unique(obj.windowdata.flies(isexp,:),'rows');
        for fliesi = 1:size(fliess,1),
          flies = fliess(fliesi,:);
          isflies = isexp & all(bsxfun(@eq,obj.windowdata.flies,flies),2);
          labelidxs = setdiff(unique(obj.windowdata.labelidx_cur(isflies)),0);
          for labelidxi = 1:numel(labelidxs),
            labelidx = labelidxs(labelidxi);
            islabel = isflies & labelidx == obj.windowdata.labelidx_cur;
            ts = sort(obj.windowdata.t(islabel));
            breaks = find(ts(1:end-1)+1~=ts(2:end));
            t1s = ts(breaks)+1;
            t0s = ts(breaks+1);
            t0s = [ts(1);t0s];%#ok<AGROW>
            t1s = [t1s;ts(end)+1];%#ok<AGROW>
            n = numel(t0s);
            trainingdata(expi).t0s(end+1:end+n,1) = t0s;
            trainingdata(expi).t1s(end+1:end+n,1) = t1s;
            trainingdata(expi).names(end+1:end+n,1) = repmat(obj.labelnames(labelidx),[1,n]);
            trainingdata(expi).flies(end+1:end+n,:) = repmat(flies,[n,1]);
          end
        end
      end

    end

    function UpdatePredictedIdx(obj)
    % UpdatePredictedIdx(obj)
    % Updates the stored predictedidx and erroridx fields to reflect
    % windowdata.predicted
      
      if obj.expi == 0,
        return;
      end
      
      n = obj.t1_curr - obj.t0_curr + 1;
      obj.predictedidx = zeros(1,n);
      obj.scoresidx = zeros(1,n);
      obj.scoresidx_old = zeros(1,n);
      obj.scoreTS = zeros(1,n);
      
      
      if isempty(obj.predictdata.exp),
        return;
      end
      
      % Overwrite by scores from windowdata.
      idxcurr = obj.FlyNdxPredict(obj.expi,obj.flies) & ...
        obj.predictdata.cur_valid;
      obj.predictedidx(obj.predictdata.t(idxcurr)-obj.t0_curr+1) = ...
        -sign(obj.predictdata.cur(idxcurr))*0.5 + 1.5;
      obj.scoresidx(obj.predictdata.t(idxcurr)-obj.t0_curr+1) = ...
        obj.predictdata.cur(idxcurr);      
      obj.scoreTS(obj.predictdata.t(idxcurr)-obj.t0_curr+1) = ...
        obj.classifierTS;      

      obj.UpdateErrorIdx();
            
    end
    
    function UpdateErrorIdx(obj)
    % UpdatePredictedIdx(obj)
    % Updates the stored erroridx and suggestedidx from predictedidx

      if obj.expi == 0,
        return;
      end
      
      n = obj.t1_curr - obj.t0_curr + 1;
      obj.erroridx = zeros(1,n);
      obj.suggestedidx = zeros(1,n);
      idxcurr = obj.predictedidx ~= 0 & obj.labelidx.vals ~= 0;
      obj.erroridx(idxcurr) = double(obj.predictedidx(idxcurr) ~= obj.labelidx.vals(idxcurr))+1;
      
      idxcurr = obj.predictedidx ~= 0 & obj.labelidx.vals == 0;
      obj.suggestedidx(idxcurr) = obj.predictedidx(idxcurr);
    end

    function Predict(obj,expi,flies,ts)
    % Predict(obj,expi,flies,ts)
    % Runs the behavior classifier on the input experiment, flies, and
    % frames. This involves first precomputing the window data for these
    % frames, then applying the classifier. 
     
      % TODO: don't store window data just because predicting. 
      
      if isempty(obj.classifier),
        return;
      end

      if isempty(ts),
        return;
      end
            
      % compute window data
%       [success,msg] = obj.PreLoadWindowData(expi,flies,ts);
%       if ~success,
%         warning(msg);
%         return;
%       end
      
      % indices into windowdata
%       idxcurr = obj.FlyNdx(expi,flies) & ...
%         ~obj.windowdata.isvalidprediction; % & ismember(obj.windowdata.t,ts);
%       
%       if ~any(idxcurr), return; end;
      
      % apply classifier
      switch obj.classifiertype,
        
        case 'ferns',
          obj.SetStatus('Applying fern classifier to %d windows',nnz(idxcurr));
          [obj.windowdata.predicted(idxcurr),...
            obj.windowdata.predicted_probs(idxcurr,:)] = ...
            fernsClfApply(obj.windowdata.X(idxcurr,:),obj.classifier);
          obj.windowdata.isvalidprediction(idxcurr) = true;

          s = exp(obj.windowdata.predicted_probs);
          s = bsxfun(@rdivide,s,sum(s,2));
          scores = max(s,[],2);
          idxcurr1 = find(idxcurr);
          idx0 = obj.windowdata.predicted(idxcurr) == 1;
          idx1 = obj.windowdata.predicted(idxcurr) > 1;
          obj.windowdata.scores(idxcurr1(idx1)) = -scores(idx1);
          obj.windowdata.scores(idxcurr1(idx0)) = scores(idx0);
          
          obj.ClearStatus();
        case 'boosting',

%           obj.SetStatus('Applying boosting classifier to %d frames',nnz(idxcurr));
%           scores = myBoostClassify(obj.windowdata.X(idxcurr,:),obj.classifier);
%           obj.windowdata.predicted(idxcurr) = -sign(scores)*0.5+1.5;
%           obj.windowdata.scores(idxcurr) = scores;
%           if ~isempty(obj.classifier_old)
%             obj.windowdata.scores_old(idxcurr) = ...
%               myBoostClassify(obj.windowdata.X(idxcurr,:),obj.classifier_old);
%           else
%             obj.windowdata.scores_old(idxcurr) = zeros(size(scores));
%           end
%             
%           obj.windowdata.isvalidprediction(idxcurr) = true;

          obj.SetStatus('Updating Predictions ...');
          obj.PredictFast(expi,flies,ts)
          obj.ApplyPostprocessing();
          obj.ClearStatus();

      end
           
      obj.UpdatePredictedIdx();
      
    end
    
    function FindFastPredictParams(obj)
      
        if isempty(obj.classifier),
          return;
        end
        
        feature_names = obj.windowdata.featurenames;
        
        % which features are actually used
        dims = [obj.classifier.dim];
        feature_names = feature_names(dims);

        % put these in with the rest of the classifiers' window features
        wfs = {};
        for j = 1:numel(feature_names),
          wfidxcurr = find(WindowFeatureNameCompare(feature_names{j},wfs),1);
          if isempty(wfidxcurr),
            wfidxcurr = numel(wfs)+1;
            wfs{wfidxcurr} = feature_names{j}; %#ok<AGROW>
          end
        end

        wf2pff = cellfun(@(x)x{1},wfs,'UniformOutput',false);
        [pffs,~,wf2pffidx] = unique(wf2pff);

        windowfeaturescellparams = struct;
        for pfi = 1:numel(pffs),
          pf = pffs{pfi};
          wfidx_cur = wf2pffidx==pfi;
          windowfeaturescellparams.(pf) = WindowFeatureName2Params(wfs(wfidx_cur));
        end

        classifiers_indexed = obj.classifier;
        for j = 1:numel(classifiers_indexed),
          classifiers_indexed(j).dim = j;
        end
        
        obj.fastPredict.classifier = classifiers_indexed;
        obj.fastPredict.windowfeaturescellparams = windowfeaturescellparams;
        obj.fastPredict.wfs = feature_names;
        obj.fastPredict.pffs = pffs;
        obj.fastPredict.ts = obj.classifierTS;
        obj.fastPredict.wfidx_valid = false;
    end
    
    
    function FindWfidx(obj,feature_names)
      
      wfs = obj.fastPredict.wfs;
      wfidx = nan(1,numel(wfs));
      for j = 1:numel(wfs),
        idxcurr = find(WindowFeatureNameCompare(wfs{j},feature_names));
        if numel(idxcurr) ~= 1,
          error('Error matching wfs for classifier with window features computed');
        end
        wfidx(j) = idxcurr;
      end
      obj.fastPredict.wfidx = wfidx;
      obj.fastPredict.wfidx_valid = true;
    end
    
    function PredictFast(obj,expi,flies,ts)
    % Predict fast by computing only the required window features.
      
        if isempty(obj.classifier) || isempty(ts),
          return;
        end
        
        if isempty(obj.fastPredict.classifier)
          obj.FindFastPredictParams();
        end
       
        idxcurr = obj.predictdata.flies == flies & ...
          obj.predictdata.exp == expi &...
          obj.predictdata.cur_valid;
        tscurr = obj.predictdata.t(idxcurr);
        missingts = setdiff(ts,tscurr);
        if numel(missingts)==0, return; end
        
%         if isempty(missingts), return; end
        
        perframeInMemory = ~isempty(obj.flies) && obj.IsCurFly(expi,flies);
        perframefile = obj.GetPerframeFiles(expi);

       while true,

          % choose a frame missing window data
          % if a predict block around that has already been used then use
          % that.
          
          t = missingts(1);
          curbs_t0 = obj.predictblocks.t0( obj.predictblocks.expi==expi & ...
            obj.predictblocks.flies == flies);
          curbs_t1 = obj.predictblocks.t1(obj.predictblocks.expi == expi & ...
            obj.predictblocks.flies == flies);
          if ~isempty(curbs_t0) && any(  (t-curbs_t0)>=0 & (t-curbs_t1)<=0)
            tempndx = find( (t-curbs_t0)>=0 & (t-curbs_t1)<=0 );
            t0 = curbs_t0(tempndx(1));
            t1 = curbs_t1(tempndx(1));
          else
            
            t = median(missingts);
            if ~ismember(t,missingts),
              t = missingts(argmin(abs(t-missingts)));
            end
            
            % bound at start and end frame of these flies
            T0 = max(obj.GetTrxFirstFrame(expi,flies));
            T1 = min(obj.GetTrxEndFrame(expi,flies));
            
            t1 = min(t+obj.predictwindowdatachunk_radius,T1);
            % go backward 2*r to find the start of the chunk
            t0 = max(t1-2*obj.predictwindowdatachunk_radius,T0);
            % go forward 2*r again to find the end of the chunk
            t1 = min(t0+2*obj.predictwindowdatachunk_radius,T1);
            
            overlap_start = find( (t0-curbs_t0)>=0 & (t0-curbs_t1)<=0);
            if ~isempty(overlap_start),
              t0 = max(curbs_t1(overlap_start))+1;
            end
            
            overlap_end = find( (t1-curbs_t0)>=0 & (t1-curbs_t1)<=0);
            if ~isempty(overlap_end),
              t1 = max(curbs_t0(overlap_end))-1;
            end
            
            obj.predictblocks.t0(end+1) = t0;
            obj.predictblocks.t1(end+1) = t1;
            obj.predictblocks.expi(end+1) = expi;
            obj.predictblocks.flies(end+1) = flies;
          end
          
          i0 = t0 - obj.GetTrxFirstFrame(expi,flies) + 1;
          i1 = t1 - obj.GetTrxFirstFrame(expi,flies) + 1;
          
          X = [];
          
          % for the parfor loop.
          perframedata_cur = obj.perframedata;
          windowfeaturescellparams = obj.fastPredict.windowfeaturescellparams;
          pffs = obj.fastPredict.pffs;
          allperframefns = obj.allperframefns;
          feature_names_list = cell(1,numel(pffs));
          x_curr_all = cell(1,numel(pffs));
          
          parfor j = 1:numel(pffs),

            fn = pffs{j};
            
            ndx = find(strcmp(fn,allperframefns));
            if perframeInMemory,
              perframedata = perframedata_cur{ndx};
            else
              perframedata = load(perframefile{ndx});
              perframedata = perframedata.data{flies(1)};
            end
            
            i11 = min(i1,numel(perframedata));
            [x_curr,cur_f] = ...
              ComputeWindowFeatures(perframedata,...
              windowfeaturescellparams.(fn){:},'t0',i0,'t1',i11);
            
            if i11 < i1,
              x_curr(:,end+1:end+i1-i11) = nan;
            end
            
            x_curr_all{j} = single(x_curr);
            feature_names_list{j} = cur_f;
          end
          
          
          for j = 1:numel(pffs),
            fn = pffs{j};
            x_curr = x_curr_all{j};
            % add the window data for this per-frame feature to X
            nold = size(X,1);
            nnew = size(x_curr,2);
            if nold > nnew,
              warning(['Number of examples for per-frame feature %s does not '...
                'match number of examples for previous features'],fn);
              x_curr(:,end+1:end+nold-nnew) = nan;
            elseif nnew > nold && ~isempty(X),
              warning(['Number of examples for per-frame feature %s does not '...
                'match number of examples for previous features'],fn);
              X(end+1:end+nnew-nold,:) = nan;
            end
            X = [X,x_curr']; %#ok<AGROW>
          end
         
          if ~obj.fastPredict.wfidx_valid,
            feature_names = {};
            for ndx = 1:numel(feature_names_list)
              fn = pffs{ndx};
              feature_names = [feature_names,cellfun(@(s) [{fn},s],...
                feature_names_list{ndx},'UniformOutput',false)]; %#ok<AGROW>
            end
            obj.FindWfidx(feature_names);
          end
          
          scores = myBoostClassify(X(:,obj.fastPredict.wfidx),obj.fastPredict.classifier);
          
          curndx = obj.FlyNdxPredict(expi,flies) & ...
            obj.predictdata.t>=t0 & ...
            obj.predictdata.t<=t1;
          
          obj.predictdata.cur(curndx) = scores;
          obj.predictdata.timestamp(curndx) = obj.classifierTS;
          obj.predictdata.cur_valid(curndx) = true;
          
          missingts(missingts >= t0 & missingts <= t1) = [];

          if isempty(missingts),
            obj.ClearStatus();
            break;
          end
        end
        
 
    end
    
    
    function allScores = PredictWholeMovie(obj,expi)
      
      if isempty(obj.classifier),
        return;
      end
           
      numFlies = obj.GetNumFlies(expi);
      scoresA = cell(1,numFlies); 
      postprocessedscoresA = cell(1,numFlies);
      tStartAll = obj.GetTrxFirstFrame(expi);
      tEndAll = obj.GetTrxEndFrame(expi);
      perframefile = obj.GetPerframeFiles(expi);
      windowfeaturescellparams = obj.fastPredict.windowfeaturescellparams;
      curperframefns = obj.fastPredict.pffs;
      allperframefns = obj.allperframefns;
      classifier = obj.fastPredict.classifier;
      
      obj.SetStatus('Classifying current movie..');
      
      if ~obj.fastPredict.wfidx_valid,
        [~,feature_names] = JLabelData.ComputeWindowDataChunkStatic(curperframefns,...
          allperframefns,perframefile,1,windowfeaturescellparams,...
          1,1);
        obj.FindWfidx(feature_names);
      end
      wfidx = obj.fastPredict.wfidx;
      
      
      parfor flies = 1:numFlies
        blockSize = 5000*2;
        tStart = tStartAll(flies);
        tEnd = tEndAll(flies);
        
        scores = nan(1,tEnd);
        
        for curt0 = tStart:blockSize:tEnd
          curt1 = min(curt0+blockSize-1,tEnd);
          X = JLabelData.ComputeWindowDataChunkStatic(curperframefns,...
            allperframefns,perframefile,flies,windowfeaturescellparams,curt0-tStart+1,curt1-tStart+1);
        
          scores(curt0:curt1) = myBoostClassify(X(:,wfidx),classifier);
        end
        scoresA{flies} = scores;
        fprintf('Prediction done for %d fly, total number of flies:%d\n',flies,numFlies);
      end
      
      allScores = struct;
      for flies = 1:numFlies
        postprocessedscoresA{flies} = nan(1,tEndAll(flies));
        postprocessedscoresA{flies}(tStartAll(flies):tEndAll(flies)) = ...
          obj.Postprocess(scoresA{flies}(tStartAll(flies):tEndAll(flies)));
      end
      allScores.scores = scoresA;
      allScores.tStart = tStartAll;
      allScores.tEnd = tEndAll;
      allScores.postprocessed = postprocessedscoresA;
      allScores.postprocessparams = obj.postprocessparams;
      
      for flies = 1:numFlies
        [i0s i1s] = get_interval_ends(allScores.scores{flies}>0);
        allScores.t0s{flies} = i0s;
        allScores.t1s{flies} = i1s;
      end
      
      obj.ClearStatus();
      
    end
    
    function PredictSaveMovie(obj,expi,sfn)
    % Predicts for the whole movie and saves the scores.
      if nargin < 3
        sfn = obj.GetFile('scores',expi);
      end
      allScores = obj.PredictWholeMovie(expi);
      obj.SaveScores(allScores,expi,sfn);
      obj.AddScores(expi,allScores,now(),'',true);
      
      idxexp = obj.predictdata.exp == expi;
      if any(obj.predictdata.loaded_valid(idxexp)),
        obj.LoadScores(expi,sfn);
      end
    end
    
    function PredictNoSaveMovie(obj,expi)
      allScores = obj.PredictWholeMovie(expi);
      obj.AddScores(expi,allScores,now(),'',true);
    end
    
% Evaluating performance
    

    function errorRates = createConfMat(obj,scores,modLabels)
      
      confMat = zeros(2*obj.nbehaviors,3);
      scoreNorm = obj.windowdata.scoreNorm;
      for ndx = 1:2*obj.nbehaviors
        if mod(ndx,2)
          curIdx = modLabels==ndx;
        else
          curIdx = modLabels > (ndx-1.5) & modLabels<(ndx+0.5);
        end
        confMat(ndx,1) = nnz(scores(curIdx)>=  (obj.confThresholds(1)*scoreNorm));
        confMat(ndx,2) = nnz(-scores(curIdx)<  (obj.confThresholds(2)*scoreNorm) & ...
                              scores(curIdx)<  (obj.confThresholds(1)*scoreNorm) );
        confMat(ndx,3) = nnz(-scores(curIdx)>= (obj.confThresholds(2)*scoreNorm));
      end
      errorRates.numbers = confMat;
      errorRates.frac = errorRates.numbers./repmat( sum(errorRates.numbers,2),[1 3]);
    end
    
    function bouts = getLabeledBouts(obj)
    % Find the bouts from window data.

    bouts = struct('ndx',[],'label',[],'timestamp',[]);
      for expNdx = 1:obj.nexps
        for flyNdx = 1:obj.nflies_per_exp(expNdx)
          curLabels = obj.GetLabels(expNdx,flyNdx);
          for boutNum = 1:numel(curLabels.t0s)
            idx =  obj.FlyNdx(expNdx,flyNdx) & ...
              obj.windowdata.t >= curLabels.t0s(boutNum) & ...
              obj.windowdata.t < curLabels.t1s(boutNum);
            if ~all(obj.windowdata.labelidx_cur(idx)), continue; end
            bouts.ndx(end+1,:) = obj.FlyNdx(expNdx,flyNdx) & ...
              obj.windowdata.t >= curLabels.t0s(boutNum) & ...
              obj.windowdata.t < curLabels.t1s(boutNum);
            bouts.label(end+1) = find(strcmp(obj.labelnames,curLabels.names{boutNum}));
            bouts.timestamp(end+1) = curLabels.timestamp(boutNum);
          end
          
        end
      end
      
    end
    
    function [success,msg,crossError,tlabels] = CrossValidate(obj,varargin)
    % Cross validate on bouts.
      
      obj.StoreLabels();
      
      [success,msg] = obj.PreLoadLabeledData();
      if ~success, 
        return;
      end
      
      [setidx,byexp] = myparse(varargin,'setidx',[],'byexp',false);

      islabeled = obj.windowdata.labelidx_cur ~= 0;
      if ~any(islabeled),                        
        crossError.numbers = zeros(4,3);
        crossError.frac = zeros(4,3);
        crossError.oldNumbers = zeros(4,3);
        crossError.oldFrac = zeros(4,3);
        tlabels = {};
        success = false;
        msg = 'No Labeled Data';
        obj.ClearStatus();
        return; 
      end
      
      if ~strcmp(obj.classifiertype,'boosting'); return; end

      obj.SetStatus('Cross validating the classifier for %d examples...',nnz(islabeled));

      obj.UpdateBoostingBins();

      bouts = obj.getLabeledBouts();
      
      if byexp && isempty(setidx),
        [~,~,setidx] = unique(obj.windowdata.exp);
      end
      
      [success,msg,crossScores, tlabels]=...
        crossValidateBout( obj.windowdata.X, ...
        obj.windowdata.labelidx_cur,bouts,obj,...
        obj.windowdata.binVals,...
        obj.classifier_params,true,setidx);
      
      if ~success, 
        crossError.numbers = zeros(4,3);
        crossError.frac = zeros(4,3);
        crossError.oldNumbers = zeros(4,3);
        crossError.oldFrac = zeros(4,3);
        tlabels = {};
        obj.ClearStatus();
        return, 
      end;

%{      
%       crossScores=...
%         crossValidate( obj.windowdata.X(islabeled,:), ...
%         obj.windowdata.labelidx_cur(islabeled,:),obj,...
%         obj.windowdata.binVals,...
%         obj.windowdata.bins(:,islabeled),obj.classifier_params);
%}
      
      obj.windowdata.scores_validated = zeros(numel(islabeled),1);
      obj.windowdata.scores_validated(islabeled) = crossScores(1,:);

      modLabels = 2*obj.windowdata.labelidx_cur(islabeled)-obj.windowdata.labelidx_imp(islabeled);

      for tndx = 1:size(crossScores,1)
        crossError(tndx) = obj.createConfMat(crossScores(tndx,:),modLabels);
      end
      
      waslabeled = false(numel(islabeled),1);
      waslabeled(1:numel(obj.windowdata.labelidx_old)) = obj.windowdata.labelidx_old~=0;
      oldSelect = waslabeled(islabeled);
      oldScores = crossScores(oldSelect);
      oldLabels = 2*obj.windowdata.labelidx_cur(waslabeled(:)&islabeled(:)) - ...
          obj.windowdata.labelidx_imp(waslabeled(:)&islabeled(:));
      oldError = obj.createConfMat(oldScores,oldLabels);
      crossError(1).oldNumbers = oldError.numbers;
      crossError(1).oldFrac = oldError.frac;
      
      obj.ClearStatus();
    end

    
    function ROCCurve(obj,JLabelHandle)
    % This now shows histogram, apt naming be damned.
    
      curNdx = obj.windowdata.labelidx_cur~=0;
      curScores = obj.windowdata.scores(curNdx);
      curLabels = obj.windowdata.labelidx_cur(curNdx);
      modLabels = ((curLabels==1)-0.5)*2;
      ShowROCCurve(modLabels,curScores,obj,JLabelHandle);
      
    end
    
    function newError = TestOnNewLabels(obj)
      obj.StoreLabels();
      newError = struct();

      prevLabeled = obj.windowdata.labelidx_cur~=0;
      Nprev = numel(prevLabeled);
      newLabels = obj.windowdata.labelidx_new ~= 0;
      tOld = newLabels(1:Nprev);
      tOld(prevLabeled) = false;
      newLabels(1:Nprev) = tOld;
      
      if ~nnz(newLabels); 
        fprintf('No new labeled data\n');
        return;
      end
      
      % Find out the index of scores with the same exp, flynum and time as
      % the newly labeled data.
      
      orderedScores = []; orderedLabels = []; orderedLabels_imp = [];
      nlexp = obj.windowdata.exp(newLabels);
      nlflies = obj.windowdata.flies(newLabels);
      nlt = obj.windowdata.t(newLabels);
      nlLabels = obj.windowdata.labelidx_new(newLabels);
      nlLabels_imp = obj.windowdata.labelidx_imp(newLabels);
      
      classifierfilename = 'None'; setClassifierfilename = 1;
      for curExp = unique(nlexp)'
        curNLexpNdx = nlexp==curExp;
        for curFly = unique(nlflies(curNLexpNdx))';
          curT = nlt( nlexp==curExp & nlflies == curFly);
          curLabels = nlLabels(nlexp==curExp & nlflies == curFly);
          curLabels_imp = nlLabels_imp(nlexp==curExp & nlflies == curFly);
          curScoreNdx = find(obj.predictdata.exp == curExp & obj.predictdata.flies==curFly);
          scoresT = obj.predictdata.t(curScoreNdx);
          [curValidScoreNdx,loc] = ismember(scoresT,curT);
          if nnz(curValidScoreNdx)~=numel(curT)
            warndlg('Scores are missing for some labeled data');
            newError = struct();
            return;
          end
          
          orderedLabels = [orderedLabels; curLabels(loc(loc~=0))];
          orderedLabels_imp = [orderedLabels_imp; curLabels_imp(loc(loc~=0))];
          orderedScores = [orderedScores; obj.predictdata.loaded(curScoreNdx(curValidScoreNdx~=0))'];
        end
%         if setClassifierfilename,
%           classifierfilename = obj.windowdata.classifierfilenames{curExp};
%           setClassifierfilename = 0;
%         elseif strcmp(classifierfilename,'multiple'),
%         elseif ~strcmp(classifierfilename,obj.windowdata.classifierfilenames{curExp}),
%           classifierfilename = 'multiple';
%         end
          
      end
      
      prediction = -sign(orderedScores)/2+1.5;
      
      modLabels = 2*orderedLabels-orderedLabels_imp;
      
      newError = obj.createConfMat(prediction,modLabels);
      newError.classifierfilename = classifierfilename;
      
    end
    
    
% Show similar frames
    

    function DoBagging(obj)

      obj.StoreLabels();
      [success,msg] = obj.PreLoadLabeledData();
      if ~success, warning(msg);return;end

      islabeled = obj.windowdata.labelidx_new ~= 0;

      if ~any(islabeled),                        return; end
      if ~strcmp(obj.classifiertype,'boosting'); return; end
      if isempty(obj.classifier), obj.Train;             end

      bouts = struct('ndx',[],'label',[],'timestamp',[]);
      for expNdx = 1:obj.nexps
        for flyNdx = 1:obj.nflies_per_exp(expNdx)
          curLabels = obj.GetLabels(expNdx,flyNdx);
          for boutNum = 1:numel(curLabels.t0s)
            bouts.ndx(end+1,:) = obj.FlyNdx(expNdx,flyNdx) & ...
              obj.windowdata.t >= curLabels.t0s(boutNum) & ...
              obj.windowdata.t < curLabels.t1s(boutNum);
            bouts.label(end+1) = find(strcmp(obj.labelnames,curLabels.names{boutNum}));
            bouts.timestamp(end+1) = curLabels.timestamp(boutNum);
          end
          
        end
      end

      obj.SetStatus('Bagging the classifier with %d examples...',nnz(islabeled));
      
      obj.windowdata.binVals = findThresholds(obj.windowdata.X(islabeled,:),obj.classifier_params);
      
      [obj.bagModels, obj.distMat] =...
        doBaggingBouts( obj.windowdata.X, ...
        obj.windowdata.labelidx_new,obj,...
        obj.windowdata.binVals,...
        obj.classifier_params,bouts);
      
      obj.windowdata.distNdx.exp = obj.windowdata.exp(islabeled);
      obj.windowdata.distNdx.flies = obj.windowdata.flies(islabeled);
      obj.windowdata.distNdx.t = obj.windowdata.t(islabeled);
      obj.windowdata.distNdx.labels = obj.windowdata.labelidx_new(islabeled);
      
      obj.ClearStatus();
    end
    
    function InitSimilarFrames(obj, HJLabel)
      obj.frameFig = showSimilarFrames;
      showSimilarFrames('SetJLabelData',obj.frameFig,obj,HJLabel);
      showSimilarFrames('CacheTracksLabeled',obj.frameFig);
      showSimilarFrames('add_prep_list', obj.frameFig);
    end
    
    function SimilarFrames(obj,curTime,JLabelHandles)
      if isempty(obj.frameFig) || ~ishandle(obj.frameFig),
          obj.InitSimilarFrames(JLabelHandles), 
      end
      
      distNdx = find( (obj.windowdata.distNdx.exp == obj.expi) & ...
        (obj.windowdata.distNdx.flies == obj.flies) & ...
        (obj.windowdata.distNdx.t == curTime) ,1);
      
      windowNdx = find( (obj.windowdata.exp == obj.expi) & ...
        (obj.windowdata.flies == obj.flies) & ...
        (obj.windowdata.t == curTime) ,1);


      if isempty(distNdx) % The example was not part of the training data.
        outOfTraining = 1;
        curX = obj.windowdata.X(windowNdx,:);
        curD = zeros(1,length(obj.bagModels)*length(obj.bagModels{1}));
        count = 1;
        for bagNo = 1:length(obj.bagModels)
          curModel = obj.bagModels{bagNo};
          for j = 1:length(curModel)
            curWk = curModel(j);
            dd = curX(curWk.dim)*curWk.dir;
            tt = curWk.tr*curWk.dir;
            curD(count) = (dd>tt)*curWk.alpha;
            count = count+1;
          end
        end
      else
        outOfTraining = 0;
        curD = obj.distMat(distNdx,:);
      end

      % Compute the distance 
      diffMat = zeros(size(obj.distMat));
      for ndx = 1:size(diffMat,2);
        diffMat(:,ndx) = abs(obj.distMat(:,ndx)-curD(ndx));
      end
      dist2train = nanmean(diffMat,2)*200;
      [rr rrNdx] = sort(dist2train,'ascend');
      
      if~outOfTraining
        rr = rr(2:end);
        curEx = rrNdx(1); rrNdx = rrNdx(2:end);
      else
        curEx = [];
      end
      
      % Find 5 closest pos and neg examples.
      % This looks complicated then it should be.
      % DEBUG: find values of actual labels 
     
      trainLabels =  obj.windowdata.distNdx.labels;
      allPos = rrNdx(trainLabels(rrNdx)>1.5);
      allNeg = rrNdx(trainLabels(rrNdx)<1.5);
      
      
      curP = zeros(1,5);
      curN = zeros(1,5);
      count = 0;
      for ex = allPos'
        if count>4; break; end;
        isClose = 0;
        if obj.windowdata.exp(windowNdx) == obj.windowdata.distNdx.exp(ex) &&...
           obj.windowdata.flies(windowNdx) == obj.windowdata.distNdx.flies(ex) && ...
           abs( obj.windowdata.t(windowNdx) - obj.windowdata.distNdx.t(ex) )<5,
           continue; 
        end
        
        for used = curP(1:count)
          if obj.windowdata.distNdx.exp(used) == obj.windowdata.distNdx.exp(ex) &&...
             obj.windowdata.distNdx.flies(used) == obj.windowdata.distNdx.flies(ex) && ...
             abs( obj.windowdata.distNdx.t(used) - obj.windowdata.distNdx.t(ex) )<5,
             isClose = 1; 
             break; 
          end
        end
        
        if isClose; continue; end;
        count = count+1;
        curP(count) = ex;
      end
      
      count = 0;
      for ex = allNeg'
        if count>4; break; end;
        isClose = 0;
        if obj.windowdata.exp(windowNdx) == obj.windowdata.distNdx.exp(ex) &&...
           obj.windowdata.flies(windowNdx) == obj.windowdata.distNdx.flies(ex) && ...
           abs(obj.windowdata.t(windowNdx) - obj.windowdata.distNdx.t(ex))<5,
           continue; 
        end
        
        for used = curN(1:count)
          if obj.windowdata.distNdx.exp(used) == obj.windowdata.distNdx.exp(ex) &&...
             obj.windowdata.distNdx.flies(used) == obj.windowdata.distNdx.flies(ex) && ...
             abs(obj.windowdata.distNdx.t(used) - obj.windowdata.distNdx.t(ex))<5,
             isClose = 1; 
             break; 
          end
        end
        
        if isClose; continue; end;
        count = count+1;
        curN(count) = ex;
      end
      
      varForSSF.curFrame.expNum = obj.windowdata.exp(windowNdx);
      varForSSF.curFrame.flyNum = obj.windowdata.flies(windowNdx);
      varForSSF.curFrame.curTime = obj.windowdata.t(windowNdx);
      
      for k = 1:4
        varForSSF.posFrames(k).expNum = obj.windowdata.distNdx.exp(curP(k));
        varForSSF.posFrames(k).flyNum = obj.windowdata.distNdx.flies(curP(k));
        varForSSF.posFrames(k).curTime = obj.windowdata.distNdx.t(curP(k));
        varForSSF.negFrames(k).expNum = obj.windowdata.distNdx.exp(curN(k));
        varForSSF.negFrames(k).flyNum = obj.windowdata.distNdx.flies(curN(k));
        varForSSF.negFrames(k).curTime = obj.windowdata.distNdx.t(curN(k));
      end
      showSimilarFrames('setFrames',obj.frameFig,varForSSF);
    end
    
    
% Fly and exp statistics
    

    function expStats = GetExpStats(obj,expi)
      % Calculates statistics such as number of labeled bouts, predicted bouts
      % and change in scores.
      
      expStats.name = obj.expnames{expi};
      expStats.nflies = obj.nflies_per_exp(expi);
      expStats.nlabeledbouts = obj.labelstats(expi).nbouts_labeled;
      expStats.nlabeledflies = obj.labelstats(expi).nflies_labeled;
      
      
      if ~isempty(obj.predictdata.exp==expi)
        expid = obj.predictdata.exp==expi;
        expStats.nscoreframes = nnz(expid);
        expStats.nscorepos = nnz(obj.predictdata.loaded(expid)>0);
%         if ~isempty(obj.predictdata.classifierfilenames) && ...
%             numel(obj.predictdata.classifierfilenames)>=expi
%           expStats.classifierfilename = obj.predictdata.classifierfilenames{expi};
%         else
%           expStats.classifierfilename = '';
%         end
      else
        expStats.nscoreframes = [];
        expStats.nscorefrac = [];
        expStats.classifierfilename = '';
      end
      
    end

    function flyStats = GetFlyStats(obj,expi,flyNum)
      % Calculates statistics such as number of labeled bouts, predicted bouts
      % and change in scores.
      obj.SetStatus('Computing stats for %s %d',obj.expnames{expi},flyNum);
      
      obj.StoreLabels();
      [ism,j] = ismember(flyNum,obj.labels(expi).flies,'rows');
      if ism,
        flyStats.nbouts = numel(obj.labels(expi).t0s{j});
        posframes = 0; negframes = 0;
        for ndx = 1:numel(obj.labels(expi).t0s{j})
          numFrames = obj.labels(expi).t1s{j}(ndx)-obj.labels(expi).t0s{j}(ndx);
          if strcmp(obj.labels(expi).names{j}{ndx},obj.labelnames{1}) 
            posframes = posframes + numFrames;
          else
            negframes = negframes + numFrames;
          end
        end
        flyStats.posframes = posframes;
        flyStats.negframes = negframes;
        flyStats.totalframes = posframes + negframes;
      else
        flyStats.nbouts = 0;
        flyStats.posframes = 0;
        flyStats.negframes = 0;
        flyStats.totalframes = 0;
      end
      
      [ism,j] = ismember(flyNum,obj.gt_labels(expi).flies,'rows');
      if ism,
        flyStats.gt_nbouts = numel(obj.gt_labels(expi).t0s{j});
        posframes = 0; negframes = 0;
        for ndx = 1:numel(obj.gt_labels(expi).t0s{j})
          numFrames = obj.gt_labels(expi).t1s{j}(ndx)-obj.gt_labels(expi).t0s{j}(ndx);
          if strcmp(obj.gt_labels(expi).names{j}{ndx},obj.labelnames{1}) 
            posframes = posframes + numFrames;
          else
            negframes = negframes + numFrames;
          end
        end
        flyStats.gt_posframes = posframes;
        flyStats.gt_negframes = negframes;
        flyStats.gt_totalframes = posframes + negframes;
      else
        flyStats.gt_nbouts = 0;
        flyStats.gt_posframes = 0;
        flyStats.gt_negframes = 0;
        flyStats.gt_totalframes = 0;
      end
      
      flyStats.endframe = obj.endframes_per_exp{expi}(flyNum);
      flyStats.firstframe = obj.firstframes_per_exp{expi}(flyNum);
      flyStats.trajLength = flyStats.endframe-flyStats.firstframe+1;
      
      if obj.hassex,
        if obj.hasperframesex,
          sexfrac = obj.GetSexFrac(expi,flyNum);
          flyStats.sexfrac = round(100*sexfrac.M);
        else
          flyStats.sexfrac = 100*strcmpi(obj.GetSex(expi,flyNum),'M');
        end
      else
        flyStats.sexfrac = [];
      end
      
      if ~isempty(obj.predictdata.exp==expi)
        idxcurr = obj.FlyNdxPredict(expi,flyNum) & obj.predictdata.loaded_valid;
        flyStats.nscoreframes_loaded = nnz(idxcurr);
        flyStats.nscorepos_loaded = nnz(obj.predictdata.loaded(idxcurr)>0);
        flyStats.nscoreneg_loaded = nnz(obj.predictdata.loaded(idxcurr)<0);
%         if ~isempty(obj.predictdata.classifierfilenames)
%           flyStats.classifierfilename = obj.predictdata.classifierfilenames{expi};
%         else
%           flyStats.classifierfilename = '';
%         end
      else
        flyStats.nscoreframes_loaded = [];
        flyStats.nscorepos_loaded = [];
        flyStats.nscoreneg_loaded = [];        
%         flyStats.classifierfilename = '';
      end
      
      if ~isempty(obj.predictdata.exp)
        curNdx = obj.FlyNdxPredict(expi,flyNum)& obj.predictdata.cur_valid;
        curWNdx = obj.FlyNdx(expi,flyNum);
      else
        curNdx = [];
      end
      
      if any(curNdx) && ~isempty(obj.classifier)
        curScores = obj.predictdata.cur(curNdx);

        if ~isempty(curWNdx)
          curWScores = myBoostClassify(obj.windowdata.X(curWNdx,:),obj.classifier);
        else
          curWScores = [];
        end
        curLabels = obj.windowdata.labelidx_new(curWNdx);
        
        curPosMistakes = nnz( curWScores<0 & curLabels ==1 );
        curNegMistakes = nnz( curWScores>0 & curLabels >1 );

        flyStats.nscoreframes = nnz(curNdx);
        flyStats.nscorepos = nnz(curScores>0);
        flyStats.nscoreneg = nnz(curScores<0);
        flyStats.errorsPos = curPosMistakes;
        flyStats.errorsNeg = curNegMistakes;
      else
        flyStats.nscoreframes = [];
        flyStats.nscorepos = [];
        flyStats.nscoreneg = [];
        flyStats.errorsPos = [];
        flyStats.errorsNeg = [];
      end
      
      flyStats.one2two = [];
      flyStats.two2one = [];
      if ~isempty(obj.classifier_old),
        curNdx = obj.FlyNdx(expi,flyNum);
        if nnz(curNdx);
          flyStats.one2two = nnz(obj.predictdata.cur(curNdx)<0 ...
            & obj.predictdata.old(curNdx)>0);
          flyStats.two2one = nnz(obj.predictdata.cur(curNdx)>0 ...
            & obj.predictdata.old(curNdx)<0);
        end
      end
      
      flyStats.validatedErrorsPos = [];
      flyStats.validatedErrorsNeg = [];
      if ~isempty(obj.windowdata.scores_validated),
        curNdx = obj.FlyNdx(expi,flyNum);
        if nnz(curNdx);
          curScores = obj.windowdata.scores_validated(curNdx);
          curLabels = obj.windowdata.labelidx_cur(curNdx);
          
          curPosMistakes = nnz( curScores(:)<0 & curLabels(:) ==1 );
          curNegMistakes = nnz( curScores(:)>0 & curLabels(:) >1 );
          
          flyStats.validatedErrorsPos = curPosMistakes;
          flyStats.validatedErrorsNeg = curNegMistakes;
        end
      end
      
      flyStats.gt_suggestion_frames = nnz(obj.GetGTSuggestionIdx(expi,flyNum));
      
      %       if ~isempty(obj.windowdata.X)
      %         idxcurr = obj.windowdata.exp==expi & obj.windowdata.flies == flyNum;
      %         flyStats.npredictframes = nnz(idxcurr);
      %         flyStats.npredictfrac = nnz(obj.windowdata.scores(idxcurr)>0)/flyStats.nscoreframes;
      %
      %       else
      %         flyStats.npredictframes = [];
      %         flyStats.npredictfrac = [];
      %       end
      
      obj.ClearStatus();
   end
    
    function scores = NormalizeScores(obj,scores)
      
      if isempty(obj.windowdata.scoreNorm) || isnan(obj.windowdata.scoreNorm)
        if ~any(obj.predictdata.cur_valid), return; end
        
        idxcurr = obj.predictdata.cur_valid;
        wScores = obj.predictdata.cur(idxcurr);
        scoreNorm = prctile(abs(wScores),80);
        obj.windowdata.scoreNorm = scoreNorm;
      end
      
      scoreNorm = obj.windowdata.scoreNorm;
      scores(scores<-scoreNorm) = -scoreNorm;
      scores(scores>scoreNorm) = scoreNorm;
      scores = scores/scoreNorm;
    end
    
    
% Status display
    

    function SetStatus(obj,varargin)
    % SetStatus(obj,<sprintf-like arguments>)
    % Update an associated status text according to the input sprintf-like
    % arguments.
      
      if isempty(obj.setstatusfn),
        fprintf(varargin{:});
        fprintf('\n');
      else
        obj.setstatusfn(sprintf(varargin{:}));
        drawnow;
      end
      
    end
    
    function ClearStatus(obj)
    % ClearStatus(obj)
    % Return an associated status text to the default.
      
      if ~isempty(obj.clearstatusfn),
        obj.clearstatusfn();
        drawnow;
      end
      
    end
    
    function SetStatusFn(obj,statusfn)
      obj.setstatusfn = statusfn;
    end
    
    function SetClearStatusFn(obj,clearfn)
      obj.clearstatusfn = clearfn;
    end

    
% Ground truthing functions.    


    function [success,msg] = SuggestRandomGT(obj,perfly,perexp)
      
      success = false; msg = '';
      
      % Do nothing if we already have suggestions with the same settings for
      % all the experiments.
      recompute = false;
      if numel(obj.randomGTSuggestions)<numel(obj.nflies_per_exp), 
        recompute = true; 
      else
        for endx = 1:obj.nexps
          prevperexp = 0;
          
          for fndx = 1:obj.nflies_per_exp(endx)
            
            if isempty(obj.randomGTSuggestions{endx}(fndx).start); continue; end
            prevperfly = obj.randomGTSuggestions{endx}(fndx).end - ...
              obj.randomGTSuggestions{endx}(fndx).start+1;
            if (prevperfly ~= perfly),
              recompute = true;
              break;
            end
            prevperexp = prevperexp+1;
          end
          
          if prevperexp ~= perexp, recompute = true; end
          
        end
      end
      if ~recompute;
        obj.GTSuggestionMode = 'Random';
        return;
      end
      
      
      for endx = 1:obj.nexps
        obj.randomGTSuggestions{endx} = repmat(struct('start',[],'end',[]),1,perexp);
        
        validflies = find( (obj.endframes_per_exp{endx} - ...
          obj.firstframes_per_exp{endx})>perfly );
        if numel(validflies)<perexp,
          msg = sprintf('Experiment %s does not have %d flies with more than %d frames',...
            obj.expnames{endx},perexp,perfly);
          success = false;
          return;
        end
        permuteValid = validflies(randperm(numel(validflies)));
        randFlies = permuteValid(1:perexp);
        
        for fndx = randFlies(:)',
          first = obj.firstframes_per_exp{endx}(fndx);
          last = obj.endframes_per_exp{endx}(fndx);
          suggestStart = first + round( (last-first-perfly)*rand(1));
          obj.randomGTSuggestions{endx}(fndx).start = suggestStart;
          obj.randomGTSuggestions{endx}(fndx).end = suggestStart+perfly-1;
        end
        
      end
      success = true;
      obj.GTSuggestionMode = 'Random';
      
    end
    
    function [success,msg] = SuggestBalancedGT(obj,intsize,numint)
      success = true; msg = '';
      
      if ~any(obj.predictdata.loaded_valid),
        msg = 'No Scores have been loaded, cannot suggest intervals for ground truthing\n';
        success = false;
        return;
      end
      
      numpos = nnz(obj.predictdata.loaded>0);
      numneg = nnz(obj.predictdata.loaded<0);
      poswt = numneg/(numneg+numpos);
      negwt = numpos/(numneg+numpos);
      
      int = struct('exp',[],'flies',[],'tStart',[],'wt',[]);
      obj.balancedGTSuggestions = {};
      for endx = 1:obj.nexps
        for flies = 1:obj.nflies_per_exp(endx)
          curidx = obj.predictdata.exp == endx & obj.predictdata.flies == flies;
          curt = obj.predictdata.t(curidx);
          if numel(curt)<intsize; continue; end
          if any(curt(2:end)-curt(1:end-1) ~= 1)
            msg = 'Scores are not in order'; 
            success = false; 
            return;
          end
          numT = nnz(curidx)-intsize+1;
          int.exp(1,end+1:end+numT) = endx;
          int.flies(1,end+1:end+numT) = flies;
          int.tStart(1,end+1:end+numT) = curt(1:end-intsize+1);
          curwt = (obj.predictdata.loaded(curidx)<0)*negwt +(obj.predictdata.loaded(curidx)>0)*poswt ;
          cumwt = cumsum(curwt);
          sumwt = cumwt(intsize+1:end)-cumwt(1:end-intsize);
          sumwt = [cumwt(intsize) sumwt];
          int.wt(1,end+1:end+numT) = sumwt;
          
        end
      end
      
      cumwt = cumsum(int.wt)/sum(int.wt);
      obj.balancedGTSuggestions = [];
      prevlocs = [];
      for ndx = 1:numint
        while true
            intlocs = rand;
            locsSel = find(cumwt<=intlocs,1,'last');
            
            % Check for overlap
            if any( abs(locsSel-prevlocs) <= intsize), continue, end;
            prevlocs(end+1) = locsSel;
            if isempty(locsSel), locsSel = numel(cumwt); end
            expi = int.exp(locsSel);
            flies = int.flies(locsSel);
            tStart = int.tStart(locsSel);
            obj.balancedGTSuggestions(ndx).start = tStart;
            obj.balancedGTSuggestions(ndx).end = tStart+intsize-1;
            obj.balancedGTSuggestions(ndx).exp = expi;
            obj.balancedGTSuggestions(ndx).flies = flies;
            break;
        end
      end
      
      obj.GTSuggestionMode = 'Balanced';
      
    end
    
    
    function SuggestLoadedGT(obj,expi,filename)
      fid = fopen(filename);
      dat = textscan(fid,'fly:%d,start:%d,end:%d');
      fclose(fid);
      fly = dat{1}; t0s = dat{2}; t1s = dat{3};
      for ndx = 1:obj.nflies_per_exp(expi)
        loc = ismember(fly,ndx);
        if ~any(loc), 
          obj.loadedGTSuggestions{expi}(ndx).start = 1;
          obj.loadedGTSuggestions{expi}(ndx).end = 0;
        else
          obj.loadedGTSuggestions{expi}(ndx).start = t0s(loc);
          obj.loadedGTSuggestions{expi}(ndx).end = t1s(loc);
        end
      end
      obj.GTSuggestionMode = 'Loaded';
    end
    
    function SaveSuggestionGT(obj,expi,filename)
      fid = fopen(filename,'w');
      switch obj.GTSuggestionMode
        
        case 'Random'
          start = obj.randomGTSuggestions{expi}(fly).start;
          last = obj.randomGTSuggestions{expi}(fly).end;
          for fly = 1:obj.nflies_per_exp(expi)
            fprintf(fid,'fly:%d,start:%d,end:%d\n',fly,start,last);
          end
          
        case 'Threshold'
          if ~isempty(obj.predictdata.loaded)
            for fly = 1:obj.nflies_per_exp(expi)
              idxcurr = obj.predictdata.exp(:) == expi & ...
                obj.predictdata.flies(:) == fly & ...
                obj.predictdata.t(:) >=T0 & ...
                obj.predictdata.t(:) <=T1;
              T0 = obj.GetTrxFirstFrame(expi,fly);
              T1 = obj.GetTrxEndFrame(expi,fly);
              suggestedidx = zeros(1,T1);
              suggestedidx( obj.predictdata.t(idxcurr)) = ...
                obj.NormalizeScores(obj.predictdata.loaded(idxcurr)) > ...
                -obj.thresholdGTSuggestions;
              [t0s t1s] = get_interval_ends(suggestedidx);
              for ndx = 1:numel(t0s)
                if t1s(ndx)< T0, continue ; end
                fprintf(fid,'fly:%d,start:%d,end:%d\n',fly,t0s(ndx),t1s(ndx));
              end
            end
          end
          
          
        case 'Balanced'

          for ndx = 1:numel(obj.balancedGTSuggestions)
            if obj.balancedGTSuggestions(ndx).exp ~= expi
              continue;
            end
            start = obj.balancedGTSuggestions(ndx).start;
            last = obj.balancedGTSuggestions(ndx).end;
            fprintf(fid,'fly:%d,start:%d,end:%d\n',...
              obj.balancedGTSuggestions(ndx).flies,start,last);
          end
          
      end
      fclose(fid);
    end
    
    function SuggestThresholdGT(obj,threshold)
      obj.thresholdGTSuggestions = threshold;
      obj.GTSuggestionMode = 'Threshold';
      
    end
    
    function suggestedidx = GetGTSuggestionIdx(obj,expi,flies,T0,T1)
    % Get the indices of gt suggestion  
      if nargin<4,
        T0 = obj.GetTrxFirstFrame(expi,flies);
        T1 = obj.GetTrxEndFrame(expi,flies);
      end
      n = T1-T0+1;
      off = 1 - T0;
      
      if isempty(obj.GTSuggestionMode)
        suggestedidx = false(1,n);
        return;
      end
      
      suggestedidx = false(1,n);
      
      switch obj.GTSuggestionMode,
        case 'Random'
          start = obj.randomGTSuggestions{expi}(flies).start;
          last = obj.randomGTSuggestions{expi}(flies).end;
          range = start+off:last+off;
          selIdx = range(range>0);
          suggestedidx(selIdx) = true;
        
        case 'Loaded'
          if numel(obj.loadedGTSuggestions)<expi || isempty(obj.loadedGTSuggestions{expi}),
            suggestedidx = false(1,n);
            return;
          end
          suggestedidx = false(1,n);
          for ndx = 1:numel(obj.loadedGTSuggestions{expi}(flies).start)
            start = obj.loadedGTSuggestions{expi}(flies).start(ndx);
            last = obj.loadedGTSuggestions{expi}(flies).end(ndx);
            range = start+off:last+off;
            selIdx = range(range>0);
            suggestedidx(selIdx) = true;
          end
          suggestedidx(n+1:end) = [];

        case 'Threshold'
          if ~isempty(obj.predictdata.loaded)
            idxcurr = obj.predictdata.exp(:) == expi & ...
              obj.predictdata.flies(:) == flies & ...
              obj.predictdata.t(:) >=T0 & ...
              obj.predictdata.t(:) <=T1;
            suggestedidx( obj.predictdata.t(idxcurr)+off) = ...
              obj.NormalizeScores(obj.predictdata.loaded(idxcurr)) > ...
              -obj.thresholdGTSuggestions;
          end

          % Should we give suggestions based on scores not calculated offline? 
          if ~any(obj.predictdata.cur_valid)
            idxcurr = obj.FlyNdxPredict(expi,flies) & ...
              obj.predictdata.t(:) >=T0 & ...
              obj.predictdata.t(:) <=T1;
            suggestedidx( obj.predictdata.t(idxcurr)+off) = ...
              obj.NormalizeScores(obj.predictdata.cur(idxcurr)) > ...
              -obj.thresholdGTSuggestions;
          end
          
        case 'Balanced'
          suggestedidx = false(1,n);
          for ndx = 1:numel(obj.balancedGTSuggestions)
            if obj.balancedGTSuggestions(ndx).exp ~= expi ||...
              obj.balancedGTSuggestions(ndx).flies ~= flies
              continue;
            end
            start = obj.balancedGTSuggestions(ndx).start;
            last = obj.balancedGTSuggestions(ndx).end;
            if start>T1 || last <T0, continue ;end
            start = max(start,T0); last = min(last,T1);
            range = start+off:last+off;
            selIdx = range(range>0);
            suggestedidx(selIdx) = true;
          end
      end
      
    end
    
    function crossError = GetGTPerformance(obj)
      % Computes the performance on the GT data.
      obj.StoreLabels();
      crossError.numbers = zeros(4,3);
      crossError.frac = zeros(4,3);
      
      for expi = 1:obj.nexps,
        for i = 1:size(obj.gt_labels(expi).flies,1),
          
          flies = obj.gt_labels(expi).flies(i,:);
          labels_curr = obj.GetLabels(expi,flies);
          ts = [];
          
          for j = 1:numel(labels_curr.t0s),
            ts = [ts,labels_curr.t0s(j):(labels_curr.t1s(j)-1)]; %#ok<AGROW>
          end
          
          % assumes that if have any loaded score for an experiment we
          % have scores for all the flies and for every frame.
          if isempty(obj.predictdata.exp) && ~nnz(obj.predictdata.exp ==expi)
            [success1,msg] = obj.PreLoadWindowData(expi,flies,ts);
            if ~success1,
              warndlg(msg);
              return;
            end
          end
          
        end
      end
      
      obj.PredictLoaded();
      
      gt_scores =[];
      gt_labels = [];
      
      for expi = 1:obj.nexps,
        for i = 1:size(obj.gt_labels(expi).flies,1),
          
          flies = obj.gt_labels(expi).flies(i,:);
          labels_curr = obj.GetLabels(expi,flies);
          
          % Find the importatnt labels
          labels_imp = [];
          for j = 1:numel(labels_curr.imp_t0s),
            t0 = labels_curr.imp_t0s(j);
            t1 = labels_curr.imp_t1s(j);
            labels_imp = [labels_imp t0:t1-1];
          end
          
          for j = 1:numel(labels_curr.t0s),
            t0 = labels_curr.t0s(j);
            t1 = labels_curr.t1s(j);
            
            curLabel = 2*repmat(find(strcmp(labels_curr.names{j},obj.labelnames)),1,t1-t0);
            curLabel(ismember(t0:t1-1,labels_imp)) = curLabel(ismember(t0:t1-1,labels_imp)) -1;
            
            gt_labels = [gt_labels curLabel];
            
            if any(obj.predictdata.loaded_valid),
              idx = obj.FlyNdxPredict(expi,flies)' &...
                obj.predictdata.t(:) >=t0 & obj.predictdata.t(:) <t1;
              ts = obj.predictdata.t(idx);
              scores = obj.predictdata.loaded(idx);
              [check,ndxInLoaded] = ismember(t0:(t1-1),ts);
              if any(check==0), warndlg('Loaded scores are missing scores for some loaded frames'); end
              gt_scores = [gt_scores scores(ndxInLoaded)];
            else
              idx = obj.FlyNdxPredict(expi,flies)' & ...
                obj.predictdata.t(:)>=t0 & obj.predictdata.t(:)<t1;
              ts = obj.predictdata.t(idx);
              scores = obj.predictdata.cur(idx);
              [check,ndxInLoaded] = ismember(t0:(t1-1),ts);
              if any(check==0), warndlg('calculated scores are missing for some labeled frames'); end
              gt_scores = [gt_scores; scores(ndxInLoaded)];
            end
          end
          
        end
      end

      crossError = obj.createConfMat(gt_scores,gt_labels);
      
    
    end

    
% Mode functions
    

    function gtMode =  IsGTMode(obj)
      gtMode = obj.gtMode;
    end
    function advancedMode = IsAdvancedMode(obj)
      advancedMode = obj.advancedMode;
    end
    function modeSet = IsModeSet(obj)
      modeSet = obj.modeSet;
    end
    function SetGTMode(obj,val)
      obj.gtMode = val;
    end
    function SetAdvancedMode(obj,val)
      obj.advancedMode = val;
    end
    function SetMode(obj)
      obj.modeSet = true;
    end
    

% Post Processing functions

    function params = GetPostprocessingParams(obj)
      params = obj.postprocessparams;
    end
    
    function SetPostprocessingParams(obj,params)
      obj.postprocessparams = params;
    end
    
    function blen = GetPostprocessedBoutLengths(obj)
      [success,msg] = obj.ApplyPostprocessing();
      blen = [];
      if ~success;return; end
      
      
      if any(obj.predictdata.cur_valid)
        
        % For predicted scores.
        for endx = 1:obj.nexps
          for flies = 1:obj.nflies_per_exp(endx)
            idx = find(obj.FlyNdxPredict(endx,flies));
            ts = obj.predictdata.t(idx);
            [sortedts, idxorder] = sort(ts);
            gaps = find((sortedts(2:end) - sortedts(1:end-1))>1);
            gaps = [1;gaps;numel(ts)+1];
            for ndx = 1:numel(gaps)-1
              curidx = idx(idxorder(gaps(ndx):gaps(ndx+1)-1));
              posts = obj.predictdata.cur_pp(curidx);
              labeled = bwlabel(posts);
              aa = regionprops(labeled,'Area');
              blen = [blen [aa.Area]];
              
              
            end
          end
        end
        
      else
        
        % For loaded scores.
        for endx = 1:obj.nexps
          for flies = 1:obj.nflies_per_exp(endx)
            curidx = obj.predictdata.exp == endx & obj.predictdata.flies == flies;
            curt = obj.predictdata.t(curidx);
            if any(curt(2:end)-curt(1:end-1) ~= 1)
              msg = 'Scores are not in order';
              success = false;
              return;
            end
            posts = obj.predictdata.loaded_pp(curidx);
            labeled = bwlabel(posts);
            aa = regionprops(labeled,'Area');
            blen = [blen [aa.Area]];
          end
        end
        
      end
    end
    
    function InitPostprocessparams(obj)
      obj.postprocessparams.method = 'Hysteresis';
      obj.postprocessparams.hystopts(1) = struct('name','High Threshold','tag','hthres','value',0);
      obj.postprocessparams.hystopts(2) = struct('name','Low Threshold','tag','lthres','value',0);
      obj.postprocessparams.filtopts(1) = struct('name','Size','tag','size','value',1);
      obj.postprocessparams.blen = 1;
    end
    
    function [success,msg] = ApplyPostprocessing(obj)
    % Applies postprocessing to loaded scores.
    % We do not apply any postprocessing to current scores.
      msg = ''; success = true;
      
      if  any(obj.predictdata.cur_valid)
        
        for endx = 1:obj.nexps
          for flies = 1:obj.nflies_per_exp(endx)
            idx = find(obj.FlyNdxPredict(endx,flies));
            ts = obj.predictdata.t(idx);
            [sortedts, idxorder] = sort(ts);
            gaps = find((sortedts(2:end) - sortedts(1:end-1))>1)+1;
            gaps = [1;gaps;numel(ts)+1];
            for ndx = 1:numel(gaps)-1
              curidx = idx(idxorder(gaps(ndx):gaps(ndx+1)-1));
              curs = obj.predictdata.cur(curidx);
              obj.predictdata.cur_pp(curidx) = obj.Postprocess(curs);
            
            end
          end
        end
        
      else

        for endx = 1:obj.nexps
          for flies = 1:obj.nflies_per_exp(endx)
            curidx = obj.predictdata.exp == endx & obj.predictdata.flies == flies & obj.predictdata.loaded_valid;
            curt = obj.predictdata.t(curidx);
            if any(curt(2:end)-curt(1:end-1) ~= 1)
              msg = 'Scores are not in order';
              success = false;
              return;
            end
            curs = obj.predictdata.loaded(curidx);
            obj.predictdata.loaded_pp(curidx) = obj.Postprocess(curs);
          end %flies
        end
        
      end
      
    end
    
    function posts = Postprocess(obj,curs)
      
      if isempty(obj.postprocessparams); 
        posts = curs;
        return; 
      end;
      
      if strcmpi(obj.postprocessparams.method,'hysteresis')
        posts = obj.ApplyHysteresis(curs,obj.postprocessparams);
      else
        posts = obj.ApplyFiltering(curs,obj.postprocessparams);
      end
      
      posts = obj.RemoveSmallBouts(posts);
      
    end
    
    function posts = RemoveSmallBouts(obj,posts)
        
      if obj.postprocessparams.blen > 1 && numel(posts)>0,
        if numel(posts)<= obj.postprocessparams.blen
          if nnz(posts>0) > numel(posts)/2,
            posts(:) = 1; 
          else
            posts(:) = -1;
          end
            return;
        end
        
        while true,
          tposts = [posts 1-posts(end)];
          ends = find(tposts(1:end-1)~=tposts(2:end));
          ends = [1 ends+1];
          blens = ends(2:end)-ends(1:end-1);
          [minblen,smallbout] = min(blens);
          if minblen>=obj.postprocessparams.blen, break; end
           posts(ends(smallbout):ends(smallbout+1)-1) = ...
              1 - posts(ends(smallbout):ends(smallbout+1)-1);
        end
      end
    end
    
    function posts = ApplyHysteresis(obj,curs,params)
      % Use imfill to find the regions.
      if isempty(curs), posts = curs; return; end
      
      
      % Select pos bouts that have at least one frame about the high
      % threshold.
      hthresh = curs > params.hystopts(1).value*obj.windowdata.scoreNorm;
      lthresh = curs > 0;
      if( nnz(hthresh)>0)
        pos = imfill(~lthresh,find(hthresh(:))) & lthresh;
        computeNeg = true;
      else
        pos = false(size(curs));
        computeNeg = false;
      end
      % Select neg bouts that have at least one frame below the low
      % threshold.
      hthresh = curs < params.hystopts(2).value*obj.windowdata.scoreNorm;
      lthresh = curs < params.hystopts(1).value*obj.windowdata.scoreNorm;
      if nnz(hthresh)>0 && computeNeg,
        neg = imfill(~lthresh,find(hthresh(:))) & lthresh;
      else
        neg = true(size(curs));
      end
      
      posts = pos | ~neg;
      
    end
    
    function posts = ApplyFiltering(obj,curs,params)
      % Use filt to find the regions.
      if isempty(curs), posts = curs; return; end
      filts = conv(curs,ones(1,params.filtopts(1).value),'same');
      posts = filts>0;
    end
    
    function [labels,labeledscores,allScores,scoreNorm] = GetAllLabelsAndScores(obj)
      if isempty(obj.windowdata.exp)
        labels = []; labeledscores = []; 
      else
        curNdx = obj.windowdata.labelidx_cur~=0;
        labeledscores = myBoostClassify(obj.windowdata.X(curNdx,:),obj.classifier);
        origlabels = obj.windowdata.labelidx_cur(curNdx);
        labels = ((origlabels==1)-0.5)*2;
      end
      allScores = obj.predictdata.cur(obj.predictdata.cur_valid>0.5);
      scoreNorm = obj.windowdata.scoreNorm;
    end
    
  end % End methods
    
end % End class




% Old commented functions

%{    
%     function AddWindowDataLabeled(obj,expi,flies)
%       
%       if numel(expi) ~= 1,
%         error('expi must be a scalar');
%       end
%       
%       if nargin < 3,
%         flies = (1:obj.nflies_per_exp(expi))';
%       end
%       
%       for i = 1:size(flies,1),
%         flies_curr = flies(i,:);
%         
%         % labels for this experiment, fly
%         labels_curr = obj.GetLabels(expi,flies_curr);
%         
%         % no labels?
%         if isempty(labels_curr.t0s),
%           continue;
%         end
%         
%         % loop through all behaviors
%         n = max(labels_curr.t1s);
%         idx = zeros(1,n);
%         for j = 1:obj.nbehaviors,
%           for k = find(strcmp(labels_curr.names,obj.labelnames{j})),
%             t0 = labels_curr.t0s(k);
%             t1 = labels_curr.t1s(k);
%             idx(t0+obj.labels_bin_off:t1+obj.labels_bin_off,j) = j;
%           end
%         end
%         m = nnz(idx);
%         if obj.expi == expi && all(flies_curr == obj.flies),
%           obj.windowdata_labeled(:,end+1:end+m) = obj.windowdata_curr(:,idx~=0);
%         else
%           windowfilenames = obj.GetFile('window',expi);
%           % TODO: make this work for multiple flies
%           windowdata_curr = load(windowfilenames{flies_curr(1)});
%           obj.windowdata_labeled(:,end+1:end+m) = windowdata_curr(:,idx~=0);
%         end
%         obj.exp_labeled(end+1:end+m) = expi;
%         obj.flies_labeled(end+1:end+m,:) = repmat(flies_curr,[m,1]);
%         obj.isintrainingset(end+1:end+m) = false;
%         obj.labelidx_labeled(end+1:end+m) = idx;
%         
%       end
%       
%     end
%}

%{
%       function UpdateWindowDataLabeled(obj)
% 
%       % indices into cached data for current experiment and flies
%       idxcurr = obj.exp_labeled' == obj.expi & all(bsxfun(@eq,obj.flies_labeled,obj.flies),2);
%       
%       % frames of current experiment and flies that have old labeled data
%       
%       % indices into cached data that are for this exp, these flies, have
%       % labelidx_cur
%       idxcurr1 = idxcurr & obj.labelidx_old_labeled ~= 0;
%       % which frames for expi, flies that have labelidx_cur ~= 0
%       tsold = obj.ts_labeled(idxcurr1);
%       % indices into labels_bin that have labelidx_cur ~= 0
%       idxold = tsold+obj.labels_bin_off;
%             
%       % keep/add windowdata if labelidx_cur ~= 0 or if labels_bin ~= 0
%       % indices into labels_bin for which new labelidx_new ~= 0
%       cacheidx = any(obj.labels_bin,2);
%       % or labelidx_cur ~= 0
%       cacheidx(idxold) = true;
%       m = nnz(cacheidx);
% 
%       % labelidx_cur for these frames
%       labelidx_cur = zeros(1,m);
%       labelidx_cur(idxold) = obj.labelidx_old_labeled(idxcurr1);
%       
%       % remove all data from the cache for the current exp, flies
%       obj.windowdata_labeled(:,idxcurr) = [];
%       obj.exp_labeled(idxcurr) = [];
%       obj.flies_labeled(idxcurr,:) = [];
%       obj.labelidx_old_labeled(idx) = [];
%       obj.labelidx_new_labeled(idx) = [];
%       obj.ts_labeled(idx) = [];
%       
%       % convert labels_bin to integer
%       n = size(obj.labels_bin,1);
%       labelidx = zeros(1,n);
%       for i = 1:obj.nbehaviors,
%         labelidx(obj.labels_bin(:,i)) = i;
%       end
%       
%       % add this data      
%       obj.windowdata_labeled(:,end+1:end+m) = obj.windowdata_curr(:,cacheidx);
%       obj.exp_labeled(end+1:end+m) = obj.expi;
%       obj.flies_labeled(end+1:end+m,:) = repmat(obj.flies,[m,1]);
%       obj.labelidx_old_labeled(end+1:end+m) = labelidx_cur;
%       obj.labelidx_new_labeled(end+1:end+m) = labelidx;
%       obj.ts_labeled(end+1:end+m) = find(cacheidx) + obj.labels_bin_off;
% 
%     end
%}

%{    
%     function RemoveFromWindowDataLabeled(obj,expi,flies)
% 
%       idx = ismember(obj.exp_labeled,expi);
%       if nargin >= 3,
%         idx = idx & ismember(obj.flies_labeled,flies,'rows');
%       end
%       obj.windowdata_labeled(:,idx) = [];
%       obj.exp_labeled(idx) = [];
%       obj.flies_labeled(idx,:) = [];
%       obj.isintrainingset(idx) = [];
%       obj.labelidx_labeled(idx) = [];
%       
%     end
%}

%{
%     % [success,msg] = LoadTrx(obj,expi)
%     % Load trajectories for input experiment. This should only be called by
%     % PreLoad()!. 
%     function [success,msg] = LoadTrx(obj,expi)
% 
%       success = false;
%       msg = '';
%       
%       if numel(expi) ~= 1,
%         error('expi must be a scalar');
%       end
% 
%       if expi < 1,
%         msg = 'expi not yet set';
%         return;
%       end
%       
%       % load trx
%       try
%         trxfilename = obj.GetFile('trx',expi);
%   
%         hwait = mywaitbar(0,'Loading trx');
%   
%         % TODO: remove this
%         global CACHED_TRX; %#ok<TLEV>
%         global CACHED_TRX_EXPNAME; %#ok<TLEV>
%         if isempty(CACHED_TRX) || isempty(CACHED_TRX_EXPNAME) || ...
%             ~strcmp(obj.expnames{expi},CACHED_TRX_EXPNAME),
%           obj.trx = load_tracks(trxfilename);
%           CACHED_TRX = obj.trx;
%           CACHED_TRX_EXPNAME = obj.expnames{expi};
%         else
%           fprintf('DEBUG: Using CACHED_TRX. REMOVE THIS\n');
%           obj.trx = CACHED_TRX;
%         end
%       catch ME,
%         msg = sprintf('Error loading trx from file %s: %s',trxfilename,getReport(ME));
%         if ishandle(hwait),
%           delete(hwait);
%           drawnow;
%         end
%         return;
%       end
% 
%       if ishandle(hwait),
%         delete(hwait);
%         drawnow;
%       end
%       success = true;
%       
%     end
%}    

%{ 
%     function [success,msg] = LoadWindowData(obj,expi,flies)
% 
%       success = false;
%       msg = '';
%       
%       windowfilenames = obj.GetFile('window',expi);
%       % TODO: make this work for multiple flies
%       try
%         obj.windowdata_curr = load(windowfilenames{flies(1)});
%       catch ME,
%         msg = getReport(ME);
%         return;
%       end
%       
%       success = true;
%       
%     end
%}    

%{
    % trx = GetTrx(obj,expi,flies,ts)
    % Returns the trajectories for the input experiment, flies, and frames.
    % If this is the currently preloaded experiment, then the preloaded
    % trajectories are used. Otherwise, the input experiment is preloaded.
    % If flies is not input, then all flies are returned. If ts is not
    % input, then all frames are returned. 
    function trx = GetTrx(obj,expi,flies,ts)
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,1);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        trx = obj.trx;
        return;
      end
      
      if nargin < 4,
        trx = obj.trx(flies);
        return;
      end
      
      nflies = numel(flies);
      c = cell(1,nflies);
      trx = struct('x',c,'y',c,'a',c,'b',c,'theta',c,'ts',c,'firstframe',c,'endframe',c);
      for i = 1:numel(flies),
        fly = flies(i);
        js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
        trx(i).x = obj.trx(fly).x(js);
        trx(i).y = obj.trx(fly).y(js);
        trx(i).a = obj.trx(fly).a(js);
        trx(i).b = obj.trx(fly).b(js);
        trx(i).theta = obj.trx(fly).theta(js);
        trx(i).ts = js-obj.trx(fly).off;
        trx(i).firstframe = trx(i).ts(1);
        trx(i).endframe = trx(i).ts(end);
      end
    end

    % x = GetTrxX(obj,expi,flies,ts)
    % Returns the x-positions for the input experiment, flies, and frames.
    % This is a cell array with an entry for each fly. If flies is not
    % input, then all flies are returned. If ts is not input, then all
    % frames are returned. 
    function x = GetTrxX(obj,expi,flies,ts)
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        if nargin < 3,
          [success,msg] = obj.PreLoad(expi,1);
        else
          [success,msg] = obj.PreLoad(expi,flies(1));
        end
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        x = {obj.trx.x};
        return;
      end
      
      if nargin < 4,
        x = {obj.trx(flies).x};
        return;
      end
      
      nflies = numel(flies);
      x = cell(1,nflies);
      for i = 1:numel(flies),
        fly = flies(i);
        js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
        x{i} = obj.trx(fly).x(js);
      end
    end
    
    % x = GetTrxX1(obj,expi,fly,ts)
    % Returns the x-positions for the input experiment, SINGLE fly, and
    % frames. If ts is not input, then all frames are returned. 
    function x = GetTrxX1(obj,expi,fly,ts)
      
      if all(expi ~= obj.expi),
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,fly);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end
      
      if nargin < 4,
        x = obj.trx(fly).x;
        return;
      end
      
      x = obj.trx(fly).x(ts + obj.trx(fly).off);

    end

    % y = GetTrxY(obj,expi,flies,ts)
    % Returns the y-positions for the input experiment, flies, and frames.
    % This is a cell array with an entry for each fly. If flies is not
    % input, then all flies are returned. If ts is not input, then all
    % frames are returned. 
    function y = GetTrxY(obj,expi,flies,ts)
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,1);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        y = {obj.trx.y};
        return;
      end
      
      if nargin < 4,
        y = {obj.trx(flies).y};
        return;
      end
      
      nflies = numel(flies);
      y = cell(1,nflies);
      for i = 1:numel(flies),
        fly = flies(i);
        js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
        y{i} = obj.trx(fly).y(js);
      end
    end
    
    % y = GetTrxY1(obj,expi,fly,ts)
    % Returns the y-positions for the input experiment, SINGLE fly, and
    % frames. If ts is not input, then all frames are returned. 
    function y = GetTrxY1(obj,expi,fly,ts)
      
      if all(expi ~= obj.expi),
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,fly);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end
      
      if nargin < 4,
        y = obj.trx(fly).y;
        return;
      end
      
      y = obj.trx(fly).y(ts + obj.trx(fly).off);

    end

    % a = GetTrxA(obj,expi,flies,ts)
    % Returns the quarter major axis lengths for the input experiment,
    % flies, and frames. This is a cell array with an entry for each fly.
    % If flies is not input, then all flies are returned. If ts is not
    % input, then all frames are returned. 
    function a = GetTrxA(obj,expi,flies,ts)
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,1);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        a = {obj.trx.a};
        return;
      end
      
      if nargin < 4,
        a = {obj.trx(flies).a};
        return;
      end
      
      nflies = numel(flies);
      a = cell(1,nflies);
      for i = 1:numel(flies),
        fly = flies(i);
        js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
        a{i} = obj.trx(fly).a(js);
      end
    end
    
    % a = GetTrxA1(obj,expi,fly,ts)
    % Returns the quarter-major-axes for the input experiment, SINGLE fly, and
    % frames. If ts is not input, then all frames are returned. 
    function a = GetTrxA1(obj,expi,fly,ts)
      
      if all(expi ~= obj.expi),
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,fly);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end
      
      if nargin < 4,
        a = obj.trx(fly).a;
        return;
      end
      
      a = obj.trx(fly).a(ts + obj.trx(fly).off);

    end
    
    % b = GetTrxB(obj,expi,flies,ts)
    % Returns the quarter minor axis lengths for the input experiment,
    % flies, and frames. This is a cell array with an entry for each fly.
    % If flies is not input, then all flies are returned. If ts is not
    % input, then all frames are returned. 
    function b = GetTrxB(obj,expi,flies,ts)
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,1);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        b = {obj.trx.b};
        return;
      end
      
      if nargin < 4,
        b = {obj.trx(flies).b};
        return;
      end
      
      nflies = numel(flies);
      b = cell(1,nflies);
      for i = 1:numel(flies),
        fly = flies(i);
        js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
        b{i} = obj.trx(fly).b(js);
      end
    end
    
    % b = GetTrxB1(obj,expi,fly,ts)
    % Returns the quarter-minor-axes for the input experiment, SINGLE fly, and
    % frames. If ts is not input, then all frames are returned. 
    function b = GetTrxB1(obj,expi,fly,ts)
      
      if all(expi ~= obj.expi),
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,fly);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end
      
      if nargin < 4,
        b = obj.trx(fly).b;
        return;
      end
      
      b = obj.trx(fly).b(ts + obj.trx(fly).off);

    end
    
    % theta = GetTrxTheta(obj,expi,flies,ts)
    % Returns the orientations for the input experiment,
    % flies, and frames. This is a cell array with an entry for each fly.
    % If flies is not input, then all flies are returned. If ts is not
    % input, then all frames are returned. 
    function theta = GetTrxTheta(obj,expi,flies,ts)
      
      if numel(expi) ~= 1,
        error('expi must be a scalar');
      end
      
      if expi ~= obj.expi,
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,1);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end

      if nargin < 3,
        theta = {obj.trx.theta};
        return;
      end
      
      if nargin < 4,
        theta = {obj.trx(flies).theta};
        return;
      end
      
      nflies = numel(flies);
      theta = cell(1,nflies);
      for i = 1:numel(flies),
        fly = flies(i);
        js = min(obj.trx(fly).nframes,max(1,ts + obj.trx(fly).off));
        theta{i} = obj.trx(fly).theta(js);
      end
    end

    % theta = GetTrxTheta1(obj,expi,fly,ts)
    % Returns the orientations for the input experiment, SINGLE fly, and
    % frames. If ts is not input, then all frames are returned. 
    function theta = GetTrxTheta1(obj,expi,fly,ts)
      
      if all(expi ~= obj.expi),
        % TODO: generalize to multiple flies
        [success,msg] = obj.PreLoad(expi,fly);
        if ~success,
          error('Error loading trx for experiment %d: %s',expi,msg);
        end
      end
      
      if nargin < 4,
        theta = obj.trx(fly).theta;
        return;
      end
      
      theta = obj.trx(fly).theta(ts + obj.trx(fly).off);

    end
    %}

%{    
%     function [success,msg] = GenerateWindowFeaturesFiles(obj,expi,doforce)
%       
%       success = false;
%       msg = '';
%       
%       if ~exist('doforce','var'),
%         doforce = false;
%       end
% 
%       hwait = mywaitbar(0,'Computing window features...');
% 
%       filenames = obj.GetFile('window',expi);
%       perframedir = obj.GetFile('perframedir',expi);
%       for fly = 1:obj.nflies_per_exp(expi),
%         filename = filenames{fly};
%         if ~doforce && exist(filename,'file'),
%           fprintf('File %s exists, skipping\n',filename);
%           continue;
%         end
%         try
%         X = [];
%         feature_names = {};
%         for j = 1:numel(obj.perframefns),
%           fn = obj.perframefns{j};
%           perframedata = load(fullfile(perframedir,[fn,'.mat']));
%           hwait = mywaitbar((fly-1+(j-1)/numel(obj.perframefns))/obj.nflies_per_exp(expi),hwait,...
%             sprintf('Computing %s window features (%d/%d) for fly %d/%d',fn,j,...
%             numel(obj.perframefns),fly,numel(perframedata.data)));
%           [x_curr,feature_names_curr] = ...
%             ComputeWindowFeatures(perframedata.data{fly},obj.windowfeaturescellparams.(fn){:});
%           nold = size(X,2);
%           nnew = size(x_curr,2);
%           if nold > nnew,
%             x_curr(:,end+1:end+nold-nnew) = nan;
%           elseif nnew > nold,
%             X(:,end+1:end+nnew-nold) = nan;
%           end
%           X = [X;x_curr]; %#ok<AGROW>
%           feature_names = [feature_names,cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false)]; %#ok<AGROW>
%         end
%         hwait = mywaitbar(fly/obj.nflies_per_exp(expi),hwait,...
%           sprintf('Saving window features for fly %d/%d to file...',fly,obj.nflies_per_exp(expi)));
%         save(filename,'X','feature_names');
%         catch ME,
%           msg = getReport(ME);
%           return;
%         end
%       end
%       if exist('hwait','var') && ishandle(hwait),
%         delete(hwait);
%       end
%       
%     end
%}

%{    
%     function [success,msg] = SetWindowFileName(obj,windowfilename)
%       
%       success = false;
%       msg = '';
% 
%       if ischar(windowfilename),
%         obj.windowfilename = windowfilename;
%         if ischar(obj.windowfilename) && strcmp(windowfilename,obj.windowfilename),
%           success = true;
%           return;
%         end
% 
%         % TODO: check window data for existing experiments, remove bad experiments, update
%         % windowdata_labeled, etc. 
%         
%         [success,msg] = obj.UpdateStatusTable('window');
%       end
%       
%     end
%}
    
