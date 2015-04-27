classdef JLabelData_inworks < handle
  
  properties (Access=public)

    % Array of experiments.
    exps = [];
    
    projconf = [];
    cache = [];
    % statistics of labeled data per experiment
    % labelstats(expi).nflies_labeled is the total number of flies labeled,
    % labelstats(expi).nbouts_labeled is the total number of bouts of
    % behaviors labeled, labelstats(expi).
    labelstats = struct('nflies_labeled',{},'nbouts_labeled',{});
    gt_labelstats = struct('nflies_labeled',{},'nbouts_labeled',{});
    
    
    % stuff cached during prediction
    predict_cache = struct;
    
    % whether there is a movie to show
    ismovie = true;
    
    
    % name of classifier file to save/load classifier from
    classifierfilename = '';
    

    % last used path for loading experiment
    defaultpath = '';
    

    % whether all necessary files for all experiments exist
    allfilesexist = true;

    % whether we can generate any missing files
    filesfixable = true;
    
    % functions for writing text to a status bar
    setstatusfn = '';
    clearstatusfn = '';
    
    % Ground truthing or not
    gtMode = false;
    advancedMode = false;
    modeSet = false;
    
    % Ground truthing suggestion
    randomGTSuggestions = {};
    thresholdGTSuggestions = [];
    loadedGTSuggestions = {};
    GTSuggestionMode = '';
    
  end
  
  methods (Access=private)
    
    
  end
    
  methods (Access=public)


% Constructor


    function obj = JLabelData_inworks(varargin)
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
      
    end
    
    
% Configuration settings.

    %%%% FIXED %%%%
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
    
    %%%% FIXED %%%%
    function defaultpath = GetDefaultPath(obj)
      defaultpath = obj.defaultpath;
    end
    
    %%%% FIXED %%%%
    function exps = GetExps(obj)
      exps = obj.exps;
    end
    
    function [success,msg] = SetClassifierFileName(obj,classifierfilename)
    % [success,msg] = SetClassifierFileName(obj,classifierfilename)
    % Sets the name of the classifier file. If the classifier file exists, 
    % it loads the data stored in the file. This involves removing all the
    % experiments and data currently loaded, setting the config file,
    % setting all the file names set in the config file, setting the
    % experiments to be those listed in the classifier file, clearing all
    % the previously computed window data and computing the window data for
    % all the labeled frames. 
      
      success = false;
      msg = '';
      
      obj.classifierfilename = classifierfilename;
      if ~isempty(classifierfilename) && exist(classifierfilename,'file'),
%         try
          obj.SetStatus('Loading classifier from %s',obj.classifierfilename);

          loadeddata = load(obj.classifierfilename); %,obj.classifiervars{:});

          % remove all experiments
          obj.RemoveExpDirs(1:obj.nexps);
          
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
          
          % featureparamsfilename
%           [success,msg] = obj.SetFeatureParamsFileName(loadeddata.featureparamsfilename);
%           if ~success,error(msg);end
          % load actual window features params instead of filename.
          if all( isfield(loadeddata,{'windowfeaturesparams','windowfeaturescellparams',...
              'basicFeatureTable','featureWindowSize'}))
            obj.UpdatePerframeParams(loadeddata.windowfeaturesparams,...
              loadeddata.windowfeaturescellparams,loadeddata.basicFeatureTable,...
              loadeddata.featureWindowSize);
          end
          
      
          % rootoutputdir
%           [success,msg] = obj.SetRootOutputDir(loadeddata.rootoutputdir);
%           if ~success,error(msg); end
           
          % set experiment directories
          [success,msg] = obj.SetExpDirs(loadeddata.expdirs,loadeddata.outexpdirs,...
            loadeddata.nflies_per_exp,loadeddata.sex_per_exp,loadeddata.frac_sex_per_exp,...
            loadeddata.firstframes_per_exp,loadeddata.endframes_per_exp);
          if ~success,error(msg); end
          
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
        
      end

    end
    
    %%%% WHAT? %%%%
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

    
% Experiment handling

    %%%% FIXED %%%%
    function [success,msg] = AddExpDir(obj,expdir)
    % [success,msg] = AddExpDir(obj,expdir)
    % Add a new experiment to the GUI. If this is the first experiment,
    % then it will be preloaded. 

      success = false; msg = '';
      
      if isnumeric(expdir), return; end
      
      if nargin < 2,
        error('Usage: obj.AddExpDirs(expdir)');
      end
      
      newexp = JExperiments();
      [success,msg] = newexp.SetExpir(expdir,obj.Projconf);
      if ~success, return; end;
      
      % preload this experiment if this is the first experiment added
      if numel(obj.exps) == 1,
        [success1,msg1] = Cache.PreLoad(1,newexp,1);
        if ~success1,
          msg = sprintf('Error caching experiment: %s',msg1);
          return;
        end
        
      else
        
        [success1,msg1] = newexp.GetTrxInfo();
        if ~success1,
          msg = sprintf('Error getting basic trx info: %s',msg1);
        end
      end
      
      [success1,msg1] = obj.projconf.UpdateStatusTable('',expH);
      if ~success1,
        msg = msg1;
        return;
      end
      
      obj.exps(end+1) = newexp;
      % save default path
      obj.defaultpath = expdir;
      
      success = true;
      
    end
   
    %%%% FIXED %%%%
    function [success,msg] = RemoveExpDirs(obj,expi)
      % [success,msg] = RemoveExpDirs(obj,expi)
    % Removes experiments in expi from the GUI. If the currently loaded
    % experiment is removed, then a different experiment may be preloaded. 

      success = false;
      msg = '';
      
      if numel(obj.exps) < expi || (expi < 1),
        msg = sprintf('expi = %s must be in the range 1 < expi < nexps = %d',mat2str(expi),obj.nexps);
        return;
      end
      
      delexpH = obj.exps(expi);
      clear('delexpH');
      
      if obj.expi >= expi
        obj.expi = obj.expi -1;
      end
      
      if obj.expi ==0 && numel(obj.exps)>0
        obj.expi =1;
      end
      
      % update current exp, flies
      if obj.expi ~= 0,
        obj.cache.Preload(obj.exps(obj.expi))
      else 
        obj.cache.ClearCachedPerExpData();
      end
      
      success = true;
      
    end

    %%%% FIXED %%%%
    function [success,msg] = PreLoadWindowData(obj)
    % [success,msg] = PreLoadWindowData(obj)
    % Compute and store the window data for all experiments
      
      for ndx = 1:numel(obj.exps)
        expH = obj.exps(ndx);
        [success,msg] = expH.UpdateWindowData();
        if ~success, return; end;
      end
      
    end
    
    function ShowSelectFeatures(obj)
      selHandle = SelectFeatures(obj);
      uiwait(selHandle);
    end

    function UpdateBoostingBins(obj)
      
      oldBinSize = size(obj.windowdata.bins,2);
      newData = size(obj.windowdata.X,1) - size(obj.windowdata.bins,2);
      if newData>0 && ~isempty(obj.windowdata.binVals)
        obj.windowdata.bins(:,end+1:end+newData) = findThresholdBins(obj.windowdata.X(oldBinSize+1:end,:),obj.windowdata.binVals);
      else
        [obj.windowdata.binVals, obj.windowdata.bins] = findThresholds(obj.windowdata.X,obj.classifier_params);
      end
      
    end

% Training and prediction.    
    

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
          
          toPredict = ~obj.windowdata.isvalidprediction;
          obj.SetStatus('Applying boosting classifier to %d windows',sum(toPredict));
          scores = myBoostClassify(obj.windowdata.X(toPredict,:),obj.classifier);
          obj.windowdata.predicted(toPredict) = -sign(scores)*0.5+1.5;
          obj.windowdata.scores(toPredict) = scores;
          obj.windowdata.isvalidprediction(toPredict) = true;
          if ~isempty(obj.classifier_old),
            obj.windowdata.scores_old(toPredict) = ...
              myBoostClassify(obj.windowdata.X(toPredict,:),obj.classifier_old);
          else
            obj.windowdata.scores_old(toPredict) = 0;
          end
          obj.ClearStatus();
          
      end
            
      % transfer to predictidx for current fly
      if ~isempty(obj.expi) && obj.expi > 0 && ~isempty(obj.flies) && all(obj.flies > 0),
        obj.UpdatePredictedIdx();
      end
      
    end
    
    %%%% WHAT? %%%%
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

    %%%% WHAT? %%%%
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
      
      
      if isempty(obj.windowdata.exp),
        return;
      end
      
      % Overwrite by scores from windowdata.
      idxcurr = obj.FlyNdx(obj.expi,obj.flies) & ...
        obj.windowdata.isvalidprediction;
      obj.predictedidx(obj.windowdata.t(idxcurr)-obj.t0_curr+1) = ...
        obj.windowdata.predicted(idxcurr);
      obj.scoresidx(obj.windowdata.t(idxcurr)-obj.t0_curr+1) = ...
        obj.windowdata.scores(idxcurr);      
      obj.scoreTS(obj.windowdata.t(idxcurr)-obj.t0_curr+1) = ...
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
            bouts.ndx(end+1,:) = obj.FlyNdx(expNdx,flyNdx) & ...
              obj.windowdata.t >= curLabels.t0s(boutNum) & ...
              obj.windowdata.t < curLabels.t1s(boutNum);
            bouts.label(end+1) = find(strcmp(obj.labelnames,curLabels.names{boutNum}));
            bouts.timestamp(end+1) = curLabels.timestamp(boutNum);
          end
          
        end
      end
      
    end
    
    function [crossError,tlabels] = CrossValidate(obj)
    % Cross validate on bouts.
    
      [success,msg] = obj.PreLoadLabeledData();
      if ~success, warning(msg);return;end

      islabeled = obj.windowdata.labelidx_cur ~= 0;
      if ~any(islabeled),                        
        crossError.numbers = zeros(4,3);
        crossError.frac = zeros(4,3);
        crossError.oldNumbers = zeros(4,3);
        crossError.oldFrac = zeros(4,3);
        tlabels = {};
        return; 
      end
      
      if ~strcmp(obj.classifiertype,'boosting'); return; end

      obj.SetStatus('Cross validating the classifier for %d examples...',nnz(islabeled));

      obj.UpdateBoostingBins();

      bouts = obj.getLabeledBouts();
      
      [crossScores, tlabels]=...
        crossValidateBout( obj.windowdata.X, ...
        obj.windowdata.labelidx_cur,bouts,obj,...
        obj.windowdata.binVals,...
        obj.windowdata.bins,obj.classifier_params);%,true);

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
      
      waslabeled = false(1,numel(islabeled));
      waslabeled(1:numel(obj.windowdata.labelidx_old)) = obj.windowdata.labelidx_old~=0;
      oldSelect = waslabeled(islabeled);
      oldScores = crossScores(oldSelect);
      oldLabels = 2*obj.windowdata.labelidx_cur(waslabeled) - obj.windowdata.labelidx_imp(waslabeled);
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
          curScoreNdx = find(obj.scoredata.exp == curExp & obj.scoredata.flies==curFly);
          scoresT = obj.scoredata.t(curScoreNdx);
          [curValidScoreNdx,loc] = ismember(scoresT,curT);
          if nnz(curValidScoreNdx)~=numel(curT)
            warndlg('Scores are missing for some labeled data');
            newError = struct();
            return;
          end
          
          orderedLabels = [orderedLabels; curLabels(loc(loc~=0))];
          orderedLabels_imp = [orderedLabels_imp; curLabels_imp(loc(loc~=0))];
          orderedScores = [orderedScores; obj.scoredata.scores(curScoreNdx(curValidScoreNdx~=0))'];
        end
        if setClassifierfilename,
          classifierfilename = obj.scoredata.classifierfilenames{curExp};
          setClassifierfilename = 0;
        elseif strcmp(classifierfilename,'multiple'),
        elseif ~strcmp(classifierfilename,obj.scoredata.classifierfilenames{curExp}),
          classifierfilename = 'multiple';
        end
          
      end
      
      prediction = -sign(orderedScores)/2+1.5;
      
      modLabels = 2*orderedLabels-orderedLabels_imp;
      
      newError = obj.createConfMat(prediction,modLabels);
      newError.classifierfilename = classifierfilename;
      
    end
    
    
    function InitSimilarFrames(obj, HJLabel)
      obj.frameFig = showSimilarFrames;
      showSimilarFrames('SetJLabelData',obj.frameFig,obj,HJLabel);
      showSimilarFrames('CacheTracksLabeled',obj.frameFig);
      showSimilarFrames('add_prep_list', obj.frameFig);
    end
    
    function SimilarFrames(obj,curTime,JLabelHandles)
      if isempty(obj.frameFig), obj.InitSimilarFrames(JLabelHandles), end
      
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
      
      
      if ~isempty(obj.scoredata.exp==expi)
        expid = obj.scoredata.exp==expi;
        expStats.nscoreframes = nnz(expid);
        expStats.nscorepos = nnz(obj.scoredata.scores(expid)>0);
        if ~isempty(obj.scoredata.classifierfilenames) && ...
            numel(obj.scoredata.classifierfilenames)>=expi
          expStats.classifierfilename = obj.scoredata.classifierfilenames{expi};
        else
          expStats.classifierfilename = '';
        end
      else
        expStats.nscoreframes = [];
        expStats.nscorefrac = [];
        expStats.classifierfilename = '';
      end
      
    end

    function flyStats = GetFlyStats(obj,expi,flyNum)
      % Calculates statistics such as number of labeled bouts, predicted bouts
      % and change in scores.
      
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
      
      if ~isempty(obj.scoredata.exp==expi)
        idxcurr = obj.scoredata.exp==expi & obj.scoredata.flies == flyNum;
        flyStats.nscoreframes_loaded = nnz(idxcurr);
        flyStats.nscorepos_loaded = nnz(obj.scoredata.scores(idxcurr)>0);
        flyStats.nscoreneg_loaded = nnz(obj.scoredata.scores(idxcurr)<0);
%         if ~isempty(obj.scoredata.classifierfilenames)
%           flyStats.classifierfilename = obj.scoredata.classifierfilenames{expi};
%         else
%           flyStats.classifierfilename = '';
%         end
      else
        flyStats.nscoreframes_loaded = [];
        flyStats.nscorepos_loaded = [];
        flyStats.nscoreneg_loaded = [];        
%         flyStats.classifierfilename = '';
      end
      
      if ~isempty(obj.windowdata.exp)
        curNdx = obj.FlyNdx(expi,flyNum);
      else
        curNdx = [];
      end
      
      if any(curNdx) && ~isempty(obj.classifier)
        curScores = obj.windowdata.scores(curNdx);
        curLabels = obj.windowdata.labelidx_cur(curNdx);
        
        curPosMistakes = nnz( curScores<0 & curLabels ==1 );
        curNegMistakes = nnz( curScores>0 & curLabels >1 );

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
          flyStats.one2two = nnz(obj.windowdata.scores(curNdx)<0 ...
            & obj.windowdata.scores_old(curNdx)>0);
          flyStats.two2one = nnz(obj.windowdata.scores(curNdx)>0 ...
            & obj.windowdata.scores_old(curNdx)<0);
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
      
    end
    
    function scores = NormalizeScores(obj,scores)
      
      if isempty(obj.windowdata.scoreNorm) || isnan(obj.windowdata.scoreNorm)
        isLabeled = obj.windowdata.labelidx_cur~=0;
        wScores = obj.windowdata.scores(isLabeled);
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
        obj.randomGTSuggestions{endx} = repmat(struct('start',[],'end',[]),1,obj.nflies_per_exp);
        
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
    
    function SuggestLoadedGT(obj,expi,filename)
      fid = fopen(filename);
      dat = textscan(fid,'fly:%d,start:%d,end:%d');
      fclose(fid);
      fly = dat{1}; t0s = dat{2}; t1s = dat{3};
      for ndx = 1:obj.nflies_per_exp(expi)
        [ism, loc] = ismember(ndx,fly);
        if ~ism, 
          obj.loadedGTSuggestions{expi}(ndx).start = 1;
          obj.loadedGTSuggestions{expi}(ndx).end = 0;
        else
          obj.loadedGTSuggestions{expi}(fly(loc)).start = t0s(loc);
          obj.loadedGTSuggestions{expi}(fly(loc)).end = t1s(loc);
        end
      end
      obj.GTSuggestionMode = 'Loaded';
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
      
      if strcmpi(obj.GTSuggestionMode,'Random')
        start = obj.randomGTSuggestions{expi}(flies).start;
        last = obj.randomGTSuggestions{expi}(flies).end;
        range = start+off:last+off;
        selIdx = range(range>0);
        suggestedidx(selIdx) = true;
        
      elseif strcmpi(obj.GTSuggestionMode,'Loaded')
        if numel(obj.loadedGTSuggestions)<expi || isempty(obj.loadedGTSuggestions{expi}),
          suggestedidx = false(1,n);
          return;
        end
        suggestedidx = false(1,n);
        start = obj.loadedGTSuggestions{expi}(flies).start;
        last = obj.loadedGTSuggestions{expi}(flies).end;
        range = start+off:last+off;
        selIdx = range(range>0);
        suggestedidx(selIdx) = true;

      elseif strcmpi(obj.GTSuggestionMode,'Threshold')
        if ~isempty(obj.scoredata.scores)
          idxcurr = obj.scoredata.exp(:) == expi & ...
            obj.scoredata.flies(:) == flies & ...
            obj.scoredata.t(:) >=T0 & ...
            obj.scoredata.t(:) <=T1;
          suggestedidx( obj.scoredata.t(idxcurr)+off) = ...
            obj.NormalizeScores(obj.scoredata.scores(idxcurr)) > ...
            -obj.thresholdGTSuggestions;
        end
        
        % Should we give suggestions based on scores not calculated offline? 
        if ~isempty(obj.windowdata.scores)
          idxcurr = obj.windowdata.exp(:) == expi & ...
            obj.windowdata.flies(:) == flies & ...
            obj.windowdata.t(:) >=T0 & ...
            obj.windowdata.t(:) <=T1;
          suggestedidx( obj.windowdata.t(idxcurr)+off) = ...
            obj.NormalizeScores(obj.windowdata.scores(idxcurr)) > ...
            -obj.thresholdGTSuggestions;
        end
      end
      
    end
    
    function crossError = GetGTPerformance(obj)
      % Computes the performance on the GT data.
      crossError.numbers = zeros(4,3);
      crossError.frac = zeros(4,3);
      
      for expi = 1:obj.nexps,
        for i = 1:size(obj.gt_labels(expi).flies,1),
          
          flies = obj.gt_labels(expi).flies(i,:);
          labels_curr = obj.GetLabels(expi,flies);
          ts = [];
          
          for j = 1:numel(labels_curr.t0s),
            ts = [ts,labels_curr.t0s(j):labels_curr.t1s(j)]; %#ok<AGROW>
          end
          
          % assumes that if have any loaded score for an experiment we
          % have scores for all the flies and for every frame.
          if isempty(obj.scoredata.exp) && ~nnz(obj.scoredata.exp ==expi)
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
            
            if ~isempty(obj.scoredata.exp) && nnz(obj.scoredata.exp ==expi)
              idx = obj.scoredata.exp(:) ==expi & obj.scoredata.flies(:) == flies &...
                obj.scoredata.t(:) >=t0 & obj.scoredata.t(:) <t1;
              ts = obj.scoredata.t(idx);
              scores = obj.scoredata.scores(idx);
              [check,ndxInLoaded] = ismember(t0:(t1-1),ts);
              if any(check==0), warndlg('Loaded scores are missing scores for some loaded frames'); end
              gt_scores = [gt_scores scores(ndxInLoaded)];
            else
              idx = obj.windowdata.exp(:) ==expi & obj.windowdata.flies(:) == flies &...
                obj.windowdata.t(:)>=t0 & obj.windowdata.t(:)<t1;
              ts = obj.windowdata.t(idx);
              scores = obj.windowdata.scores(idx);
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
    
  end % End methods
    
end % End class

