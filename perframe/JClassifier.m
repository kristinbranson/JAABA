classdef JClassifier < handle
  
  properties (Access = public)
    confThresholds = [];
    classifier = [];
    classifier_old = [];
    lastFullClassifierTrainingSize = 0;
    
    % Classifiers Time Stamp
    classifierTS = 0;
    % type of classifier to use
    classifiertype = 'boosting';
    
    % parameters to learning the classifier. struct fields depend on type
    % of classifier.
    % TODO
    classifier_params = struct('iter',100,'iter_updates',10,...
      'numSample',2500,'numBins',30,'CVfolds',7,...
      'baseClassifierTypes',{'Decision Stumps'},'baseClassifierSelected',1);
    
    % constant: stuff stored in classifier mat file
    classifiervars = {'expdirs','outexpdirs','expnames',...
      'nflies_per_exp','sex_per_exp','frac_sex_per_exp',...
      'firstframes_per_exp','endframes_per_exp',...
      'moviefilename','trxfilename','labelfilename','perframedir','clipsdir',...,'featureparamsfilename',...
      'configfilename','rootoutputdir','classifiertype','classifier','trainingdata','classifier_params',...
      'classifierTS','confThresholds','scoreNorm','windowfeaturesparams','windowfeaturescellparams',...
      'basicFeatureTable','featureWindowSize'};
    
    % data for show similar frames.
    frameFig = [];
    distMat = [];
    bagModels = {};
    
    % Retrain properly
    doUpdate = true;
    
    
  end
  
  methods (Access = public)
    
    function obj = JClassifier(projconf)
      nbehaviors = projconf.nbehaviors;
      obj.confThresholds = zeros(1,nbehaviors);
    end
    
    function res = istrained(obj)
      res = isempty(obj.classifier);
    end
    
    function [success,msg] = SetClassifierFileNameWoExp(obj,classifierfilename)
      
      success = false;
      msg = '';
      
      obj.classifierfilename = classifierfilename;
      if ~isempty(classifierfilename) && exist(classifierfilename,'file'),
        %         try
        obj.SetStatus('Loading classifier from %s',obj.classifierfilename);
        
        loadeddata = load(obj.classifierfilename,obj.classifiervars{:});
        
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
          obj.UpdatePerframeParams(loadeddata.windowfeaturesparams,...
            loadeddata.windowfeaturescellparams,loadeddata.basicFeatureTable,...
            loadeddata.featureWindowSize);
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
      end
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
    
    function LoadClassifier()
      
    end
    
    function Train(obj,doFastUpdates)
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
      if ~any(islabeled),
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
            idx_relabel = islabeled & waslabeled & (obj.windowdata.labelidx_new(1:Nprev) ~= obj.windowdata.labelidx_cur(1:Nprev));
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
              [obj.classifier] = fernsClfRemoveTrainingData(obj.windowdata.labelidx_cur(waslabeled), idx_remove(waslabeled), obj.classifier );
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
            [obj.windowdata.binVals, obj.windowdata.bins] = findThresholds(obj.windowdata.X,obj.classifier_params);
            [obj.classifier, ~] =...
              boostingWrapper( obj.windowdata.X(islabeled,:), ...
              obj.windowdata.labelidx_new(islabeled),obj,...
              obj.windowdata.binVals,...
              obj.windowdata.bins(:,islabeled),obj.classifier_params);
            obj.lastFullClassifierTrainingSize = nnz(islabeled);
            
          else
            oldNumPts = nnz(obj.windowdata.labelidx_cur ~= 0 & obj.windowdata.labelidx_imp );
            newNumPts = nnz(obj.windowdata.labelidx_new ~= 0 & obj.windowdata.labelidx_imp );
            newData = newNumPts - oldNumPts;
            obj.SetStatus('Updating boosting classifier with %d examples...',newData);
            
            oldBinSize = size(obj.windowdata.bins,2);
            newData = size(obj.windowdata.X,1) - size(obj.windowdata.bins,2);
            
            if newData>0
              obj.windowdata.bins(:,end+1:end+newData) = findThresholdBins(obj.windowdata.X(oldBinSize+1:end,:),obj.windowdata.binVals);
            end
            
            obj.classifier_old = obj.classifier;
            [obj.classifier, ~] = boostingUpdate(obj.windowdata.X(islabeled,:),...
              obj.windowdata.labelidx_new(islabeled),...
              obj.classifier,obj.windowdata.binVals,...
              obj.windowdata.bins(:,islabeled),obj.classifier_params);
          end
          obj.classifierTS = now();
          obj.windowdata.labelidx_old = obj.windowdata.labelidx_cur;
          obj.windowdata.labelidx_cur = obj.windowdata.labelidx_new;
          
          obj.windowdata.scoreNorm = [];
          % To later find out where each example came from.
          
          obj.windowdata.isvalidprediction = false(numel(islabeled),1);
          obj.PredictLoaded();
      end
      
      obj.ClearStatus();
      
    end
    
    function res = DoFullTraining(obj,doFastUpdates)
      % Check if we should do fast updates or not.
      res = true;
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
    
    function UpdateBoostingBins(obj)
      
      oldBinSize = size(obj.bins,2);
      newData = size(obj.X,1) - size(obj.bins,2);
      if newData>0 && ~isempty(obj.binVals)
        obj.bins(:,end+1:end+newData) = findThresholdBins(obj.X(oldBinSize+1:end,:),obj.binVals);
      else
        [obj.binVals, obj.bins] = findThresholds(obj.X,obj.classifier_params);
      end
      
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
      [success,msg] = obj.PreLoadWindowData(expi,flies,ts);
      if ~success,
          warning(msg);
          return;
      end
      
      
      % indices into windowdata
      idxcurr = obj.FlyNdx(expi,flies) & ...
        ~obj.windowdata.isvalidprediction & ismember(obj.windowdata.t,ts);
      
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
          
          obj.SetStatus('Applying boosting classifier to %d frames',nnz(idxcurr));
          scores = myBoostClassify(obj.windowdata.X(idxcurr,:),obj.classifier);
          obj.windowdata.predicted(idxcurr) = -sign(scores)*0.5+1.5;
          obj.windowdata.scores(idxcurr) = scores;
          if ~isempty(obj.classifier_old)
            obj.windowdata.scores_old(idxcurr) = ...
              myBoostClassify(obj.windowdata.X(idxcurr,:),obj.classifier_old);
          else
            obj.windowdata.scores_old(idxcurr) = zeros(size(scores));
          end
          
          obj.windowdata.isvalidprediction(idxcurr) = true;
          obj.ClearStatus();
          
      end
      
      obj.UpdatePredictedIdx();
      
    end
    
    function SetConfidenceThreshold(obj,thresholds,ndx)
      obj.confThresholds(ndx) = thresholds;
    end
    
    function thresholds = GetConfidenceThreshold(obj,ndx)
      thresholds =obj.confThresholds(ndx) ;
    end
    
    
    function DoBagging(obj)
      [success,msg] = obj.PreLoadLabeledData();
      
      if ~success, warning(msg);return;end
      
      islabeled = obj.windowdata.labelidx_new ~= 0;
      
      if ~any(islabeled),                        return; end
      if ~strcmp(obj.classifiertype,'boosting'); return; end
      if isempty(obj.classifier), obj.Train;             end
      
      obj.SetStatus('Bagging the classifier with %d examples...',nnz(islabeled));
      
      oldBinSize = size(obj.windowdata.bins,2);
      newData = size(obj.windowdata.X,1) - size(obj.windowdata.bins,2);
      if newData>0 && ~isempty(obj.windowdata.binVals)
        obj.windowdata.bins(:,end+1:end+newData) = findThresholdBins(obj.windowdata.X(oldBinSize+1:end,:),obj.windowdata.binVals);
      else
        [obj.windowdata.binVals, obj.windowdata.bins] = findThresholds(obj.windowdata.X,obj.classifier_params);
      end
      
      [obj.bagModels, obj.distMat] =...
        doBagging( obj.windowdata.X(islabeled,:), ...
        obj.windowdata.labelidx_new(islabeled),obj,...
        obj.windowdata.binVals,...
        obj.windowdata.bins(:,islabeled),obj.classifier_params);
      
      obj.windowdata.distNdx.exp = obj.windowdata.exp(islabeled);
      obj.windowdata.distNdx.flies = obj.windowdata.flies(islabeled);
      obj.windowdata.distNdx.t = obj.windowdata.t(islabeled);
      obj.windowdata.distNdx.labels = obj.windowdata.labelidx_new(islabeled);
      
      obj.ClearStatus();
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
    
    
  end
  
end