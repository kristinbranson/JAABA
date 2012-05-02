classdef JWindowdata < handles
  
  properties ( Access = public)
    X = [];
    flies = [];
    bins = [];
    t = [];
    labelidx_cur=[];
    labelidx_new=[];
    labelidx_old=[];
    labelidx_imp=[];
    featurenames={{}};
    predicted=[];
    distNdx=[]; %todo
    binVals=[]; %todo
    scores=[];
    scoreNorm=[];
    scores_old=[];
    scores_validated=[];

  end
  
  methods (Access = public)
    
    function obj = JWindowdata()
      
    end
    
    function AddX(obj)
    
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
      obj.windowdata.scoreNorm=[];
      obj.windowdata.binVals=[];
      obj.windowdata.bins=[];
      
      obj.UpdatePredictedIdx();

    end
  
    function TrimWindowData(obj)
      % If the size of windowdata is too large, removes windowdata for
      % unlabeled examples.
      sizeLimit = 8e9; % 5GB.
      classSize = 4;
      ratioLimit = 0.2;
      
      numUnlabeled = nnz(obj.windowdata.labelidx_new==0);
      numLabeled = nnz(obj.windowdata.labelidx_new);
      
      if numel(obj.windowdata.X)*classSize < sizeLimit || numUnlabeled/numLabeled<ratioLimit;
        return;
      end
      
      idx2remove = obj.windowdata.labelidx_new==0 & ...
        ~obj.FlyNdx(obj.expi,obj.flies);
      if ~any(idx2remove); return; end
      obj.windowdata.X(idx2remove,:) = [];
      obj.windowdata.exp(idx2remove,:) = [];
      obj.windowdata.flies(idx2remove,:) = [];
      obj.windowdata.t(idx2remove,:) = [];
      obj.windowdata.labelidx_cur(idx2remove,:) = [];
      obj.windowdata.labelidx_new(idx2remove,:) = [];
      obj.windowdata.labelidx_imp(idx2remove,:) = [];
      obj.windowdata.labelidx_old(idx2remove,:) = [];
      obj.windowdata.predicted(idx2remove,:) = [];
      obj.windowdata.scores(idx2remove,:) = [];
      obj.windowdata.scores_old(idx2remove,:) = [];
      obj.windowdata.scores_validated(idx2remove,:) = [];
      obj.windowdata.isvalidprediction(idx2remove,:) = [];
      obj.windowdata.binVals = [];
      obj.windowdata.bins = [];
      
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

    
  end
  
end