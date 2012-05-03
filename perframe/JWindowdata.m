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
    expH = [];
  end
  
  methods (Access = public)
    
    function obj = JWindowdata()
      
    end
    
    function AddMissingTs(fly)
      %% Stopped here.. offsets.. augh..
      [labelidxStruct,labelidx_off] = expH.GetLabelIdx(fly);
      ts = find(labelidxStruct.vals ~=0);
      missingts = obj.GetMissingTs(target,ts);

    end
    
    function AddX(obj)
      
    end
    
    function idx = FlyNdx(obj,fly)
      idx = obj.flies == fly;
    end
    
    function ts = GetT(obj,fly)
      ts = obj.t(obj.FlyNdx(fly));
    end
    
    function ts = GetMissingTs(obj,fly,ts)
      ts = setdiff(ts,obj.GetT(obj,fly));
    end
    
    function ClearWindowData(obj)
      % Clears window features and predictions for a clean start when selecting
      % features.
      obj.X = [];
      obj.flies=[];
      obj.t=[];
      obj.labelidx_cur=[];
      obj.labelidx_new=[];
      obj.labelidx_imp=[];
      obj.labelidx_old=[];
      obj.featurenames={{}};
      obj.predicted=[];
      obj.predicted_probs=[];
      obj.distNdx=[];
      obj.scores=[];
      obj.scores_old=[];
      obj.scores_validated=[];
      obj.scoreNorm=[];
      obj.binVals=[];
      obj.bins=[];
      
    end
    
    function TrimWindowData(obj)
      % If the size of windowdata is too large, removes windowdata for
      % unlabeled examples.
      sizeLimit = 8e9; % 5GB.
      classSize = 4;
      ratioLimit = 0.2;
      
      numUnlabeled = nnz(obj.labelidx_new==0);
      numLabeled = nnz(obj.labelidx_new);
      
      if numel(obj.X)*classSize < sizeLimit || numUnlabeled/numLabeled<ratioLimit;
        return;
      end
      
      idx2remove = obj.labelidx_new==0 & ...
        ~obj.FlyNdx(obj.expi,obj.flies);
      if ~any(idx2remove); return; end
      obj.X(idx2remove,:) = [];
      obj.exp(idx2remove,:) = [];
      obj.flies(idx2remove,:) = [];
      obj.t(idx2remove,:) = [];
      obj.labelidx_cur(idx2remove,:) = [];
      obj.labelidx_new(idx2remove,:) = [];
      obj.labelidx_imp(idx2remove,:) = [];
      obj.labelidx_old(idx2remove,:) = [];
      obj.predicted(idx2remove,:) = [];
      obj.scores(idx2remove,:) = [];
      obj.scores_old(idx2remove,:) = [];
      obj.scores_validated(idx2remove,:) = [];
      obj.isvalidprediction(idx2remove,:) = [];
      obj.binVals = [];
      obj.bins = [];
      
    end
    
    
    function PredictLoaded(obj)
      % PredictLoaded(obj)
      % Runs the classifier on all preloaded window data.
      
      if isempty(obj.classifier),
        return;
      end
      
      % apply classifier
      switch obj.classifiertype,
        
        case 'boosting',
          
          scores = myBoostClassify(obj.X,obj.classifier);
          obj.predicted = -sign(scores)*0.5+1.5;
          obj.scores = scores;
          if ~isempty(obj.classifier_old),
            obj.scores_old = myBoostClassify(obj.X,obj.classifier_old);
          else
            obj.scores_old = [];
          end
          obj.ClearStatus();
          
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
        isexp = obj.exp == i;
        for j = 1:numel(trainingdata(i).t0s),
          t0 = trainingdata(i).t0s(j);
          t1 = trainingdata(i).t1s(j);
          l = labelidx(j);
          flies = trainingdata(i).flies(j,:);
          isflies = isexp & all(bsxfun(@eq,obj.flies,flies),2);
          ist = isflies & obj.t >= t0 & obj.t < t1;
          if nnz(ist) ~= (t1-t0),
            error('Sanity check: number of training examples does not match windowdata');
          end
          obj.labelidx_cur(ist) = l;
        end
      end
      
    end
    
    function trainingdata = SummarizeTrainingData(obj)
      % trainingdata = SummarizeTrainingData(obj)
      % Summarize labelidx_cur into trainingdata, which is similar to the
      % form of labels.
      
      trainingdata = struct('t0s',{},'t1s',{},'names',{},'flies',{});
      waslabeled = obj.labelidx_cur;
      for expi = 1:obj.nexps,
        trainingdata(expi) = struct('t0s',[],'t1s',[],'names',{{}},'flies',[]);
        isexp = waslabeled & obj.exp == expi;
        if ~any(isexp),
          continue;
        end
        fliess = unique(obj.flies(isexp,:),'rows');
        for fliesi = 1:size(fliess,1),
          flies = fliess(fliesi,:);
          isflies = isexp & all(bsxfun(@eq,obj.flies,flies),2);
          labelidxs = setdiff(unique(obj.labelidx_cur(isflies)),0);
          for labelidxi = 1:numel(labelidxs),
            labelidx = labelidxs(labelidxi);
            islabel = isflies & labelidx == obj.labelidx_cur;
            ts = sort(obj.t(islabel));
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
      % predicted
      
      if obj.expi == 0,
        return;
      end
      
      n = obj.t1_curr - obj.t0_curr + 1;
      obj.predictedidx = zeros(1,n);
      obj.scoresidx = zeros(1,n);
      obj.scoresidx_old = zeros(1,n);
      obj.scoreTS = zeros(1,n);
      
      
      if isempty(obj.exp),
        return;
      end
      
      % Overwrite by scores from 
      idxcurr = obj.FlyNdx(obj.expi,obj.flies) & ...
        obj.isvalidprediction;
      obj.predictedidx(obj.t(idxcurr)-obj.t0_curr+1) = ...
        obj.predicted(idxcurr);
      obj.scoresidx(obj.t(idxcurr)-obj.t0_curr+1) = ...
        obj.scores(idxcurr);
      obj.scoreTS(obj.t(idxcurr)-obj.t0_curr+1) = ...
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
        ~obj.isvalidprediction & ismember(obj.t,ts);
      
      % apply classifier
      switch obj.classifiertype,
        
        case 'ferns',
          obj.SetStatus('Applying fern classifier to %d windows',nnz(idxcurr));
          [obj.predicted(idxcurr),...
            obj.predicted_probs(idxcurr,:)] = ...
            fernsClfApply(obj.X(idxcurr,:),obj.classifier);
          obj.isvalidprediction(idxcurr) = true;
          
          s = exp(obj.predicted_probs);
          s = bsxfun(@rdivide,s,sum(s,2));
          scores = max(s,[],2);
          idxcurr1 = find(idxcurr);
          idx0 = obj.predicted(idxcurr) == 1;
          idx1 = obj.predicted(idxcurr) > 1;
          obj.scores(idxcurr1(idx1)) = -scores(idx1);
          obj.scores(idxcurr1(idx0)) = scores(idx0);
          
          obj.ClearStatus();
        case 'boosting',
          
          obj.SetStatus('Applying boosting classifier to %d frames',nnz(idxcurr));
          scores = myBoostClassify(obj.X(idxcurr,:),obj.classifier);
          obj.predicted(idxcurr) = -sign(scores)*0.5+1.5;
          obj.scores(idxcurr) = scores;
          if ~isempty(obj.classifier_old)
            obj.scores_old(idxcurr) = ...
              myBoostClassify(obj.X(idxcurr,:),obj.classifier_old);
          else
            obj.scores_old(idxcurr) = zeros(size(scores));
          end
          
          obj.isvalidprediction(idxcurr) = true;
          obj.ClearStatus();
          
      end
      
      obj.UpdatePredictedIdx();
      
    end
    
    
  end
  
end