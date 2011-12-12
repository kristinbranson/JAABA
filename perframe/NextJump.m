classdef NextJump < handle
  
  properties (Access=public)
    
    
  end
  
  methods (Access = public, Static = true)
    
    function t = Automatic_bout_start(data,expi,flies,ts,t0,t1,seek_behaviors_go)
      
      t = [];
      
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      j = find(predictedidx ~= predictedidx(1),1);
      if isempty(j), return; end
      
      k = find(ismember(predictedidx(j:end),seek_behaviors_go),1);
      if isempty(k), return; end
      
      t = ts + j - 1 + k - 1;
    end
    
    
    function t = Automatic_bout_end(data,expi,flies,ts,t0,t1,seek_behaviors_go)

      t = [];
      if t0 >= ts, return; end

      t1 = min(max(ts,t0),t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      j = find(predictedidx ~= predictedidx(end),1,'last');
      if isempty(j), return; end
      
      k = find(ismember(predictedidx(1:j),seek_behaviors_go),1,'last');
      if isempty(k), return; end
      
      t = t0 + k - 1;
    end

    
    function t = Manual_bout_start(data,expi,flies,ts,t0,t1,seek_behaviors_go)
      
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      j = find(labelidx ~= labelidx(1),1);
      if isempty(j), return; end
      
      k = find(ismember(labelidx(j:end),seek_behaviors_go),1);
      if isempty(k),return;end
      
      t = ts + j - 1 + k - 1;
    end
   
    
    function t = Manual_bout_end(data,expi,flies,ts,t0,t1,seek_behaviors_go)
      
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      j = find(labelidx ~= labelidx(end),1,'last');
      if isempty(j), return; end
      k = find(ismember(labelidx(1:j),seek_behaviors_go),1,'last');
      if isempty(k), return; end
      t = t0 + k - 1;
    end
    
    function t = Error_bout_start(data,expi,flies,t0,t1,seek_behaviors_go)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      erroridx = labelidx ~=predictedidx;
      
      j = find(erroridx ~= erroridx(1),1);
      if isempty(j), return; end
      
      k = find(ismember(labelidx(j:end),seek_behaviors_go),1);
      if isempty(k),return;end
      
      t = ts + j - 1 + k - 1;
    end
    
    function t = Error_bout_end(data,expi,flies,ts,t0,t1,seek_behaviors_go)
      
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      erroridx = labelidx ~=predictedidx;
      j = find(erroridx~= erroridx(end),1,'last');
      if isempty(j), return; end
      
      k = find(ismember(labelidx(1:j),seek_behaviors_go),1,'last');
      if isempty(k), return; end
      t = t0 + k - 1;
    end
    
    function t = Lowconf_bout_start(data,expi,flies,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      scores = data.NormalizeScores(prediction.scoresidx);
      lowconfidx = false(size(scores));
      
      for behaviori = 1:handles.data.nbehaviors
        idxScores = (predictedidx == behaviori) & ...
          (abs(scores)<GetConfidenceThreshold(behaviori));
        lowconfidx(idxScores) = true;
      end

      lowconfCandidates = lowconfidx(2:end)~=lowconfidx(1:end-1) & ...
                          lowconfidx(2:end)==true;
      
      j = find(lowconfCandidates,1);
      if isempty(j), return; end

      t = ts + j;
    end
    
    function t = Lowconf_bout_end(data,expi,flies,ts,t0,t1)
      
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      scores = data.NormalizeScores(prediction.scoresidx);
      lowconfidx = false(size(scores));
      for behaviori = 1:handles.data.nbehaviors
        idxScores = (predictedidx == behaviori) & ...
          (abs(scores)<GetConfidenceThreshold(behaviori));
        lowconfidx(idxScores) = true;
      end
      
      lowconfCandidates = lowconfidx(1:end-1)~=lowconfidx(2:end) & ...
                          lowconfidx(1:end-1)==true;
      
      
      j = find(lowconfidx,1,'last');
      if isempty(j), return; end
      t = t0 + j - 1;
    end
    

    
  end % End methods
  
end % End classdef