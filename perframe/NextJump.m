classdef NextJump < handle
  
  properties (Access=public)
    curType = 'Errors';
    allTypes = {'Automatic','Errors','High Confidence Errors',...
        'Low Confidence','Thresholds'};
    seek_behaviors_go = [];
    
  end
  
  methods (Access = public, Static = true)
    
  end % End Static methods
  
  methods (Access= public)
  
    function str = GetCurrentType(obj)
      str = obj.curType;
    end
    
    function str = GetAllTypes(obj)
      str = obj.allTypes;
    end
    
    function SetCurrentType(obj,str)
      if any(strcmp(str,obj.allTypes))
        obj.curType = str;
      else
        warning('Invalid Jump to choice..');
      end
    end
    
    function SetSeekBehaviorsGo(obj,seek_behaviors_go)
      obj.seek_behaviors_go = seek_behaviors_go;
    end
    
    function seek_behaviors_go = GetSeekBehaviorsGo(obj)
      seek_behaviors_go = obj.seek_behaviors_go;
    end
    
    function t = JumpToStart(obj,data,expi,flies,ts,t0,t1)
      
      switch obj.curType
        case 'Automatic'
          t = obj.Automatic_bout_start(data,expi,flies,ts,t0,t1);
        case 'Errors'
          t = obj.Error_bout_start(data,expi,flies,ts,t0,t1);
        case 'Low Confidence'
          t = obj.Lowconf_bout_start(data,expi,flies,ts,t0,t1);
        otherwise
          t = ts;
      end
    end
    
    function t = JumpToEnd(obj,data,expi,flies,ts,t0,t1)
      
      switch obj.curType
        case 'Automatic'
          t = obj.Automatic_bout_end(data,expi,flies,ts,t0,t1);
        case 'Errors'
          t = obj.Error_bout_end(data,expi,flies,ts,t0,t1);
        case 'Low Confidence'
          t = obj.Lowconf_bout_end(data,expi,flies,ts,t0,t1);
        otherwise
          t = ts;
         
      end
    end
    
    function t = Automatic_bout_start(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      
      scores = data.NormalizeScores(prediction.scoresidx);
      lowconfidx = false(size(scores));
      for behaviori = 1:data.nbehaviors
        idxScores = (predictedidx == behaviori) & ...
          (abs(scores)<data.GetConfidenceThreshold(behaviori));
        lowconfidx(idxScores) = true;
      end

      j = find( (predictedidx ~= predictedidx(1)) & ~lowconfidx,1);
      if isempty(j), return; end
      
      k = find(ismember(predictedidx(j:end)&~lowconfidx(j:end),obj.seek_behaviors_go),1);
      if isempty(k), return; end
      
      t = ts + j - 1 + k - 1;
    end
    
    function t = Automatic_bout_end(obj,data,expi,flies,ts,t0,t1)

      t = [];
      if t0 >= ts, return; end

      t1 = min(max(ts,t0),t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      
      scores = data.NormalizeScores(prediction.scoresidx);
      lowconfidx = false(size(scores));
      for behaviori = 1:data.nbehaviors
        idxScores = (predictedidx == behaviori) & ...
          (abs(scores)<data.GetConfidenceThreshold(behaviori));
        lowconfidx(idxScores) = true;
      end
      
      j = find( (predictedidx ~= predictedidx(end)) & ~lowconfidx,1,'last');
      if isempty(j), return; end
      
      k = find(ismember(predictedidx(1:j)&~lowconfidx(1:j),obj.seek_behaviors_go),1,'last');
      if isempty(k), return; end
      
      t = t0 + k - 1;
    end
    
    function t = Manual_bout_start(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      j = find(labelidx ~= labelidx(1),1);
      if isempty(j), return; end
      
      k = find(ismember(labelidx(j:end),obj.seek_behaviors_go),1);
      if isempty(k),return;end
      
      t = ts + j - 1 + k - 1;
    end
    
    function t = Manual_bout_end(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      j = find(labelidx ~= labelidx(end),1,'last');
      if isempty(j), return; end
      k = find(ismember(labelidx(1:j),obj.seek_behaviors_go),1,'last');
      if isempty(k), return; end
      t = t0 + k - 1;
    end
    
    function t = Error_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      erroridx = labelidx ~=predictedidx;
      
      j = find(erroridx ~= erroridx(1),1);
      if isempty(j), return; end
      
      k = find(ismember(labelidx(j:end),obj.seek_behaviors_go)&erroridx(j:end),1);
      if isempty(k),return;end
      
      t = ts + j - 1 + k - 1;
    end
    
    function t = Error_bout_end(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      erroridx = labelidx ~=predictedidx;
      j = find(erroridx~= erroridx(end),1,'last');
      if isempty(j), return; end
      
      k = find(ismember(labelidx(1:j),obj.seek_behaviors_go)&erroridx(1:j),1,'last');
      if isempty(k), return; end
      t = t0 + k - 1;
    end
    
    function t = Lowconf_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      scores = data.NormalizeScores(prediction.scoresidx);
      lowconfidx = false(size(scores));
      
      for behaviori = 1:data.nbehaviors
        idxScores = (predictedidx == behaviori) & ...
          (abs(scores)<data.GetConfidenceThreshold(behaviori));
        lowconfidx(idxScores) = true;
      end

      lowconfCandidates = lowconfidx(2:end)~=lowconfidx(1:end-1) & ...
                          lowconfidx(2:end)==true;
      
      j = find(lowconfCandidates,1);
      if isempty(j), return; end

      t = ts + j;
    end
    
    function t = Lowconf_bout_end(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      scores = data.NormalizeScores(prediction.scoresidx);
      lowconfidx = false(size(scores));
      for behaviori = 1:data.nbehaviors
        idxScores = (predictedidx == behaviori) & ...
          (abs(scores)<data.GetConfidenceThreshold(behaviori));
        lowconfidx(idxScores) = true;
      end
      
      lowconfCandidates = lowconfidx(1:end-1)~=lowconfidx(2:end) & ...
                          lowconfidx(1:end-1)==true;
      
      
      j = find(lowconfCandidates,1,'last');
      if isempty(j), return; end
      t = t0 + j - 1;
    end
    
    
    
  end % End methods
end % End classdef