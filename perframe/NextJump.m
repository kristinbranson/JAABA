classdef NextJump < handle
  
  properties (Access=public)
    curType = 'Bouts in current scores';
    allTypes = {'Bouts in current scores',...
      'Bouts in imported scores',...
      'Bouts in postprocessed imported scores',...
      'Errors in current scores',...
      'Errors in validated scores',...
      'Errors in imported scores',...      
      'Postprocessing Changes',...      
      'High Confidence Errors',...
      'Low Confidence',...
      'Thresholds on perframe values',...
      'Ground Truth Suggestions',...
      'Jump To Similar Frames'};
    seek_behaviors_go = [];
    perframefns = {};
    perframeSelFeatures = [];
    perframeSelThresholds = [];
    perframeComparisonType = [];
    hthresh = 0;
    lthresh = 0;
    compareFramesHandle = [];
  end
  
  methods (Access = public, Static = true)
    
  end % End Static methods
  
  methods (Access= public)
    
    function state = GetState(obj)
      state.curType = obj.curType;
      state.perframeSelFeatures = obj.perframeSelFeatures;
      state.perframeSelThresholds = obj.perframeSelThresholds;
      state.perframeComparisonType = obj.perframeComparisonType;
      state.hthresh = obj.hthresh;
      state.lthresh = obj.lthresh;
    end
    
    function SetState(obj,state)
      if isstruct(state) && isfield(state,'curType') && any(strcmp(state.curType,obj.curType)),
        obj.curType = state.curType;
      end
      if isstruct(state) && isfield(state,'perframeSelFeatures') && ...
          all(state.perframeSelFeatures <= numel(obj.perframefns)) ,
        obj.perframeSelFeatures = state.perframeSelFeatures;
        obj.perframeSelThresholds = state.perframeSelThresholds;
        obj.perframeComparisonType = state.perframeComparisonType;
      else
        obj.perframeSelFeatures = 1;
        obj.perframeSelThresholds = 0;
        obj.perframeComparisonType = 1;
      end
      if isstruct(state) && isfield(state,'hthresh'),
        obj.hthresh = state.hthresh;
        obj.lthresh = state.lthresh;
      end
    end
    
    
    function SetPerframefns(obj,perframefns)
      obj.perframefns = perframefns;
    end
    
    function perframefns = GetPerframefns(obj)
       perframefns = obj.perframefns;
    end
    
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
    
    function SetLowThresh(obj,lthresh)
      obj.lthresh = lthresh;
    end
    
    function SetHighThresh(obj,hthresh)
      obj.hthresh = hthresh;
    end
    
    function t = JumpToStart(obj,data,expi,flies,ts,t0,t1)
      
      switch obj.curType
        case obj.allTypes{1}
          t = obj.Automatic_bout_start(data,expi,flies,ts,t0,t1);
        case obj.allTypes{2}
          t = obj.AutomaticLoaded_bout_start(data,expi,flies,ts,t0,t1);          
        case obj.allTypes{3}
          t = obj.Postprocessed_bout_start(data,expi,flies,ts,t0,t1);          
        case obj.allTypes{4}
          t = obj.Error_bout_start(data,expi,flies,ts,t0,t1);
        case obj.allTypes{5}
          t = obj.AutomaticValidated_error_start(data,expi,flies,ts,t0,t1);          
        case obj.allTypes{6}
          t = obj.Loaded_error_start(data,expi,flies,ts,t0,t1);
        case obj.allTypes{7}
          t = obj.Postprocessed_change_start(data,expi,flies,ts,t0,t1);          
        case obj.allTypes{8}
          t = obj.HighconfError_bout_start(data,expi,flies,ts,t0,t1);
        case obj.allTypes{9}
          t = obj.Lowconf_bout_start(data,expi,flies,ts,t0,t1);
        case obj.allTypes{10}
          t = obj.Threshold_bout_start(data,expi,flies,ts,t0,t1);
        case obj.allTypes{11}
          t = obj.GT_Suggestion_start(data,expi,flies,ts,t0,t1);
        case obj.allTypes{12}
          t = obj.SimilarFramesNext(data,expi,flies,ts,t0,t1);
        otherwise
          t = ts;
      end
    end
    
    function t = JumpToEnd(obj,data,expi,flies,ts,t0,t1)
      
      switch obj.curType
        case obj.allTypes{1}
          t = obj.Automatic_bout_end(data,expi,flies,ts,t0,t1);
        case  obj.allTypes{2}
          t = obj.AutomaticLoaded_bout_end(data,expi,flies,ts,t0,t1);          
        case obj.allTypes{3}
          t = obj.Postprocessed_bout_end(data,expi,flies,ts,t0,t1);          
        case  obj.allTypes{4}
          t = obj.Error_bout_end(data,expi,flies,ts,t0,t1);
        case  obj.allTypes{5}
          t = obj.AutomaticValidated_error_end(data,expi,flies,ts,t0,t1);          
        case obj.allTypes{6}
          t = obj.Loaded_error_end(data,expi,flies,ts,t0,t1);
        case obj.allTypes{7}
          t = obj.Postprocessed_change_end(data,expi,flies,ts,t0,t1);          
        case  obj.allTypes{8}
          t = obj.HighconfError_bout_end(data,expi,flies,ts,t0,t1);
        case  obj.allTypes{9}
          t = obj.Lowconf_bout_end(data,expi,flies,ts,t0,t1);
        case  obj.allTypes{10}
          t = obj.Threshold_bout_end(data,expi,flies,ts,t0,t1);
        case obj.allTypes{11}
          t = obj.GT_Suggestion_end(data,expi,flies,ts,t0,t1);
        case obj.allTypes{12}
          t = obj.SimilarFramesPrevious(data,expi,flies,ts,t0,t1);
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
      
      toUse = predictedidx(j:end);
      toUse(lowconfidx(j:end)) = 0;
      k = find(ismember(toUse,obj.seek_behaviors_go),1);
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
      
      toUse = predictedidx(1:j);
      toUse(lowconfidx(1:j)) = 0;
      k = find(ismember(toUse,obj.seek_behaviors_go),1,'last');
      if isempty(k), return; end
      
      t = t0 + k - 1;
    end

    function t = AutomaticLoaded_bout_start(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      [~,predictedidx] = data.GetLoadedScores(expi,flies,t0,t1);
%       predictedidx = double(scores~=0);
%       predictedidx(scores>0) = 1;
%       predictedidx(scores<0) = 2;
% 
%       lowconfidx = false(size(scores));
%       for behaviori = 1:data.nbehaviors
%         idxScores = (predictedidx == behaviori) & ...
%           (abs(scores)<data.GetConfidenceThreshold(behaviori));
%         lowconfidx(idxScores) = true;
%       end
% 
%       j = find( (predictedidx ~= predictedidx(1)) & ~lowconfidx,1);
      j = find( (predictedidx ~= predictedidx(1)) ,1);
      if isempty(j), return; end
      
      toUse = predictedidx(j:end);
%       toUse(lowconfidx(j:end)) = 0;
      k = find(ismember(toUse,obj.seek_behaviors_go),1);
      if isempty(k), return; end
      
      t = ts + j - 1 + k - 1;
    end
    
    function t = AutomaticLoaded_bout_end(obj,data,expi,flies,ts,t0,t1)

      t = [];
      if t0 >= ts, return; end

      t1 = min(max(ts,t0),t1);
      [~,predictedidx] = data.GetLoadedScores(expi,flies,t0,t1);
%       predictedidx = double(scores~=0);
%       predictedidx(scores>0) = 1;
%       predictedidx(scores<0) = 2;
% 
%       lowconfidx = false(size(scores));
%       for behaviori = 1:data.nbehaviors
%         idxScores = (predictedidx == behaviori) & ...
%           (abs(scores)<data.GetConfidenceThreshold(behaviori));
%         lowconfidx(idxScores) = true;
%       end
%       
%       j = find( (predictedidx ~= predictedidx(end)) & ~lowconfidx,1,'last');
      j = find( (predictedidx ~= predictedidx(end)) ,1,'last');
      if isempty(j), return; end
      
      toUse = predictedidx(1:j);
%       toUse(lowconfidx(1:j)) = 0;
      k = find(ismember(toUse,obj.seek_behaviors_go),1,'last');
      if isempty(k), return; end
      
      t = t0 + k - 1;
    end

    function t = Postprocessed_bout_start(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      [~,predictedidx] = data.GetPostprocessedScores(expi,flies,t0,t1);

      j = find( (predictedidx ~= predictedidx(1)),1);
      if isempty(j), return; end
      
      toUse = predictedidx(j:end);
      k = find(ismember(toUse,obj.seek_behaviors_go),1);
      if isempty(k), return; end
      
      t = ts + j - 1 + k - 1;
    end
    
    function t = Postprocessed_bout_end(obj,data,expi,flies,ts,t0,t1)

      t = [];
      if t0 >= ts, return; end

      t1 = min(max(ts,t0),t1);
      [~,predictedidx] = data.GetPostprocessedScores(expi,flies,t0,t1);

      j = find( (predictedidx ~= predictedidx(end)) ,1,'last');
      if isempty(j), return; end
      
      toUse = predictedidx(1:j);
      k = find(ismember(toUse,obj.seek_behaviors_go),1,'last');
      if isempty(k), return; end
      
      t = t0 + k - 1;
    end
    
    function t = Postprocessed_change_start(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      pscores = data.GetPostprocessedScores(expi,flies,t0,t1);
      lscores = data.GetLoadedScores(expi,flies,t0,t1);
      predictedidx  = sign(pscores)~=sign(lscores);

      j = find(predictedidx ==0,1);
      if isempty(j), return; end
      k = find(predictedidx(j:end)==1,1);
      if isempty(k), return; end
      
      t = ts + j - 1 + k -1 ;
    end
    
    function t = Postprocessed_change_end(obj,data,expi,flies,ts,t0,t1)

      t = [];
      if t0 >= ts, return; end

      t1 = min(max(ts,t0),t1);
      pscores = data.GetPostprocessedScores(expi,flies,t0,t1);
      lscores = data.GetLoadedScores(expi,flies,t0,t1);
      predictedidx  = sign(pscores)~=sign(lscores);

      j = find( (predictedidx == 0) ,1,'last');
      if isempty(j), return; end
      
      k = find( (predictedidx(1:j) == 1) ,1,'last');
      if isempty(k), return; end
 
      t = t0 + k - 1;
    end
    
    function t = AutomaticValidated_error_start(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      scores = data.GetValidatedScores(expi,flies,t0,t1);
      predictedidx = double(scores~=0);
      predictedidx(scores>0) = 1;
      predictedidx(scores<0) = 2;
      
      erroridx = labelidx.vals ~=predictedidx;
      erroridx(~labelidx.imp) = 0;
      
      j = find(erroridx ~= erroridx(1),1);
      if isempty(j), return; end
      
      k = find(ismember(labelidx.vals(j:end),obj.seek_behaviors_go)&erroridx(j:end),1);
      if isempty(k),return;end
      
      t = ts + j - 1 + k - 1;
    end
    
    function t = AutomaticValidated_error_end(obj,data,expi,flies,ts,t0,t1)

      t = [];
      if t0 >= ts, return; end

      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      scores = data.GetValidatedScores(expi,flies,t0,t1);
      predictedidx = double(scores~=0);
      predictedidx(scores>0) = 1;
      predictedidx(scores<0) = 2;
      erroridx = labelidx.vals ~=predictedidx;
      erroridx(~labelidx.imp) = 0;
      j = find(erroridx~= erroridx(end),1,'last');
      if isempty(j), return; end
      
      k = find(ismember(labelidx.vals(1:j),obj.seek_behaviors_go)&erroridx(1:j),1,'last');
      if isempty(k), return; end
      t = t0 + k - 1;
    end
    
    function t = Loaded_error_start(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      scores = data.GetLoadedScores(expi,flies,t0,t1);
      predictedidx = double(scores~=0);
      predictedidx(scores>0) = 1;
      predictedidx(scores<0) = 2;
      
      erroridx = labelidx.vals ~=predictedidx;
      
      j = find(erroridx ~= erroridx(1),1);
      if isempty(j), return; end
      
      k = find(ismember(labelidx.vals(j:end),obj.seek_behaviors_go)&erroridx(j:end),1);
      if isempty(k),return;end
      
      t = ts + j - 1 + k - 1;
    end
    
    function t = Loaded_error_end(obj,data,expi,flies,ts,t0,t1)

      t = [];
      if t0 >= ts, return; end

      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      scores = data.GetLoadedScores(expi,flies,t0,t1);
      predictedidx = double(scores~=0);
      predictedidx(scores>0) = 1;
      predictedidx(scores<0) = 2;
      erroridx = labelidx.vals ~=predictedidx;
      j = find(erroridx~= erroridx(end),1,'last');
      if isempty(j), return; end
      
      k = find(ismember(labelidx.vals(1:j),obj.seek_behaviors_go)&erroridx(1:j),1,'last');
      if isempty(k), return; end
      t = t0 + k - 1;
    end


    function t = Manual_bout_start(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      j = find(labelidx.vals ~= labelidx.vals(1),1);
      if isempty(j), return; end
      
      k = find(ismember(labelidx.vals(j:end),obj.seek_behaviors_go),1);
      if isempty(k),return;end
      
      t = ts + j - 1 + k - 1;
    end
    
    function t = Manual_bout_end(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      j = find(labelidx.vals ~= labelidx.vals(end),1,'last');
      if isempty(j), return; end
      k = find(ismember(labelidx.vals(1:j),obj.seek_behaviors_go),1,'last');
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
      erroridx = (labelidx.vals ~=predictedidx) & labelidx.vals;
      
      j = find(erroridx ~= erroridx(1),1);
      if isempty(j), return; end
      
      k = find(ismember(labelidx.vals(j:end),obj.seek_behaviors_go)&erroridx(j:end),1);
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
      erroridx = (labelidx.vals ~=predictedidx) & labelidx.vals;
      j = find(erroridx~= erroridx(end),1,'last');
      if isempty(j), return; end
      
      k = find(ismember(labelidx.vals(1:j),obj.seek_behaviors_go)&erroridx(1:j),1,'last');
      if isempty(k), return; end
      t = t0 + k - 1;
    end
    
    function t = Lowconf_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      scores = data.NormalizeScores(prediction.scoresidx);
      lowconfidx = false(size(scores));
      
      idxScores = scores>obj.lthresh & scores<obj.hthresh;
      lowconfidx(idxScores) = true;

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
      scores = data.NormalizeScores(prediction.scoresidx);
      lowconfidx = false(size(scores));
      idxScores = scores>obj.lthresh & scores<obj.hthresh;
      lowconfidx(idxScores) = true;
      
      lowconfCandidates = lowconfidx(1:end-1)~=lowconfidx(2:end) & ...
                          lowconfidx(1:end-1)==true;
      
      
      j = find(lowconfCandidates,1,'last');
      if isempty(j), return; end
      t = t0 + j - 1;
    end

    function t = HighconfError_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      erroridx = labelidx.vals ~=predictedidx;
      scores = data.NormalizeScores(prediction.scoresidx);
      highconfidx = false(size(scores));
      idxScores = scores<obj.lthresh | scores>obj.hthresh;
      
      highconfidx(idxScores) = true;
      highconfError = highconfidx & erroridx;
      highconfError(labelidx.vals==0) = 0;

      highconfCandidates = highconfError(2:end)~=highconfError(1:end-1) & ...
              highconfError(2:end)==true;
      
      highconfCandidates(labelidx.vals(2:end)==0) = false;
      j = find(highconfCandidates,1);
      if isempty(j), return; end

      t = ts + j;
    end
    
    function t = HighconfError_bout_end(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      erroridx = labelidx.vals ~=predictedidx;
      scores = data.NormalizeScores(prediction.scoresidx);
      highconfidx = false(size(scores));
      
      idxScores = scores<obj.lthresh | scores>obj.hthresh;
      highconfidx(idxScores) = true;

      highconfError = highconfidx & erroridx;
      highconfError(labelidx.vals==0) = 0;
      
      highconfCandidates = highconfError(1:end-1)~=highconfError(2:end) & ...
                          highconfError(1:end-1)==true;
      highconfCandidates(labelidx.vals(1:end-1)==0) = false;
    
      j = find(highconfCandidates,1,'last');
      if isempty(j), return; end
      t = t0 + j - 1;
    end

    function t = Threshold_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      valid = [];
      for ndx = 1:numel(obj.perframeSelFeatures)
        perframedata = data.GetPerFrameData(expi,flies,...
          obj.perframeSelFeatures(ndx),t0,t1);
        if isempty(valid)
          if obj.perframeComparisonType(ndx) == 1 % Less than.
            valid = perframedata<obj.perframeSelThresholds(ndx);
          else
            valid = perframedata>obj.perframeSelThresholds(ndx);
          end
        else
          if obj.perframeComparisonType(ndx) == 1 % Less than.
            valid = valid & perframedata < obj.perframeSelThresholds(ndx);
          else
            valid = valid & perframedata > obj.perframeSelThresholds(ndx);
          end
          
        end
      end
      validCandidates = valid(2:end)~=valid(1:end-1) & ...
                          valid(2:end)==true;
      
      j = find(validCandidates,1);
      if isempty(j), return; end
      t = ts + j;
    end
    
    function t = Threshold_bout_end(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      valid = [];
      for ndx = 1:numel(obj.perframeSelFeatures)
        perframedata = data.GetPerFrameData(expi,flies,...
          obj.perframeSelFeatures(ndx),t0,t1);
        if isempty(valid)
          if obj.perframeComparisonType(ndx) == 1 % Less than.
            valid = perframedata<obj.perframeSelThresholds(ndx);
          else
            valid = perframedata>obj.perframeSelThresholds(ndx);
          end
        else
          if obj.perframeComparisonType(ndx) == 1 % Less than.
            valid = valid & perframedata < obj.perframeSelThresholds(ndx);
          else
            valid = valid & perframedata > obj.perframeSelThresholds(ndx);
          end
          
        end
      end
      
      validCandidates = valid(1:end-1)~=valid(2:end) & ...
                          valid(1:end-1)==true;
      
      
      j = find(validCandidates,1,'last');
      if isempty(j), return; end
      t = t0 + j - 1;
    end
    
    function t = GT_Suggestion_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      suggestidx = data.GetGTSuggestionIdx(expi,flies,t0,t1);
      
      j = find(suggestidx(2:end) & ~suggestidx(1:end-1),1)+1;
      if isempty(j), return; end
      
      t = ts + j -1;
    end
    
    function t = GT_Suggestion_end(obj,data,expi,flies,ts,t0,t1)
      
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      suggestidx = data.GetGTSuggestionIdx(expi,flies,t0,t1);
      j = find(suggestidx(2:end-1) & ~suggestidx(1:end-2),1,'last')+1;
      if isempty(j), return; end
      
      t = t0 + j - 1;
    end
    
    function t = SimilarFramesNext(obj,data,expi,flies,ts,t0,t1)
      t  = [];
      if isempty(obj.compareFramesHandle) || ~ishandle(obj.compareFramesHandle),
        return;
      end
      
      eventdata.Key = 'uparrow';
      CompareFrames('figure1_WindowKeyPressFcn',obj.compareFramesHandle,...
        eventdata,guidata(obj.compareFramesHandle));
      
    end
    
    function t = SimilarFramesPrevious(obj,data,expi,flies,ts,t0,t1)
      t  = [];
      if isempty(obj.compareFramesHandle) || ~ishandle(obj.compareFramesHandle),
        return;
      end
    
      eventdata.Key = 'downarrow';
      CompareFrames('figure1_WindowKeyPressFcn',obj.compareFramesHandle,...
        eventdata,guidata(obj.compareFramesHandle));
      
      
    end
    
    function SetCompareFramesHandle(obj,handles)
      obj.compareFramesHandle = handles;
    end
    
    function ResetCompareFramesHandle(obj)
      obj.compareFramesHandle = [];
    end
    
  end % End methods
end % End classdef
