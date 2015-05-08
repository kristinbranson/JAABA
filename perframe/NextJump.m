classdef NextJump < handle
  
  properties (Access=public)
    curType = 'Bouts in current scores';
    allTypes = {'Bouts in current scores',...
      'Bouts in imported scores',...
      'Bouts in postprocessed scores',...
      'Errors in current scores',...
      'Errors in validated scores',...
      'Errors in imported scores',...      
      'Postprocessing Changes',...      
      'High Confidence Errors',...
      'Low Confidence',...
      'Thresholds on perframe values',...
      'Ground Truth Suggestions',...
      'Jump To Similar Frames'};
    seek_behaviors_go = []; % vector of label/behavior indices currently included in seek list
    perframefns = {};
    perframeSelFeatures = [];
    perframeSelThresholds = [];
    perframeComparisonType = [];
    hthresh = 0;
    lthresh = 0;
    compareFramesHandle = [];
    
    % AL 20150226: Seems like NextJump can hold the JLabelData it is
    % associated with
  end
  
  methods (Access=public, Static=true)
    
  end % End Static methods
  
  methods (Access=public)
    
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
    
    function SetCompareFramesHandle(obj,handles)
      obj.compareFramesHandle = handles;
    end
    
    function ResetCompareFramesHandle(obj)
      obj.compareFramesHandle = [];
    end    
    
    function [t,flies,expi] = JumpToStart(obj,data,expi,flies,ts,t0,t1)
      
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
          [t,flies,expi] = obj.GT_Suggestion_start(data,expi,flies,ts,t0,t1);
        case obj.allTypes{12}
          t = obj.SimilarFramesNext(data,expi,flies,ts,t0,t1);
        otherwise
          t = ts;
      end
    end
    
    function [t,flies,expi] = JumpToEnd(obj,data,expi,flies,ts,t0,t1)
      
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
          [t,flies,expi] = obj.GT_Suggestion_end(data,expi,flies,ts,t0,t1);
        case obj.allTypes{12}
          t = obj.SimilarFramesPrevious(data,expi,flies,ts,t0,t1);
        otherwise
          t = ts;
         
      end
    end
    
    function [iLblsSeek,iClsSeek] = seekIndices(obj,data)
      iLblsSeek = obj.seek_behaviors_go;
      iClsSeek = data.iLbl2iCls(iLblsSeek);
      iClsSeek = unique(iClsSeek);
    end
    
    function behaviorsRemoved(obj,iLbls)
      % iLbls: label/behavior indices for removed behaviors
      obj.seek_behaviors_go = setdiff(obj.seek_behaviors_go,iLbls);      
    end
    
    function behaviorsAdded(obj,iLbls)
      % iLbls: label/behavior indices for added behaviors
      sbg = union(obj.seek_behaviors_go,iLbls);
      obj.seek_behaviors_go = sbg(:)';      
    end
    
  end
  
  methods
      
    function t = Automatic_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      prediction = data.GetPredictedIdx(expi,flies,t0,t1); 
      scores = data.NormalizeScores(prediction.scoresidx);
      predictedidx = prediction.predictedidx;
      predictedidx = data.PredictedIdxExpandValues(predictedidx);
      
      lowconfidx = NextJump.hlpLowConfidence(data,scores,predictedidx);
      assert(isequal(size(predictedidx),size(lowconfidx)));
      predictedidx(lowconfidx) = 0;

      [iLblsSeek,iClsSeek] = obj.seekIndices(data);
      [~,j] = NextJump.hlpFirstStartNewBout(predictedidx,2,iLblsSeek,iClsSeek);
      if isempty(j), return; end

      t = ts + j - 1;
    end
    
    function t = Automatic_bout_end(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if t0 >= ts, return; end      
      t1 = min(max(ts,t0),t1);
      
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      scores = data.NormalizeScores(prediction.scoresidx);
      predictedidx = prediction.predictedidx;
      predictedidx = data.PredictedIdxExpandValues(predictedidx);
      
      nSamp = size(predictedidx,2);
      if nSamp==1, return; end
      
      lowconfidx = obj.hlpLowConfidence(data,scores,predictedidx);
      assert(isequal(size(predictedidx),size(lowconfidx)));
      predictedidx(lowconfidx) = 0;
      
      [iLblsSeek,iClsSeek] = obj.seekIndices(data);
      [~,j] = NextJump.hlpLastEndNewBout(predictedidx,nSamp-1,iLblsSeek,iClsSeek);
      if isempty(j), return; end      

      t = t0 + j - 1;
    end

    function t = AutomaticLoaded_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      [~,predictedidx] = data.GetLoadedScores(expi,flies,t0,t1);
      predictedidx = data.PredictedIdxExpandValues(predictedidx);

      [iLblsSeek,iClsSeek] = obj.seekIndices(data);
      [~,j] = obj.hlpFirstStartNewBout(predictedidx,2,iLblsSeek,iClsSeek);
      if isempty(j), return; end

      t = ts + j - 1;
    end
    
    function t = AutomaticLoaded_bout_end(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if t0 >= ts, return; end

      t1 = min(max(ts,t0),t1);

      [~,predictedidx] = data.GetLoadedScores(expi,flies,t0,t1);
      predictedidx = data.PredictedIdxExpandValues(predictedidx);
      nSamp = size(predictedidx,2);
      if nSamp==1, return; end

      [iLblsSeek,iClsSeek] = obj.seekIndices(data);
      [~,j] = obj.hlpLastEndNewBout(predictedidx,nSamp-1,iLblsSeek,iClsSeek);
      if isempty(j), return; end

      t = t0 + j - 1;
    end

    function t = Postprocessed_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      [~,predictedidx] = data.GetPostprocessedScores(expi,flies,t0,t1);
      predictedidx = data.PredictedIdxExpandValues(predictedidx);

      [iLblsSeek,iClsSeek] = obj.seekIndices(data);
      [~,j] = obj.hlpFirstStartNewBout(predictedidx,2,iLblsSeek,iClsSeek);
      if isempty(j), return; end
      
      t = ts + j - 1;
    end
    
    function t = Postprocessed_bout_end(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if t0 >= ts, return; end

      t1 = min(max(ts,t0),t1);
      [~,predictedidx] = data.GetPostprocessedScores(expi,flies,t0,t1);
      predictedidx = data.PredictedIdxExpandValues(predictedidx);
      nSamp = size(predictedidx,2);
      if nSamp==1, return; end
      
      [iLblsSeek,iClsSeek] = obj.seekIndices(data);
      [~,j] = obj.hlpLastEndNewBout(predictedidx,nSamp-1,iLblsSeek,iClsSeek);
      if isempty(j), return; end
 
      t = t0 + j - 1;
    end
    
    function t = Postprocessed_change_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      pscores = data.GetPostprocessedScores(expi,flies,t0,t1);
      lscores = data.GetLoadedScores(expi,flies,t0,t1);
      predictedidx  = sign(pscores)~=sign(lscores);

      [~,iClsSeek] = obj.seekIndices(data);
      [~,j] = obj.hlpFirstStartNewBout(predictedidx,2,1,iClsSeek);
      if isempty(j), return; end      
      
      t = ts + j - 1;
    end
    
    function t = Postprocessed_change_end(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);
      
      pscores = data.GetPostprocessedScores(expi,flies,t0,t1);
      lscores = data.GetLoadedScores(expi,flies,t0,t1);
      predictedidx  = sign(pscores)~=sign(lscores);
      nSamp = size(predictedidx,2);
      if nSamp==1, return; end
      
      [~,iClsSeek] = obj.seekIndices(data);
      [~,j] = obj.hlpLastEndNewBout(predictedidx,nSamp-1,1,iClsSeek);
      if isempty(j), return; end   

      t = t0 + j - 1;
    end
    
    function t = AutomaticValidated_error_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      scores = data.GetValidatedScores(expi,flies,t0,t1);
      [~,j] = obj.hlpFirstStartErrBout(labelidx.vals,scores,data,labelidx.imp,false);
      if isempty(j), return; end
      
      t = ts + j - 1;
    end
    
    function t = AutomaticValidated_error_end(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      scores = data.GetValidatedScores(expi,flies,t0,t1);
      [~,j] = obj.hlpLastEndErrBout(labelidx.vals,scores,data,labelidx.imp,false);
      if isempty(j), return; end   
 
      t = t0 + j - 1;
    end
    
    function t = Loaded_error_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      scores = data.GetLoadedScores(expi,flies,t0,t1);
      [~,j] = obj.hlpFirstStartErrBout(labelidx.vals,scores,data,labelidx.imp,false);
      if isempty(j), return; end
      
      t = ts + j - 1;
    end
    
    function t = Loaded_error_end(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      scores = data.GetLoadedScores(expi,flies,t0,t1);
      [~,j] = obj.hlpLastEndErrBout(labelidx.vals,scores,data,labelidx.imp,false);
      if isempty(j), return; end

      t = t0 + j - 1;
    end
        
    function t = Manual_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      [iLblsSeek,iClsSeek] = obj.seekIndices(data);
      [~,j] = obj.hlpFirstStartNewBout(labelidx.vals,2,iLblsSeek,iClsSeek);
      if isempty(j), return; end
      
      t = ts + j - 1;
    end
    
    function t = Manual_bout_end(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      nSamp = size(labelidx.vals,2);
      if nSamp==1, return; end

      [iLblsSeek,iClsSeek] = obj.seekIndices(data);
      [~,j] = obj.hlpLastEndNewBout(labelidx.vals,nSamp-1,iLblsSeek,iClsSeek);
      if isempty(j), return; end

      t = t0 + j - 1;
    end    
    
    function t = Error_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      [~,j] = obj.hlpFirstStartErrBout(labelidx.vals,predictedidx,data,labelidx.vals,true);
      if isempty(j), return; end
                 
      t = ts + j - 1;
    end
    
    function t = Error_bout_end(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);      
      predictedidx = prediction.predictedidx;
      [~,j] = obj.hlpLastEndErrBout(labelidx.vals,predictedidx,data,labelidx.vals,true);
      if isempty(j), return; end                
      
      t = t0 + j - 1;
    end
    
    function t = Lowconf_bout_start(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      scores = data.NormalizeScores(prediction.scoresidx);
      
      idxScores = scores>obj.lthresh & scores<obj.hthresh; % same thresholds applied to all classifiers
      lowconfidx = false(size(scores));
      lowconfidx(idxScores) = true;
      
      [~,iClsSeek] = obj.seekIndices(data);
      [~,j] = NextJump.hlpFirstStartNewBout(lowconfidx,2,1,iClsSeek);
      if isempty(j), return; end

      t = ts + j - 1;
    end
    
    function t = Lowconf_bout_end(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      scores = data.NormalizeScores(prediction.scoresidx);
      nSamp = size(scores,2);
      if nSamp==1, return; end
      
      idxScores = scores>obj.lthresh & scores<obj.hthresh;
      lowconfidx = false(size(scores));
      lowconfidx(idxScores) = true;
      
      [~,iClsSeek] = obj.seekIndices(data);
      [~,j] = NextJump.hlpLastEndNewBout(lowconfidx,nSamp-1,1,iClsSeek);
      
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
      predictedidx = data.PredictedIdxExpandValues(predictedidx);
      erroridx = labelidx.vals~=predictedidx & labelidx.vals~=0;
      
      scores = data.NormalizeScores(prediction.scoresidx);
      idxScores = scores<obj.lthresh | scores>obj.hthresh;
      highconfidx = false(size(scores));
      highconfidx(idxScores) = true;
      highconfError = highconfidx & erroridx;
      
      [~,iClsSeek] = obj.seekIndices(data);
      [~,j] = NextJump.hlpFirstStartNewBout(highconfError,2,1,iClsSeek);
      if isempty(j), return; end
      
      t = ts + j - 1;
    end
    
    function t = HighconfError_bout_end(obj,data,expi,flies,ts,t0,t1)
      t = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      labelidx = data.GetLabelIdx(expi,flies,t0,t1);
      prediction = data.GetPredictedIdx(expi,flies,t0,t1);
      predictedidx = prediction.predictedidx;
      predictedidx = data.PredictedIdxExpandValues(predictedidx);
      erroridx = labelidx.vals~=predictedidx & labelidx.vals~=0;
      
      scores = data.NormalizeScores(prediction.scoresidx);
      idxScores = scores<obj.lthresh | scores>obj.hthresh;
      highconfidx = false(size(scores));      
      highconfidx(idxScores) = true;
      highconfError = highconfidx & erroridx;
      
      nSamp = size(scores,2);
      [~,iClsSeek] = obj.seekIndices(data);
      [~,j] = NextJump.hlpLastEndNewBout(highconfError,nSamp-1,1,iClsSeek);
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
    
    function [t,flies1,expi1] = GT_Suggestion_start(obj,data,expi,flies,ts,t0,t1) %#ok<INUSL>
      t = []; flies1 = []; expi1 = [];
      if ts >= t1, return; end
      t0 = min(max(ts,t0),t1);
      
      assert(data.nclassifiers==1,'Only supported for single-classifier projects.');
      suggestidx = data.GetGTSuggestionIdx(expi,flies,t0,t1);
      
      j = find(suggestidx(2:end) & ~suggestidx(1:end-1),1)+1;
      if ~isempty(j),         
        flies1 = flies;
        expi1 = expi;
        t = ts + j -1;
        return;        
      else        
        for curexpi = expi:data.nexps
          if curexpi == expi,
            flyStart = flies+1;
          else
            flyStart = 1;
          end
          for curfly = flyStart:data.nflies_per_exp(curexpi)
            suggestidx = data.GetGTSuggestionIdx(curexpi,curfly);
            
            j = find(suggestidx(2:end) & ~suggestidx(1:end-1),1)+1;
            if ~isempty(j),
              t = data.GetTrxFirstFrame(curexpi,curfly) + j -1;
              flies1 = curfly;
              expi1 = curexpi;
              return;
            end
          end
        end        
      end      
    end
    
    function [t,flies1,expi1] = GT_Suggestion_end(obj,data,expi,flies,ts,t0,t1) %#ok<INUSL>
      t = [];flies1 = []; expi1 = [];
      if t0 >= ts, return; end
      t1 = min(max(ts,t0),t1);

      assert(data.nclassifiers==1,'Only supported for single-classifier projects.');
      suggestidx = data.GetGTSuggestionIdx(expi,flies,t0,t1);
      j = find(suggestidx(2:end-1) & ~suggestidx(1:end-2),1,'last')+1;

      if ~isempty(j),        
        flies1 = flies;
        expi1 = expi;
        t = t0 + j -1;
        return;        
      else        
        for curexpi = expi:-1:1
          if curexpi == expi,
            flyend = flies-1;
          else
            flyend = data.nflies_per_exp(curexpi);
          end
          for curfly = flyend:-1:1
            suggestidx = data.GetGTSuggestionIdx(curexpi,curfly);
            
            j = find(suggestidx(2:end) & ~suggestidx(1:end-1),1,'last')+1;
            if ~isempty(j),
              t = data.GetTrxFirstFrame(curexpi,curfly) + j -1;
              flies1 = curfly;
              expi1 = curexpi;
              return;
            end
          end
        end        
      end
      if isempty(j), return; end
      
      t = t0 + j - 1;
    end
    
    function t = SimilarFramesNext(obj,data,~,~,~,~,~)
      t  = [];
      if isempty(obj.compareFramesHandle) || ~ishandle(obj.compareFramesHandle),
        return;
      end
      
      assert(data.nclassifiers==1,'Only supported for single-classifier projects.');
      
      eventdata.Key = 'uparrow';
      CompareFrames('figure1_WindowKeyPressFcn',obj.compareFramesHandle,...
        eventdata,guidata(obj.compareFramesHandle));      
    end
    
    function t = SimilarFramesPrevious(obj,data,~,~,~,~,~)
      t  = [];
      if isempty(obj.compareFramesHandle) || ~ishandle(obj.compareFramesHandle),
        return;
      end
      
      assert(data.nclassifiers==1,'Only supported for single-classifier projects.');
    
      eventdata.Key = 'downarrow';
      CompareFrames('figure1_WindowKeyPressFcn',obj.compareFramesHandle,...
        eventdata,guidata(obj.compareFramesHandle));
    end
        
    function [i,j] = hlpFirstStartErrBout(obj,labelidxvals,scores,data,errmask,tfPred12NotScores)
      nSamp = size(scores,2);
      if nSamp==1
        i = [];
        j = [];
        return; 
      end

      if tfPred12NotScores
        pred12 = scores;
      else
        pred12 = double(scores~=0);
        pred12(scores>0) = 1;
        pred12(scores<0) = 2;
      end        
      predictedidx = data.PredictedIdxExpandValues(pred12);
      
      erroridx = labelidxvals~=predictedidx & errmask;
      predictedidx(~erroridx) = 0; % now: predictedIndices-where-there-are-errs
      
      [iLblsSeek,iClsSeek] = obj.seekIndices(data);
      [i,j] = obj.hlpFirstStartNewBout(predictedidx,2,iLblsSeek,iClsSeek);
    end
    
    function [i,j] = hlpLastEndErrBout(obj,labelidxvals,scores,data,errmask,tfPred12NotScores)
      nSamp = size(scores,2);
      if nSamp==1
        i = [];
        j = [];
        return;
      end
      
      if tfPred12NotScores
        pred12 = scores;
      else
        pred12 = double(scores~=0);
        pred12(scores>0) = 1;
        pred12(scores<0) = 2;
      end
      predictedidx = data.PredictedIdxExpandValues(pred12);
      
      erroridx = labelidxvals~=predictedidx & errmask;
      predictedidx(~erroridx) = 0;
      
      [iLblsSeek,iClsSeek] = obj.seekIndices(data);
      [i,j] = obj.hlpLastEndNewBout(predictedidx,nSamp-1,iLblsSeek,iClsSeek);
    end
    
  end
  
  methods (Static)
    
    function lowconfidx = hlpLowConfidence(data,scores,predIdxExpanded)
      % predIdxExpanded: has values in {0,1,... nLbl}

      lowconfidx = false(size(scores));
      nTL = data.ntimelines;
      assert(isequal(nTL,size(scores,1),size(predIdxExpanded,1)));
      for iTL = 1:nTL
        iLbls = data.iCls2iLbl{iTL};
        confThreshs = data.GetConfidenceThreshold(iLbls);
        for i12 = [1 2]
          idxScores = predIdxExpanded(iTL,:)==iLbls(i12) & ...
            abs(scores(iTL,:))<confThreshs(i12);
          lowconfidx(iTL,idxScores) = true;
        end
      end
    end

    function [i,j] = hlpFirstStartNewBout(m,j0,macc,irows)
      % m: abstract label array, integer values.
      % j0: column at which to start (inclusive) search
      % macc: accepted values for m
      % irows: row indices of m to consider
      %
      % This searches m from left-to-right to find the left-most entry
      % (i,j) which satisfies:
      % * i is in irows
      % * m(i,j) is in macc
      % * m(i,j) differs from m(i,j-1)
      %
      % The first entry that is found is returned in [i,j]. If no entry is
      % found, i and j are returned as [].            
      
      [nTL,nSamp] = size(m);
      % Note: j0>size(m,2) is allowed, and will be handed correctly.      
      
      tf = ismember(m,macc);
      if ~exist('irows','var')
        irows = 1:nTL;
      else
        irows = irows(:)';
      end

      for j = j0:nSamp
        for i = irows
          if tf(i,j) && m(i,j)~=m(i,j-1)
            return;
          end
        end
      end
      
      i = [];
      j = [];
    end
    
    function [i,j] = hlpLastEndNewBout(m,jN,macc,irows)
      % See hlpFirstStartNewBout
      
      [nTL,nSamp] = size(m);  %#ok<NASGU>
      % Note: jN<1 is allowed, and will be handed correctly.
      
      tf = ismember(m,macc);
      if ~exist('irows','var')
        irows = 1:nTL;
      else
        irows = irows(:)';
      end
      
      for j = jN:-1:1
        for i = irows
          if tf(i,j) && m(i,j)~=m(i,j+1)
            return;
          end
        end
      end
      
      i = [];
      j = [];
    end
        
  end
  
end % End classdef
