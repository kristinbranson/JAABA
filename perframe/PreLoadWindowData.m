function [success,msg,predictblocks,windowdata] = PreLoadWindowData(...
  object,perframefn,perframedata,missingts,labelIdxVals,labelIdxImp,labelIdxT0)
% PreLoadWindowData 
% Compute window data for experiment, flies, classifier as specified in
% dummy-object object to fill in missing times missingts.
%
% object: dummy object created by JLabelData.createPreLoadWindowDataObj
%   using specific expi, flies, iClassifier
% perframefn: perframe feature name
% perframedata: perframe data for perframefn
% missingts: vector of times over which to compute window data, except any 
%   elements that are also elements of object.windowdata_t_flyndx will not be
%   computed (as those times already have window data). missingts does NOT
%   have to be sorted
% labelIdxVals/Imp: vector, labelidx.vals and labelidx.imp for this
%   expi/fly/classifier
% labelIdxT0: T0 for labelIdx
%
% Start with missingts; in a loop, find the first frame that is missing
% window data, and computes window data for all frames in a chunk of size
% 2*obj.windowdatachunk_radius + 1. Iterate until all elements of missingts
% are filled in.

% MERGESTUPDATED

  success = false; 
  msg = '';
  predictblocks.t0 = [];
  predictblocks.t1 = [];
  windowdata.X = [];
  windowdata.t = [];
  windowdata.labelidx_new = [];
  windowdata.labelidx_imp = [];

  % If there are no per-frame features or no frames missing data, declare victory.
  if  isempty(missingts) || object.isempty_fieldnames_windowfeaturesparams
    success = true;
    return;
  end

  while true

    % choose a frame missing window data
    t = median(missingts);
    if ~ismember(t,missingts)
      t = missingts(argmin(abs(t-missingts)));
    end

    % compute window data for a chunk starting at t
    T0 = max(object.GetTrxFirstFrame);
    T1 = min(object.GetTrxEndFrame);

    % go forward r to find the end of the chunk
    t1 = min(t+object.windowdatachunk_radius,T1);
    % go backward 2*r to find the start of the chunk
    t0 = max(t1-2*object.windowdatachunk_radius,T0);
    % go forward 2*r again to find the end of the chunk
    t1 = min(t0+2*object.windowdatachunk_radius,T1);

    % find a continuous interval that covers all uncomputed ts between t0
    % and t1
    off = 1-t0;
    n = t1-t0+1;
    docompute = true(1,n);
    if object.not_isempty_windowdata_exp
      tscomputed = object.windowdata_t_flyndx;
      tscomputed = tscomputed(tscomputed >= t0 & tscomputed <= t1);
      docompute(tscomputed+off) = false;
    end

    if ~any(docompute),
      %t1 = t0-1;
      success = true;
      return;
    end

    t0 = find(docompute,1,'first') - off;
    t1 = find(docompute,1,'last') - off;
    i0 = t0 - object.GetTrxFirstFrame + 1;
    i1 = t1 - object.GetTrxFirstFrame + 1;

%     fn = object.curperframefns{perframei};

    i11 = min(i1,numel(perframedata));
    [x_curr,~] = ...
      ComputeWindowFeatures(perframedata,object.windowfeaturescellparams.(perframefn){:},'t0',i0,'t1',i11);

    if i11 < i1,
      x_curr(:,end+1:end+i1-i11) = nan;
    end

    X = single(x_curr');    
    
    % only store window data that isn't already cached
    tsnew = t0:t1;  % frame indices in the new chunk
    idxnew = (~ismember(tsnew,object.windowdata_t_flyndx)) & ismember(tsnew,missingts);
      % a boolean array the same size as tsnew, each element true iff
      % that element is not in tscurr, and is in missingts
      %
      % AL 20141120: For current single callsite (JLabelData.PreLoadPeri) I 
      % think the first clause is automatic, ie if you are an element of 
      % missingts then you cannot be an element of object.windowdata_t_flyndx.      
    m = nnz(idxnew);  % the number of frames for which we now have window data, but we didn't before
    if m==0
      % if we didn't make progress, return, signaling failure
      return;
    end  

    % Add this chunk to predictblocks, which lists all the chunks for
    % which we do prediction, when we do prediction
    [i0s,i1s] = get_interval_ends(idxnew); 
    i1s = i1s-1;
    for j = 1:numel(i0s)
      predictblocks.t0(end+1) = t0+i0s(j)-1;
      predictblocks.t1(end+1) = t0+i1s(j)-1;
    end

    % add to windowdata
    windowdata.X(end+1:end+m,:) = X(idxnew,:);
    windowdata.t(end+1:end+m,1) = tsnew(idxnew);
    tempLabelsNew = labelIdxVals(t0-labelIdxT0+1:t1-labelIdxT0+1);
    windowdata.labelidx_new(end+1:end+m,1) = tempLabelsNew(idxnew);
    tempLabelsImp = labelIdxImp(t0-labelIdxT0+1:t1-labelIdxT0+1);
    windowdata.labelidx_imp(end+1:end+m,1) = tempLabelsImp(idxnew);

    % remove from missingts all ts that were computed in this chunk
    missingts(missingts >= t0 & missingts <= t1) = [];

    % stop if we're done
    if isempty(missingts),
      break;
    end

  end

  % store feature_names -- these shouldn't really change
  % windowdata.featurenames = feature_names;

  success = true;

end