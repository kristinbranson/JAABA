function [success,msg,predictblocks,windowdata] = PreLoadWindowData(object,perframei,perframedata,missingts,labelidxStruct,t0_labelidx)
% [success,msg] = PreLoadWindowData(obj,expi,flies,ts)
% Compute and store the window data for experiment expi, flies flies,
% and all frames ts. 
% This function finds all frames that currently do not have window data
% cached. In a loop, it finds the first frame that is missing window
% data, and computes window data for all frames in a chunk of size
% 2*obj.windowdatachunk_radius + 1 after this frame using the function
% ComputeWindowDataChunk. Then, it updates the frames that are missing
% window data. It proceeds in this loop until there are no frames
% in the input ts that lack window data. 

  success = false; msg = '';
  predictblocks.t0=[];
  predictblocks.t1=[];
  windowdata.X=[];
  windowdata.t=[];
  windowdata.labelidx_new=[];
  windowdata.labelidx_imp=[];

  % If there are no per-frame features, declare victory.
  if isempty(fieldnames(object.windowfeaturesparams)) ,
    success=true;
    return
  end

  % no frames missing data?
  if isempty(missingts),
    success = true;
    return;
  end

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
      t1 = t0-1;
      success = true;
      return;
    end

    t0 = find(docompute,1,'first') - off;
    t1 = find(docompute,1,'last') - off;
    i0 = t0 - object.GetTrxFirstFrame + 1;
    i1 = t1 - object.GetTrxFirstFrame + 1;

    fn = object.curperframefns{perframei};

    i11 = min(i1,numel(perframedata));
    [x_curr,~] = ...
      ComputeWindowFeatures(perframedata,object.windowfeaturescellparams.(fn){:},'t0',i0,'t1',i11);  %#ok

    if i11 < i1,
      x_curr(:,end+1:end+i1-i11) = nan;
    end

    X = single(x_curr');    
    
    % only store window data that isn't already cached
    tsnew = t0:t1;  % frame indices in the new chunk
    idxnew = (~ismember(tsnew,object.windowdata_t_flyndx)) & ismember(tsnew,missingts);
      % a boolean array the same size as tsnew, each element true iff
      % that element is not in tscurr, and is in missingts
    m = nnz(idxnew);  % the number of frames for which we now have window data, but we didn't before
    if m==0; return; end  % if we didn't make progress, return, signalling failure

    % Add this chunk to predictblocks, which lists all the chunks for
    % which we do prediction, when we do prediction
    [i0s,i1s] = get_interval_ends(idxnew); i1s = i1s-1;
    for j = 1:numel(i0s),
      predictblocks.t0(end+1) = t0+i0s(j)-1;
      predictblocks.t1(end+1) = t0+i1s(j)-1;
    end

    % add to windowdata
    windowdata.X(end+1:end+m,:) = X(idxnew,:);
    windowdata.t(end+1:end+m,1) = tsnew(idxnew);
    tempLabelsNew = labelidxStruct.vals(t0-t0_labelidx+1:t1-t0_labelidx+1);
    windowdata.labelidx_new(end+1:end+m,1) = tempLabelsNew(idxnew);
    tempLabelsImp = labelidxStruct.imp(t0-t0_labelidx+1:t1-t0_labelidx+1);        
    windowdata.labelidx_imp(end+1:end+m,1) = tempLabelsImp(idxnew);        

    % remove from missingts all ts that were computed in this chunk
    missingts(missingts >= t0 & missingts <= t1) = [];

    % stop if we're done
    if isempty(missingts),
      break;
    end

  end

  % store feature_names -- these shouldn't really change
%   windowdata.featurenames = feature_names;

  success = true;

end  % method