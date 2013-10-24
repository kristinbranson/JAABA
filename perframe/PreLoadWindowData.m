function [success,msg,predictblocks,windowdata] = PreLoadWindowData(object,missingts,labelidxStruct,t0_labelidx)
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

%       % Check that the given experiment index and target indices are valid
%       % obj.CheckExp(expi); obj.CheckFlies(flies);
%       
%       [labelidxStruct,t0_labelidx] = obj.GetLabelIdx(expi,flies);
% 
%       % which frames don't have window data, which do
%       if isempty(obj.windowdata.exp),
%         missingts = ts;
%         tscurr = [];
%       else      
%         idxcurr = obj.FlyNdx(expi,flies);
%         tscurr = obj.windowdata.t(idxcurr);
%         obj.windowdata.labelidx_new(idxcurr) = labelidxStruct.vals(tscurr-t0_labelidx+1);
%         obj.windowdata.labelidx_imp(idxcurr) = labelidxStruct.imp(tscurr-t0_labelidx+1);
%         missingts = setdiff(ts,tscurr);
%       end
  % tscurr: frame indices that do have window data for the whole track
  %         (not just ts)
  % missingts: frame indices in ts that do not have window data

  % no frames missing data?
  if isempty(missingts),
    success = true;
    return;
  end

  % get labels for current flies -- will be used when filling in
  % windowdata

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

    % update the status
    %obj.SetStatus('Computing windowdata for exp %s, target %d: %d%% done...',...
    %  obj.expnames{expi},flies,round(100*(nts0-numel(missingts))/nts0));

    % compute window data for a chunk starting at t
    if(~exist('feature_names','var'))
      [success1,msg,t0,t1,X,feature_names] = ComputeWindowDataChunk(object,t,'center');
    else
      [success1,msg,t0,t1,X] = ComputeWindowDataChunk(object,t,'center');
    end
    if ~success1, warning(msg); return; end

    % only store window data that isn't already cached
    tsnew = t0:t1;  % frame indices in the new chunk
    idxnew = (~ismember(tsnew,object.windowdata_t_flyndx)) & ismember(tsnew,missingts);
      % a boolean array the same size as tsnew, each element true iff
      % that element is not in tscurr, and is in missingts
    m = nnz(idxnew);  % the number of frames for which we now have window data, but we didn't before
    if m==0; return; end  % if we didn't make progress, return, signalling failure

    % Add this chunk to predictblocks, which lists all the chunks for
    % which we do prediction, when we do prediction
%         obj.predictblocks.expi(end+1) = expi;
%         obj.predictblocks.flies(end+1) = flies;
%         obj.predictblocks.t0(end+1) = t0;
%         obj.predictblocks.t1(end+1) = t1;
    [i0s,i1s] = get_interval_ends(idxnew); i1s = i1s-1;
    for j = 1:numel(i0s),
%           obj.predictblocks.expi(end+1) = expi;
%           obj.predictblocks.flies(end+1) = flies;
%           obj.predictblocks.t0(end+1) = t0+i0s(j)-1;
%           obj.predictblocks.t1(end+1) = t0+i1s(j)-1;
      predictblocks.t0(end+1) = t0+i0s(j)-1;
      predictblocks.t1(end+1) = t0+i1s(j)-1;
    end

    % add to windowdata
%         obj.windowdata.X(end+1:end+m,:) = X(idxnew,:);
%         obj.windowdata.exp(end+1:end+m,1) = expi;
%         obj.windowdata.flies(end+1:end+m,:) = repmat(flies,[m,1]);
%         obj.windowdata.t(end+1:end+m,1) = tsnew(idxnew);
%         obj.windowdata.labelidx_cur(end+1:end+m,1) = 0;
%         tempLabelsNew = labelidxStruct.vals(t0-t0_labelidx+1:t1-t0_labelidx+1);
%         obj.windowdata.labelidx_new(end+1:end+m,1) = tempLabelsNew(idxnew);
%         tempLabelsImp = labelidxStruct.imp(t0-t0_labelidx+1:t1-t0_labelidx+1);        
%         obj.windowdata.labelidx_imp(end+1:end+m,1) = tempLabelsImp(idxnew);        
%         obj.windowdata.labelidx_old(end+1:end+m,1) = 0;
%         obj.windowdata.predicted(end+1:end+m,1) = 0;
%         obj.windowdata.scores(end+1:end+m,1) = 0;
%         obj.windowdata.scores_old(end+1:end+m,1) = 0;   
%         obj.windowdata.scores_validated(end+1:end+m,1) = 0;           
%         obj.windowdata.postprocessed(end+1:end+m,1) = 0;           
%         obj.windowdata.isvalidprediction(end+1:end+m,1) = false;
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
%           obj.ClearStatus();
      break;
    end

  end

  % Clean the window data.
%       obj.CleanWindowData();

  % store feature_names -- these shouldn't really change
  windowdata.featurenames = feature_names;

  success = true;
%       obj.TrimWindowData();

end  % method