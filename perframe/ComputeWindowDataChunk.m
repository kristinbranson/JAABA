function [success,msg,t0,t1,X,feature_names] = ComputeWindowDataChunk(object,t,mode,forceCalc)
% [success,msg,t0,t1,X,feature_names] = ComputeWindowDataChunk(obj,expi,flies,t)
% Computes a chunk of windowdata near frame t for experiment expi and
% flies flies. if mode is 'start', then the chunk will start at t. if
% it is 'center', the chunk will be centered at t. if mode is 'end',
% the chunk will end at t. by default, mode is 'center'. 
% t0 and t1 define the bounds of the chunk of window data computed. X
% is the nframes x nfeatures window data, feature_names is a cell array
% of length nfeatures containing the names of each feature. 
%
% This function first chooses an interval of frames around t, depending 
% on the mode. it then chooses a subinterval of this interval that
% covers all frames in this interval that do not have window data. This
% defines t0 and t1. 
% 
% It then loops through all the per-frame features, and calls
% ComputeWindowFeatures to compute all the window data for that
% per-frame feature. 
%
% To predict over the whole movie we use forceCalc which
% forces the function to recalculate all the features even though they
% were calculated before.

  success = false; msg = '';  %#ok

  if ~exist('mode','var'), mode = 'center'; end
  if ~exist('forceCalc','var'), forceCalc = false; end

  % Check if the features have been configured.
  % I really don't like this.  The JLabelData is pretty close to being a
  % model in the MVC sense.  As such, it shouldn't be creating a view.
  % Not clear to me what the best way to fix this is, though.
  % --ALT, Feb 4, 2013
  if isempty(fieldnames(object.windowfeaturesparams))
    figureJLabel=findall(0,'tag','figure_JLabel');
    figureSelectFeatures=SelectFeatures(figureJLabel);
    uiwait(figureSelectFeatures);
    if isempty(fieldnames(object.windowfeaturesparams))
      error('No features selected!');
    end
  end

  % choose frames to compute:

  % bound at start and end frame of these flies
  T0 = max(object.GetTrxFirstFrame);
  T1 = min(object.GetTrxEndFrame);

  switch lower(mode),
    case 'center',
      % go forward r to find the end of the chunk
      t1 = min(t+object.windowdatachunk_radius,T1);
      % go backward 2*r to find the start of the chunk
      t0 = max(t1-2*object.windowdatachunk_radius,T0);
      % go forward 2*r again to find the end of the chunk
      t1 = min(t0+2*object.windowdatachunk_radius,T1);
    case 'start',
      t0 = max(t,T0);
      t1 = min(t0+2*object.windowdatachunk_radius,T1);
    case 'end',
      t1 = min(t,T1);
      t0 = max(t1-2*object.windowdatachunk_radius,T0);
    otherwise
      error('Unknown mode %s',mode);
  end

  % find a continuous interval that covers all uncomputed ts between t0
  % and t1
  off = 1-t0;
  n = t1-t0+1;
  docompute = true(1,n);
  if object.not_isempty_windowdata_exp && ~forceCalc,
    tscomputed = object.windowdata_t_flyndx;
    tscomputed = tscomputed(tscomputed >= t0 & tscomputed <= t1);
    docompute(tscomputed+off) = false;
  end

  X = single([]);
  feature_names = {};
  if ~any(docompute),
    t1 = t0-1;
    success = true;
    return;
  end

  t0 = find(docompute,1,'first') - off;
  t1 = find(docompute,1,'last') - off;
  i0 = t0 - object.GetTrxFirstFrame + 1;
  i1 = t1 - object.GetTrxFirstFrame + 1;

  %       try

%   curperframefns = obj.curperframefns;
%   allperframefns = obj.allperframefns;
%   perframeInMemory = ~isempty(obj.flies) && obj.IsCurFly(expi,flies);
%   perframedata_all = obj.perframedata;
%   perframefile = obj.GetPerframeFiles(expi);
  x_curr_all = cell(1,numel(object.curperframefns));
  feature_names_all = cell(1,numel(object.curperframefns));
%   windowfeaturescellparams = obj.windowfeaturescellparams;

  % loop through per-frame fields to check that they exist.
%   for j = 1:numel(object.curperframefns),
%     fn = object.curperframefns{j};        
% 
%     % get per-frame data
%     ndx = find(strcmp(fn,object.allperframefns));
%     if isempty(ndx),
%       success = false;
%       msg = sprintf('Internal error: There is at least one per-frame feature in the vocabulary (%s) that is not in the subdialect.',fn);
%       return;
%     end
% 
%     if ~exist(perframefile{ndx},'file'),
%       if ~isdeployed 
%         if isempty(obj.GetGenerateMissingFiles())
%           res = questdlg(sprintf(['Experiment %s is missing some perframe files '...
%             '(%s, possibly more). Generate now?'],obj.expnames{expi},perframefile{ndx}),...
%             'Generate missing files?','Yes','Cancel','Yes');
%           if strcmpi(res,'Yes');
%             obj.SetGenerateMissingFiles();
%           end
%         else 
%           res = fif(obj.GetGenerateMissingFiles(),'Yes','No');
%         end
%       else
%         res = 'Yes';
%       end
% 
%       if strcmpi(res,'Yes'),
%         for ndx = 1:obj.nexps  
%           [success1,msg1] = obj.GenerateMissingFiles(ndx);
%           if ~success1,
%             success = success1; msg = msg1;
%             return;
%           end
%         end
% 
%       else
%         success = false;
%         msg = sprintf('Cannot compute window data for %s ',obj.expnames{expi});
%         return;
%       end
% 
% 
%     end
% 
%   end

%   parfor j = 1:numel(object.curperframefns),
  for j = 1:numel(object.curperframefns),
    fn = object.curperframefns{j};
%        fprintf('Computing window data for per-frame feature %d: %s\n',j,fn);

    % get per-frame data
    ndx = find(strcmp(fn,object.allperframefns));
%     if object.perframeInMemory,
%       perframedata = object.perframedata_all{ndx};  %#ok
%     else
%       perframedata = load(object.perframefile{ndx});  %#ok
%       perframedata = perframedata.data{flies(1)};  %#ok
%     end
    perframedata = object.perframedata{ndx};

    i11 = min(i1,numel(perframedata));
    [x_curr,feature_names_curr] = ...
      ComputeWindowFeatures(perframedata,object.windowfeaturescellparams.(fn){:},'t0',i0,'t1',i11);  %#ok
    if any(imag(x_curr(:)))
      fprintf('Feature values are complex, check input\n');
    end

    if i11 < i1,
      x_curr(:,end+1:end+i1-i11) = nan;
    end

    x_curr_all{j} = single(x_curr);
    feature_names_all{j} = feature_names_curr;

  end

  feature_names=cell(1,numel(object.curperframefns));
  for j = 1:numel(object.curperframefns),
    fn = object.curperframefns{j};
    x_curr = x_curr_all{j};
    feature_names_curr = feature_names_all{j};
    % add the window data for this per-frame feature to X
    nold = size(X,1);
    nnew = size(x_curr,2);
    if nold > nnew,
      warning(['Number of examples for per-frame feature %s does not match number '...
        'of examples for previous features'],fn);
      x_curr(:,end+1:end+nold-nnew) = nan;
    elseif nnew > nold && ~isempty(X),
      warning(['Number of examples for per-frame feature %s does not match number '...
        'of examples for previous features'],fn);
      X(end+1:end+nnew-nold,:) = nan;
    end
    X(:,end+1:end+size(x_curr,1)) = x_curr';
    % add the feature names
    if nargout>5
      feature_names{j} = cellfun(@(s) [{fn},s],feature_names_curr,'UniformOutput',false); %#ok<AGROW>
    end
  end
  feature_names=[feature_names{:}];
  %       catch ME,
  %         msg = getReport(ME);
  %         return;
  %       end
  %X = single(X);
  success = true;

end  % method