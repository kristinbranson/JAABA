% [y,feature_names,cache] = ComputeHistWindowFeatures(x,...)
% 
% Computes the histogram window features y for the input one-dimensional
% per-frame time series data x. y(i,t) will correspond to histogram
% window feature i and frame t, and is the fraction of (the transformation
% of) the per-frame data within the window defined by i and t that falls in
% the range specified by the histogram bin.
% 
% For window feature i, let b_i be the bin index, r_i be the window radius,
% and off_i be the window offset. If the transformation type for feature i
% is 'none', the window feature y(i,t) is the fraction of per-frame data
% within the window from t-r_i+off_i through t+r_i+off_i that falls within
% the ranges for bin b_i. If the transformation type for feature i is
% 'abs', then we are histogramming the absolute value of the per-frame
% data. If the transformation type for feature i is 'flip', then we are
% histogramming the sign of x(t) times the per-frame data in the window.
%
% Input:
%
% x: 1 x NFRAMES array of per-frame data. 
% 
% Output: 
%
% y: NWINDOWFEATURES x NFRAMES matrix of window data, where y(i,t)
% corresponds to window feature i and frame t. i indexes the radius and
% offset of the window as well as the transformation type. 
% feature_names: 1 x NWINDOWFEATURES cell in which feature_names{i}
% describes the ith window feature computed. feature_names{i} is itself a
% cell that can be interpreted as pairs of a string description followed by
% a value, e.g. 
% {'stat','hist','trans','abs','radius',1,'offset',1,'bin',3,'lim',[.3,.4]}
% cache: if DOCACHE is set to true, then the non-offset window mins will
% be cached to potentially be used in future window feature computations. 
%
% Optional inputs:
%
% 'hist_edges': The 1 x NBINS+1 edges to use when histogramming. 
%
% DEFAULT WINDOW LOCATIONS:
% These window locations are used if window locations are not specified on
% a per-feature type basis. Default default values set by
% SetDefaultWindowParameters.
% Inputs interpreted by SetWindowParameters. 
% The same window parameter interpretation is done in all
% Compute*WindowFeatures functions. 
% 'windows': window offsets and radii to try ([radius1,offset1];...;[radiusn,offsetn])
% if empty, then window_radii and window_offsets are used to specify
% windows. if the window is [r,off], then the window feature at time t
% will be computed from the window at t-r+off to t+r+off. default value: [].
% 'window_radii': if windows is empty, then the cross-product of window_radii
% and window_offsets is used to set windows. if empty, then
% min_window_radius, max_window_radius, and nwindow_radii are used to set
% window_radii. default value: [].
% 'window_offsets': if windows is empty, then the cross-product of  window_radii
% and window_offsets is used to set windows. window_offsets are relative to
% radius, so the window corresponding to radius r_i and offset off_j_rel
% will be [r_i,off_i=r_i*off_j_rel]. default value: [-1,0,1].
% 'min_window_radius': if windows and window_radii are both empty, then
% window_radii is set to nwindow_radii evenly spaced radii between
% min_window_radii and max_window_radii. default value: 0
% 'max_window_radius': see 'min_window_radius'. default value: 20.
% 'nwindow_radii: see 'min_window_radius'. default value: 5. 
%
% 'trans_types': default types of transformations to apply for each feature
% type, if not otherwise specified. Options include 'abs','flip', and
% 'none'. 'abs' corresponds to the absolute value, flip corresponds to
% flipping the sign of the feature if the per-frame feature at the frame t
% is negative, and 'none' corresponds to no transformation. 
%
% 'sanitycheck': whether to compute all the features a second time in the
% obvious way to make sure that the optimized computations are correct.
% default value: false. 
%
% 'docache': whether to cache computations that might be useful, e.g. the
% min computations are useful when computing the diff_neighbor_min.
% default value: true.
%
% 'cache': input cached computations that might be useful in this
% computation. 

function [y,feature_names,cache] = ComputeHistWindowFeatures(x,varargin)

x = x(:)';
N = numel(x);
y = nan(0,N);
feature_names = {};

%% default parameters

[windows,window_offsets,...
  window_radii,min_window_radius,...
  max_window_radius,nwindow_radii] = ...
  SetDefaultWindowParameters();

% use all transformation types by default
trans_types = 'all';

% for debugging purposes
SANITY_CHECK = true;

% whether to cache results
DOCACHE = true;

% initialize empty cache
cache = InitializeCache();

% initialize feature_types already computed to empty
feature_types = {};

% edges of histogram bins
hist_edges = [];

relativeParams = [];
%% parse parameters

[...
  windows,...
  window_radii,window_offsets,...
  min_window_radius,max_window_radius,nwindow_radii,...
  feature_types,...
  trans_types,...
  SANITY_CHECK,...
  DOCACHE,...
  cache,...
  hist_edges,...
  relativeParams...
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'feature_types',feature_types,...
  'trans_types',trans_types,...
  'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE,...
  'cache',cache,...
  'hist_edges',hist_edges,...
  'relativeParams',relativeParams); %#ok<ASGLU>

%% whether we've specified to use all trans types by default
if ischar(trans_types) && strcmpi(trans_types,'all'),
  trans_types = {'none','abs','flip','relative'};
end

%% select default windows from various ways of specifying windows

[windows,window_radii,windowi2radiusi,nradii] = ...
  SetWindowParameters(...
  windows,window_offsets,...
  window_radii,min_window_radius,...
  max_window_radius,nwindow_radii);

%% compute per-frame transformations 
[x_trans,IDX,ntrans] = ComputePerFrameTrans(x,trans_types);

if ismember('relative',trans_types)
  if DOCACHE && ~isempty(cache.relX)
    modX = cache.relX;
  else
    modX = convertToRelative(x,relativeParams);
    cache.relX = modX;
  end
  x_trans(end+1,:) = modX;
  IDX.rel = size(x_trans,1);
  ntrans = size(x_trans,1);
else
  IDX.rel = 0;
end

%% Bin Edges for relative histogram

rel_hist_edges = linspace(0,100,length(hist_edges));

%% main computation

if isempty(hist_edges),
  warning('hist_edges is empty, but hist features are supposed to be computed. Skipping these features.');
  return;
end

% find bin for each value -- notice the transpose
nbins = numel(hist_edges)-1;
bin = [];
[~,bin(:,IDX.orig)] = histc(x_trans(IDX.orig,:)',hist_edges);

% * No abs and flip transforms for hist. * 
%
% if ismember('none',trans_types)
%   [~,bin(:,IDX.orig)] = histc(x_trans(IDX.orig,:)',hist_edges);
% end
% if IDX.abs>0
%   [~,bin(:,IDX.abs)] = histc(x_trans(IDX.abs,:)',hist_edges);
% end
% if IDX.flip>0
%     [~,bin(:,IDX.flip)] = histc(x_trans(IDX.flip,:)',hist_edges);
% end
% if IDX.rel>0
%   [~,bin(:,IDX.rel)] = histc(x_trans(IDX.rel,:)',rel_hist_edges);
% end

bin(bin > nbins) = nbins;

for radiusi = 1:nradii,
  r = window_radii(radiusi);
  w = 2*r+1;
  
  res = zeros(nbins,N+2*r,ntrans);
  for j = 1:ntrans,
    % initialize first value:
    % t contains [max(1,t-r),min(N,t+r)]
    % t = 1 will contain 1:r+1
    if bin(1,j) > 0 && bin(1,j) <= nbins,
      res(bin(1,j),1,j) = 1;
    end
    for n = 1-r+1:N+r,
      res(:,n+r,j) = res(:,n+r-1,j);
      if n > r+1 && bin(n-r-1,j) > 0 && bin(n-r-1,j) <= nbins,
        res(bin(n-r-1,j),n+r,j) = res(bin(n-r-1,j),n+r,j) - 1;
      end
      if n <= N-r && bin(n+r,j) > 0 && bin(n+r,j) <= nbins
        res(bin(n+r,j),n+r,j) = res(bin(n+r,j),n+r,j) + 1;
      end
    end
  end
  
  % normalize
  res(:,w:end-w+1,:) = res(:,w:end-w+1,:) / w;
  % boundary conditions
  res(:,1:w-1,:) = bsxfun(@rdivide,res(:,1:w-1,:),1:w-1);
  res(:,end-w+2:end,:) = bsxfun(@rdivide,res(:,end-w+2:end,:),w-1:-1:1);
  
  % offset
  windowis = find(windowi2radiusi == radiusi);
  for windowi = windowis',
    off = windows(windowi,2);
    % frame t, radius r, offset off:
    % [t-r+off, t+r+off]
    % so for r = 0, off = 1, we want [t+1,t+1]
    % which corresponds to res(t+r+off)
    % so we want to grab for 1+r+off through N+r+off
    res1 = padgrab(res,nan,1,nbins,1+r+off,N+r+off,1,ntrans);
    y(end+1:end+nbins,:) = res1(:,:,IDX.orig);
    for bini = 1:nbins,
      feature_names{end+1} = {'stat','hist','trans','none','radius',r,'offset',off,'bin',bini,'lim',hist_edges(bini:bini+1)}; %#ok<*AGROW>
    end
    
%     if ismember('none',trans_types),
%       
%       y(end+1:end+nbins,:) = res1(:,:,IDX.orig);
%       for bini = 1:nbins,
%         feature_names{end+1} = {'stat','hist','trans','none','radius',r,'offset',off,'bin',bini,'lim',hist_edges(bini:bini+1)}; %#ok<*AGROW>
%       end
%     end
%     
%     if IDX.abs > 0,
%       y(end+1:end+nbins,:) = res1(:,:,IDX.abs);
%       for bini = 1:nbins,
%         feature_names{end+1} = {'stat','hist','trans','abs','radius',r,'offset',off,'bin',bini,'lim',hist_edges(bini:bini+1)};
%       end
%     end
%     
%     if IDX.flip > 0 && ~( (r == 0) && (off == 0) && (IDX.abs > 0) ),
%       for bini = 1:nbins,
%         y(end+1,:) = res1(bini,:,IDX.orig);
%         y(end,x<0) = res1(bini,x<0,IDX.flip);
%         feature_names{end+1} = {'stat','hist','trans','flip','radius',r,'offset',off,'bin',bini,'lim',hist_edges(bini:bini+1)};
%       end
%     end
%     
%     if IDX.rel>0,
%       for bini = 1:nbins,
%         y(end+1,:) = res1(bini,:,IDX.rel);
%         feature_names{end+1} = {'stat','hist','trans','relative','radius',r,'offset',off,'bin',bini,'lim',hist_edges(bini:bini+1)};
%       end
%     end
    
    if SANITY_CHECK,
      
      if ismember('none',trans_types),
        fastY = res1(:,:,IDX.orig);
        res_real = fastY;
        res_dumb = nan(nbins,N);
        for n_dumb = 1:N,
          tmp = padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off);
          Z_dumb = nnz(~isnan(tmp));
          if all(isnan(tmp)),
            res_dumb(:,n_dumb) = nan;
          else
            tmp = histc(tmp,hist_edges);
            tmp = [tmp(1:end-2),tmp(end-1)+tmp(end)];
            tmp = tmp / Z_dumb;
            res_dumb(:,n_dumb) = tmp;
          end
        end
        checkSanity(res_real(:),res_dumb(:),r,off,'hist','none');
      end
      
%       if IDX.abs > 0,
%         fastY = res1(:,:,IDX.abs);
%         res_real = fastY;
%         res_dumb = nan(nbins,N);
%         for n_dumb = 1:N,
%           tmp = abs(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
%           Z_dumb = nnz(~isnan(tmp));
%           if all(isnan(tmp)),
%             res_dumb(:,n_dumb) = nan;
%           else
%             tmp = histc(tmp,hist_edges);
%             tmp = [tmp(1:end-2),tmp(end-1)+tmp(end)];
%             res_dumb(:,n_dumb) = tmp/Z_dumb;
%           end
%         end
%         checkSanity(res_real(:),res_dumb(:),r,off,'hist','abs');        
%       end
%       
%       if IDX.flip > 0 && ~( (r == 0) && (off == 0) && (IDX.abs > 0) ),
%         fastY = [];
%         for bini = 1:nbins,
%           fastY(bini,:) = res1(bini,:,IDX.orig);
%           fastY(bini,x<0) = res1(bini,x<0,IDX.flip);
%         end
%         res_real = fastY;
%         res_dumb = nan(nbins,N);
%         for n_dumb = 1:N,
%           m = 1;
%           if x(n_dumb) < 0,
%             m = -1;
%           end
%           tmp = m*padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off);
%           Z_dumb = nnz(~isnan(tmp));
%           if all(isnan(tmp)),
%             res_dumb(:,n_dumb) = nan;
%           else
%             tmp = histc(tmp,hist_edges);
%             tmp = [tmp(1:end-2),tmp(end-1)+tmp(end)];
%             res_dumb(:,n_dumb) = tmp/Z_dumb;
%           end
%         end
%         checkSanity(res_real(:),res_dumb(:),r,off,'hist','flip');
%       end
      
    end
  end
end
