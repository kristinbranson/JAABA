% [y,feature_names,cache] = ComputeChangeWindowFeatures(x,...)
% 
% Computes the change window features y for the input one-dimensional
% per-frame time series data x. y(i,t) will correspond to change
% window feature i and frame t, and is the (transformation of) the
% difference between the average of the per-frame data at the end and start
% of the window defined by i and t. 
% 
% For window feature i, let rho_i be the change window radius, r_i be the
% window radius, and off_i be the window offset. If the transformation type
% for feature i is 'none', the window feature y(i,t) is the difference
% between the average per-frame data at the end of the window (from
% (t+r_i+off_i)-rho_i to (t+r_i+off_i)+rho_i) and at the start of the
% window (from (t-r_i+off_i)-rho_i to (t-r_i+off_i)+rho_i). If the
% transformation type for feature i is 'abs', then we compute the absolute
% value of this difference. If the transformation type for feature
% i is 'flip', then we compute the sign of x(t) times this difference.
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
% {'stat','change','trans','abs','radius',1,'offset',1,'change_window_radius',0}
% cache: if DOCACHE is set to true, then computations useful in other
% window feature calculations may be saved in the output cache. 
%
% Optional inputs:
%
% 'change_window_radii': The 1 X NCHANGEWINDOWRADII change window radii to
% use (this defines the size of the sub-window averaged over at the end and
% start of the original window). default value: 0. 
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
% mean computations are useful when computing the change window features.
% default value: true.
%
% 'cache': input cached computations that might be useful in this
% computation. 

function [y,feature_names,cache] = ComputeChangeWindowFeatures(x,varargin)

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
%feature_types = {};

% change window radii (width = 2*radius + 1)
change_window_radii = 0;

%% parse parameters

[...
  windows,...
  window_radii,window_offsets,...
  min_window_radius,max_window_radius,nwindow_radii,...
  trans_types,...
  SANITY_CHECK,...
  DOCACHE,...
  cache,...
  change_window_radii...
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'trans_types',trans_types,...
  'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE,...
  'cache',cache,...
  'change_window_radii',change_window_radii); 
%   'feature_types',feature_types,...

%% whether we've specified to use all trans types by default
if ischar(trans_types) && strcmpi(trans_types,'all'),
  trans_types = {'none','abs','flip'};
end

%% select default windows from various ways of specifying windows

[windows,window_radii,windowi2radiusi,nradii] = ...
  SetWindowParameters(...
  windows,window_offsets,...
  window_radii,min_window_radius,...
  max_window_radius,nwindow_radii);

%% main computation

for change_r_i = 1:numel(change_window_radii),
  
  % take the mean across all windows of radius change_r
  change_r = change_window_radii(change_r_i);
  change_w = 2*change_r+1;
  
  if DOCACHE && ismember(change_r,cache.mean.radii),
    cache_i = find(change_r == cache.mean.radii,1);
    res_mean = cache.mean.data{cache_i};
  else
    % average filter
    fil = ones(1,change_w);
    % full: res_mean(t+change_r) corresponds to frame t
    res_mean = imfilter(x,fil,'full',0);
    
    % normalize
    res_mean(change_w:end-change_w+1) = res_mean(change_w:end-change_w+1) / change_w;
    % boundary conditions
    res_mean(1:change_w-1) = bsxfun(@rdivide,res_mean(1:change_w-1),1:change_w-1);
    res_mean(end-change_w+2:end) = bsxfun(@rdivide,res_mean(end-change_w+2:end),change_w-1:-1:1);
    
    % store for future computations
    if DOCACHE,
      cache.mean.radii(end+1) = change_r;
      cache.mean.data{end+1} = res_mean;
    end
    
  end
  
  % loop over window radii
  for radiusi = 1:nradii,
    r = window_radii(radiusi);
    % don't do for r <= change_r -- the windows will overlap, cancel, ...
    if r <= change_r,
      continue;
    end
    w = 2*r+1;
    
    % take the difference between the end of the window and the start of
    % the window
    % res(t) corresponds to t+r-change_r
    % so res(t-r+change_r) corresponds to frame t
    res = res_mean(w:end) - res_mean(1:end-w+1);
    
    % offset
    windowis = find(windowi2radiusi == radiusi);
    for windowi = windowis',
      
      off = windows(windowi,2);
      % res(t-r+change_r) corresponds to frame t,
      % and we want to grab from [1+off,N+off] which is
      % [1+off-r+change_r,N+off-r+change_r], relative to res
      res1 = padgrab(res,nan,1,1,1+off-r+change_r,N+off-r+change_r)/r;
      
      if ismember('none',trans_types),
        
        y(end+1,:) = res1;
        feature_names{end+1} = {'stat','change','trans','none','radius',r,'offset',off,'change_window_radius',change_r}; %#ok<*AGROW>
        
        if SANITY_CHECK,
          
          res_real = y(end,:);
          res_dumb = nan(1,N);
          for n_dumb = 1:N,
            tmp1 = padgrab(x,nan,1,1,n_dumb+off-r-change_r,n_dumb+off-r+change_r);
            tmp2 = padgrab(x,nan,1,1,n_dumb+off+r-change_r,n_dumb+off+r+change_r);
            res_dumb(n_dumb) = (nanmean(tmp2) - nanmean(tmp1))/r;
          end
          if any(isnan(res_real(:)) ~= isnan(res_dumb(:))),
            fprintf('SANITY CHECK: change, trans = none, r = %d, off = %d, change_r = %d, nan mismatch\n',r,off,change_r);
          else
            fprintf('SANITY CHECK: change, trans = none, r = %d, off = %d, change_r = %d, max error = %f\n',r,off,change_r,max(abs(res_real(:)-res_dumb(:))));
          end
          
        end
        
      end
      
      if ismember('abs',trans_types),
        y(end+1,:) = abs(res1);
        feature_names{end+1} = {'stat','change','trans','abs','radius',r,'offset',off,'change_window_radius',change_r};
        
        if SANITY_CHECK,
          
          res_real = y(end,:);
          res_dumb = nan(1,N);
          for n_dumb = 1:N,
            tmp1 = padgrab(x,nan,1,1,n_dumb+off-r-change_r,n_dumb+off-r+change_r);
            tmp2 = padgrab(x,nan,1,1,n_dumb+off+r-change_r,n_dumb+off+r+change_r);
            res_dumb(n_dumb) = abs(nanmean(tmp2) - nanmean(tmp1))/r;
          end
          if any(isnan(res_real(:)) ~= isnan(res_dumb(:))),
            fprintf('SANITY CHECK: change, trans = abs, r = %d, off = %d, change_r = %d, nan mismatch\n',r,off,change_r);
          else
            fprintf('SANITY CHECK: change, trans = abs, r = %d, off = %d, change_r = %d, max error = %f\n',r,off,change_r,max(abs(res_real(:)-res_dumb(:))));
          end
          
        end
        
      end
      
      if ismember('flip',trans_types)
        res2 = res1;
        res2(x<0) = -res2(x<0);
        y(end+1,:) = res2;
        feature_names{end+1} = {'stat','change','trans','flip','radius',r,'offset',off,'change_window_radius',change_r};
        
        if SANITY_CHECK,
          
          res_real = y(end,:);
          res_dumb = nan(1,N);
          for n_dumb = 1:N,
            tmp1 = padgrab(x,nan,1,1,n_dumb+off-r-change_r,n_dumb+off-r+change_r);
            tmp2 = padgrab(x,nan,1,1,n_dumb+off+r-change_r,n_dumb+off+r+change_r);
            res_dumb(n_dumb) = (nanmean(tmp2) - nanmean(tmp1))/r;
            if x(n_dumb) < 0,
              res_dumb(n_dumb) = -res_dumb(n_dumb);
            end
          end
          if any(isnan(res_real(:)) ~= isnan(res_dumb(:))),
            fprintf('SANITY CHECK: change, trans = abs, r = %d, off = %d, change_r = %d, nan mismatch\n',r,off,change_r);
          else
            fprintf('SANITY CHECK: change, trans = abs, r = %d, off = %d, change_r = %d, max error = %f\n',r,off,change_r,max(abs(res_real(:)-res_dumb(:))));
          end
          
        end
        
      end
      
    end
  end
  
end
