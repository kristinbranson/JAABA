% [y,feature_names,cache] = ComputeMeanWindowFeatures(x,...)
% 
% Computes the mean window features y for the input one-dimensional
% per-frame time series data x. y(i,t) will correspond to mean
% window feature i and frame t, and is the (transformation of) the mean of
% all the per-frame data within the window defined by i at t. 
% 
% Let r_i be the window radius for feature i and off_i be the window offset
% for feature i. If the transformation for feature i is 'none', the window
% feature y(i,t) is the mean of the per-frame data within the window from
% t-r_i+off_i through t+r_i+off_i. If the transformation for feature i is
% 'abs', then the window feature y(i,t) is the absolute value of the mean
% of the per-frame data within this same window. If the transformation for
% feature i is 'flip', then the window feature y(i,t) is the sign of x(t)
% times the mean of the per-frame data within the window. This is
% equivalent to the absolute value feature if radius == 0 and offset == 0,
% so this feature will not be added if it is guaranteed to be totally
% redundant.
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
% {'stat','mean','trans','abs','radius',1,'offset',1}
% cache: if DOCACHE is set to true, then the non-offset window means will
% be cached to potentially be used in future window feature computations. 
%
% Optional inputs:
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
% mean computations are useful when computing the diff_neighbor_mean.
% default value: true.
%
% 'cache': input cached computations that might be useful in this
% computation. 

function [y,feature_names,cache] = ComputeMeanWindowFeatures(x,varargin)

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
trans_types = uint8(15);

% for debugging purposes
SANITY_CHECK = true;

% whether to cache results
DOCACHE = true;

% initialize empty cache
cache = InitializeCache();

relativeParams = [];

% initialize feature_types already computed to empty
%feature_types = {};

%% parse parameters

[...
  windows,...
  window_radii,window_offsets,...
  min_window_radius,max_window_radius,nwindow_radii,...
  trans_types,...
  SANITY_CHECK,...
  DOCACHE,...
  cache,...
  relativeParams,...
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'trans_types',trans_types,...
  'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE,...
  'cache',cache,...
  'relativeParams',relativeParams);
%  'feature_types',feature_types,...

%% whether we've specified to use all trans types by default
%if ischar(trans_types) && strcmpi(trans_types,'all'),
%  trans_types = {'none','abs','flip','relative'};
%end

%% select default windows from various ways of specifying windows

[windows,window_radii,windowi2radiusi,nradii] = ...
  SetWindowParameters(...
  windows,window_offsets,...
  window_radii,min_window_radius,...
  max_window_radius,nwindow_radii);

%% main computation

%if ismember('relative',trans_types)
if bitand(8,trans_types)
  if DOCACHE && ~isempty(cache.relX)
    modX = cache.relX;
  else
    modX = convertToRelative(x,relativeParams);
    cache.relX = modX;
  end
end

% loop over window radii
for radiusi = 1:nradii,

  r = window_radii(radiusi);
  w = 2*r+1;

%   % no need to do mean for r == 0 if we've done min or max
%   if r == 0 && any(ismember({'min','max'},feature_types)),
%     continue;
%   end

  % is this already in the cache?
  if DOCACHE && ismember(r,cache.mean.radii),
    cache_i = find(r == cache.mean.radii,1);
    res = cache.mean.data{cache_i};
  else
    res = MeanWindowCore(x,w);
    % store for future computations
    if DOCACHE,
      cache.mean.radii(end+1) = r;
      cache.mean.data{end+1} = res;
    end
  end
  
%  if ismember('relative',trans_types)
  if bitand(8,trans_types)
    if DOCACHE && ismember(r,cache.meanRel.radii),
      cache_i = find(r == cache.meanRel.radii,1);
      resRel = cache.meanRel.data{cache_i};
    else
      resRel = MeanWindowCore(modX,w);
      % store for future computations
      if DOCACHE,
        cache.meanRel.radii(end+1) = r;
        cache.meanRel.data{end+1} = resRel;
      end
    end
  end
  
  % all offsets for this radius
  windowis = find(windowi2radiusi == radiusi);
  for windowi = windowis',
    off = windows(windowi,2);
    % frame t, radius r, offset off:
    % [t-r+off, t+r+off]
    % so for r = 0, off = 1, we want [t+1,t+1]
    % which corresponds to res(t+r+off)
    % so we want to grab for 1+r+off through N+r+off
    res1 = padgrab2(res,nan,1,1,1+r+off,N+r+off);
%    if ismember('none',trans_types),
    if bitand(1,trans_types),
      y(end+1,:) = res1; %#ok<*AGROW>
      feature_names{end+1} = {'stat','mean','trans','none','radius',r,'offset',off};
    end
    
%    if ismember('abs',trans_types),
    if bitand(2,trans_types),
      y(end+1,:) = abs(res1);
      feature_names{end+1} = {'stat','mean','trans','abs','radius',r,'offset',off};
    end
    
%    if ismember('flip',trans_types) && ~( (r == 0) && (off == 0) && ismember('abs',trans_types) ),
    if bitand(4,trans_types) && ~( (r == 0) && (off == 0) && bitand(2,trans_types) ),
      y(end+1,:) = res1;
      y(end,x<0) = -res1(x<0);
      feature_names{end+1} = {'stat','mean','trans','flip','radius',r,'offset',off};
    end
    
%    if ismember('relative',trans_types),
    if bitand(8,trans_types),
      resRel1 = padgrab2(resRel,nan,1,1,1+r+off,N+r+off);
      y(end+1,:) = resRel1;
      feature_names{end+1} = {'stat','mean','trans','relative','radius',r,'offset',off};
    end
    
    if SANITY_CHECK,
      funcType = 'mean';
      
%      if ismember('none',trans_types),
      if bitand(1,trans_types),
        fastY = res1;
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = nanmean(padgrab2(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
        end
        checkSanity(fastY,res_dumb,r,off,funcType);
      end  
      
%      if ismember('abs',trans_types),
      if bitand(2,trans_types),
        fastY = abs(res1);
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = abs(nanmean(padgrab2(x,nan,1,1,n_dumb-r+off,n_dumb+r+off)));
        end
        checkSanity(fastY,res_dumb,r,off,funcType);
      end
      
      % flip is redundant with abs if r = 0 && off = 0
%      if ismember('flip',trans_types) && ~( (r == 0) && (off == 0) && ismember('abs',trans_types) ),
      if bitand(4,trans_types) && ~( (r == 0) && (off == 0) && bitand(2,trans_types) ),
        fastY = res1;
        fastY(x<0) = -res1(x<0);
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = nanmean(padgrab2(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
          if sign(x(n_dumb)) < 0,
            res_dumb(n_dumb) = -res_dumb(n_dumb);
          end
        end
        checkSanity(fastY,res_dumb,r,off,funcType);
        
      end
      
    end % End Sanity Check
  end
end
