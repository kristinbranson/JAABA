% [y,feature_names,cache] = ComputeStdWindowFeatures(x,...)
% 
% Computes the standard deviation window features y for the input
% one-dimensional per-frame time series data x. y(i,t) will correspond to
% standard deviation window feature i and frame t, and is standard
% deviation of all the per-frame data in the window defined by i and t. 
% 
% For window feature i, let r_i be the window radius, and off_i be the
% window offset. The window feature y(i,t) is the standard deviation of the
% per-frame data between (t-r_i+off_i) and (t+r_i+off_i). The
% transformation types are ignored for this feature. 
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
% {'stat','std','trans','none','radius',1,'offset',1}
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
% type, if not otherwise specified. This input is actually ignored by this
% function. 
%
% 'sanitycheck': whether to compute all the features a second time in the
% obvious way to make sure that the optimized computations are correct.
% default value: false. 
%
% 'docache': whether to cache computations that might be useful, e.g. the
% mean computations are useful when computing the std window features.
% default value: true.
%
% 'cache': input cached computations that might be useful in this
% computation. 

function [y,feature_names,cache] = ComputeStdWindowFeatures(x,varargin)

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
%trans_types = 'all';
trans_types = uint8(15);

% for debugging purposes
SANITY_CHECK = true;

% whether to cache results
DOCACHE = true;

% initialize empty cache
cache = InitializeCache();

% initialize feature_types already computed to empty
feature_types = {};

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
  'relativeParams',relativeParams); %#ok<ASGLU>

%% whether we've specified to use all trans types by default
%if ischar(trans_types) && strcmpi(trans_types,'all'),
%  trans_types = {'none','abs','flip','relative'}; %#ok<NASGU>
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

for radiusi = 1:nradii,
  r = window_radii(radiusi);
  w = 2*r+1;
  
  % std will always be 0 for r == 0
  if r == 0,
    continue;
  end
  
  % average
  
  if DOCACHE && ismember(r,cache.std.radii),
    cache_i = find(r == cache.std.radii,1);
    res = cache.std.data{cache_i};
  else
    
    if DOCACHE && ismember(r,cache.mean.radii),
      cache_i = find(r == cache.mean.radii,1);
      res_mean = cache.mean.data{cache_i};
    else
      
      res_mean = MeanWindowCore(x,w);
%{      
%       % full: res(t+r) corresponds to frame t
%       res_mean = imfilter(x,fil,'full',0);
%       % normalize
%       res_mean(w:end-w+1) = res_mean(w:end-w+1) / w;
%       % boundary conditions
%       res_mean(1:w-1) = bsxfun(@rdivide,res_mean(1:w-1),1:w-1);
%       res_mean(end-w+2:end) = bsxfun(@rdivide,res_mean(end-w+2:end),w-1:-1:1);
%}
      
      if DOCACHE,
        cache.mean.radii(end+1) = r;
        cache.mean.data{end+1} = res_mean;
      end
      
    end
    
    res = MeanWindowCore(x.^2,w);
%{    
    % full: res(t+r) corresponds to frame t
%     res = imfilter(x.^2,fil,'full',0);
%     res(w:end-w+1) = res(w:end-w+1) / w;
%     % boundary conditions
%     res(1:w-1) = bsxfun(@rdivide,res(1:w-1),1:w-1);
%     res(end-w+2:end) = bsxfun(@rdivide,res(end-w+2:end),w-1:-1:1);
%}    
    % combine to get standard deviation
    diff = res - res_mean.^2;
    diff(diff<0) = 0;
    res = sqrt(diff);
    
    if DOCACHE,
      cache.std.radii(end+1) = r;
      cache.std.data{end+1} = res;
    end
    
  end

%  if ismember('relative',trans_types),
  if bitand(8,trans_types),
    
    if DOCACHE && ismember(r,cache.stdRel.radii),
      cache_i = find(r == cache.stdRel.radii,1);
      resRel = cache.stdRel.data{cache_i};
    else
      
      if DOCACHE && ismember(r,cache.meanRel.radii),
        cache_i = find(r == cache.meanRel.radii,1);
        resRel_mean = cache.meanRel.data{cache_i};
      else
        resRel_mean = MeanWindowCore(modX,w);
        if DOCACHE,
          cache.meanRel.radii(end+1) = r;
          cache.meanRel.data{end+1} = resRel_mean;
        end
      end
      
      resRel = MeanWindowCore(modX.^2,w);
      % combine to get standard deviation
      diff = resRel - resRel_mean.^2;
      diff(diff<0) = 0;
      resRel = sqrt(diff);
      
      if DOCACHE,
        cache.stdRel.radii(end+1) = r;
        cache.stdRel.data{end+1} = resRel;
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
    y(end+1,:) = res1; %#ok<*AGROW>
    feature_names{end+1} = {'stat','std','trans','none','radius',r,'offset',off};
    
%    if ismember('relative',trans_types)
    if bitand(8,trans_types)
      resRel1 = padgrab2(resRel,nan,1,1,1+r+off,N+r+off);
      y(end+1,:) = resRel1; %#ok<*AGROW>
      feature_names{end+1} = {'stat','std','trans','relative','radius',r,'offset',off};
    end
    
    if SANITY_CHECK,
      res_dumb = nan(1,N);
      for n_dumb = 1:N,
        res_dumb(n_dumb) = nanstd(padgrab2(x,nan,1,1,n_dumb-r+off,n_dumb+r+off),1);
      end
      checkSanity(y(end,:),res_dumb(:),r,off,'std','none');
    end
    
  end
end
