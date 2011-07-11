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
trans_types = 'all';

% for debugging purposes
SANITY_CHECK = true;

% whether to cache results
DOCACHE = true;

% initialize empty cache
cache = InitializeCache();

% initialize feature_types already computed to empty
feature_types = {};

%% parse parameters

[...
  windows,...
  window_radii,window_offsets,...
  min_window_radius,max_window_radius,nwindow_radii,...
  feature_types,...
  trans_types,...
  SANITY_CHECK,...
  DOCACHE,...
  cache...
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'feature_types',feature_types,...
  'trans_types',trans_types,...
  'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE,...
  'cache',cache); %#ok<ASGLU>

%% whether we've specified to use all trans types by default
if ischar(trans_types) && strcmpi(trans_types,'all'),
  trans_types = {'none','abs','flip'}; %#ok<NASGU>
end

%% select default windows from various ways of specifying windows

[windows,window_radii,windowi2radiusi,nradii] = ...
  SetWindowParameters(...
  windows,window_offsets,...
  window_radii,min_window_radius,...
  max_window_radius,nwindow_radii);

%% main computation
  

for radiusi = 1:nradii,
  r = window_radii(radiusi);
  w = 2*r+1;
  
  % std will always be 0 for r == 0
  if r == 0,
    continue;
  end
  
  % average
  fil = ones(1,w);
  
  
  if DOCACHE && ismember(r,cache.std.radii),
    cache_i = find(r == cache.std.radii,1);
    res = cache.std.data{cache_i};
  else
    
    if DOCACHE && ismember(r,cache.mean.radii),
      cache_i = find(r == cache.mean.radii,1);
      res_mean = cache.mean.data{cache_i};
    else
      
      % full: res(t+r) corresponds to frame t
      res_mean = imfilter(x,fil,'full',0);
      % normalize
      res_mean(w:end-w+1) = res_mean(w:end-w+1) / w;
      % boundary conditions
      res_mean(1:w-1) = bsxfun(@rdivide,res_mean(1:w-1),1:w-1);
      res_mean(end-w+2:end) = bsxfun(@rdivide,res_mean(end-w+2:end),w-1:-1:1);
      
      if DOCACHE,
        cache.mean.radii(end+1) = r;
        cache.mean.data{end+1} = res_mean;
      end
      
    end
    
    % full: res(t+r) corresponds to frame t
    res = imfilter(x.^2,fil,'full',0);
    res(w:end-w+1) = res(w:end-w+1) / w;
    % boundary conditions
    res(1:w-1) = bsxfun(@rdivide,res(1:w-1),1:w-1);
    res(end-w+2:end) = bsxfun(@rdivide,res(end-w+2:end),w-1:-1:1);
    
    % combine to get standard deviation
    res = sqrt(res - res_mean.^2);
    
    if DOCACHE,
      cache.std.radii(end+1) = r;
      cache.std.data{end+1} = res;
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
    res1 = padgrab(res,nan,1,1,1+r+off,N+r+off);
    y(end+1,:) = res1; %#ok<*AGROW>
    feature_names{end+1} = {'stat','std','trans','none','radius',r,'offset',off};
    
    if SANITY_CHECK,
      
      res_dumb = nan(1,N);
      for n_dumb = 1:N,
        res_dumb(n_dumb) = nanstd(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off),1);
      end
      
      if any(isnan(y(end,:)) ~= isnan(res_dumb)),
        fprintf('SANITY CHECK: std, trans = none, r = %d, off = %d, nan mismatch\n',r,off);
      else
        fprintf('SANITY CHECK: std, trans = none, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
      end
      
    end
    
  end
end