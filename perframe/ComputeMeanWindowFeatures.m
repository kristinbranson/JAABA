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
  'cache',cache);

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

% loop over window radii
for radiusi = 1:nradii,

  r = window_radii(radiusi);
  w = 2*r+1;

  % no need to do mean for r == 0 if we've done min or max
  if r == 0 && any(ismember({'min','max'},feature_types)),
    continue;
  end

  % is this already in the cache?
  if DOCACHE && ismember(r,cache.mean.radii),
    cache_i = find(r == cache.mean.radii,1);
    res = cache.mean.data{cache_i};
  else
    
    % average
    fil = ones(1,w);
    % full: res(t+r) corresponds to frame t
    res = imfilter(x,fil,'full',0);
    % normalize
    res(w:end-w+1) = res(w:end-w+1) / w;
    % boundary conditions
    res(1:w-1) = bsxfun(@rdivide,res(1:w-1),1:w-1);
    res(end-w+2:end) = bsxfun(@rdivide,res(end-w+2:end),w-1:-1:1);
    
    % store for future computations
    if DOCACHE,
      cache.mean.radii(end+1) = r;
      cache.mean.data{end+1} = res;
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
    feature_names{end+1} = {'stat','mean','trans','none','radius',r,'offset',off};    
    
    if SANITY_CHECK,
      
      res_dumb = nan(1,N);
      for n_dumb = 1:N,
        res_dumb(n_dumb) = nanmean(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
      end
      
      if any(isnan(y(end,:)) ~= isnan(res_dumb)),
        fprintf('SANITY CHECK: mean, trans = none, r = %d, off = %d, nan mismatch\n',r,off);
      else
        fprintf('SANITY CHECK: mean, trans = none, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
      end
      
    end
    
    if ismember('abs',trans_types),
      y(end+1,:) = abs(res1);
      feature_names{end+1} = {'stat','mean','trans','abs','radius',r,'offset',off};
      
      if SANITY_CHECK,
        
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = abs(nanmean(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off)));
        end
        
        if any(isnan(y(end,:)) ~= isnan(res_dumb)),
          fprintf('SANITY CHECK: mean, trans = abs, r = %d, off = %d, nan mismatch\n',r,off);
        else
          fprintf('SANITY CHECK: mean, trans = abs, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
        end
      end
      
      
    end
    % flip is redundant with abs if r = 0 && off = 0
    if ismember('flip',trans_types) && ~( (r == 0) && (off == 0) && ismember('abs',trans_types) ),
      y(end+1,:) = res1;
      y(end,x<0) = -res1(x<0);
      feature_names{end+1} = {'stat','mean','trans','flip','radius',r,'offset',off};
      
      if SANITY_CHECK,
        
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = nanmean(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
          if sign(x(n_dumb)) < 0,
            res_dumb(n_dumb) = -res_dumb(n_dumb);
          end
        end
        
        if any(isnan(y(end,:)) ~= isnan(res_dumb)),
          fprintf('SANITY CHECK: mean, trans = none, r = %d, off = %d, nan mismatch\n',r,off);
        else
          fprintf('SANITY CHECK: mean, trans = flip, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
        end
      end
      
    end
  end
end