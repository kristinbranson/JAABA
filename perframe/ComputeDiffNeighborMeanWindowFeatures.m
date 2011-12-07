function [y,feature_names,cache] = ComputeDiffNeighborMeanWindowFeatures(x,varargin)

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
  relativeParams,...
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
if ischar(trans_types) && strcmpi(trans_types,'all'),
  trans_types = {'none','abs','flip','relative'};
end

%% select default windows from various ways of specifying windows

[windows,window_radii,windowi2radiusi,nradii] = ...
  SetWindowParameters(...
  windows,window_offsets,...
  window_radii,min_window_radius,...
  max_window_radius,nwindow_radii);

%% main computation

if ismember('relative',trans_types)
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
  
  % doesn't make sense for r == 0
  if r == 0,
    continue;
  end
  
  if DOCACHE && ismember(r,cache.mean.radii),
    cache_i = find(r == cache.mean.radii,1);
    res = cache.mean.data{cache_i};
  else
    res = MeanWindowCore(x,w);
    if DOCACHE,
      cache.mean.radii(end+1) = r;
      cache.mean.data{end+1} = res;
    end
  end
  
  if ismember('relative',trans_types)
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
    res1 = x - padgrab(res,nan,1,1,1+r+off,N+r+off);

    if ismember('none',trans_types),
      y(end+1,:) = res1; %#ok<*AGROW>
      feature_names{end+1} = {'stat','diff_neighbor_mean','trans','none','radius',r,'offset',off};
    end
    
    if ismember('abs',trans_types),
      y(end+1,:) = abs(res1);
      feature_names{end+1} = {'stat','diff_neighbor_mean','trans','abs','radius',r,'offset',off};
    end
    
    if ismember('flip',trans_types),
      res2 = res1.*sign(x);
      y(end+1,:) = res2;
      feature_names{end+1} = {'stat','diff_neighbor_mean','trans','flip','radius',r,'offset',off};
    end
    
    if ismember('relative',trans_types),
      resRel1 = modX - padgrab(resRel,nan,1,1,1+r+off,N+r+off);
      y(end+1,:) = resRel1;
      feature_names{end+1} = {'stat','diff_neighbor_mean','trans','relative','radius',r,'offset',off};
    end
    
    if SANITY_CHECK,
      funcType = 'DiffNeighborMean';
      if ismember('none',trans_types),
        fastY = res1; %#ok<*AGROW>
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = x(n_dumb) - nanmean(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
        end
        checkSanity(fastY,res_dumb,r,off,funcType,'none');
      end
    
      if ismember('abs',trans_types),
        fastY = abs(res1);
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = abs(x(n_dumb) - nanmean(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off)));
        end
        checkSanity(fastY,res_dumb,r,off,funcType,'abs');
      end
      
      if ismember('flip',trans_types),
        res2 = res1; res2(x<0) = -res2(x<0); fastY = res2;
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = x(n_dumb) - nanmean(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
          if x(n_dumb) < 0,
            res_dumb(n_dumb) = -res_dumb(n_dumb);
          end
        end
        checkSanity(fastY,res_dumb,r,off,funcType,'flip');
      end
      
    end
    
  end
end
