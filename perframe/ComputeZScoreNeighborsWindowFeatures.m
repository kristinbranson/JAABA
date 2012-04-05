function [y,feature_names,cache] = ComputeZscoreNeighborsWindowFeatures(x,varargin)
  
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

for radiusi = 1:nradii,
  r = window_radii(radiusi);
  
  % doesn't make sense for r == 0
  if r == 0,
    continue;
  end
  
  w = 2*r+1;
  
  if DOCACHE && ismember(r,cache.mean.radii),
    cache_i = find(r == cache.mean.radii,1);
    res_mean = cache.mean.data{cache_i};
  else
    res_mean = MeanWindowCore(x,w);
    if DOCACHE,
      cache.mean.radii(end+1) = r;
      cache.mean.data{end+1} = res_mean;
    end
  end
  
  if DOCACHE && ismember(r,cache.std.radii),
    cache_i = find(r == cache.std.radii,1);
    res_std = cache.std.data{cache_i};
  else
    % full: res_std(t+r) corresponds to frame t
    res_std = MeanWindowCore(x.^2,w);    
    % combine to get standard deviation
    res_std = sqrt(res_std - res_mean.^2);
    if DOCACHE,
      cache.std.radii(end+1) = r;
      cache.std.data{end+1} = res_std;
    end
  end
  res_std(res_std==0) = 1;
  
  
  %if ismember('relative',trans_types),
  if bitand(8,trans_types),
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
    
    if DOCACHE && ismember(r,cache.stdRel.radii),
      cache_i = find(r == cache.stdRel.radii,1);
      resRel_std = cache.stdRel.data{cache_i};
    else
      % full: res_std(t+r) corresponds to frame t
      resRel_std = MeanWindowCore(modX.^2,w);
      % combine to get standard deviation
      resRel_std = sqrt(resRel_std - resRel_mean.^2);
      if DOCACHE,
        cache.stdRel.radii(end+1) = r;
        cache.stdRel.data{end+1} = resRel_std;
      end
    end
    resRel_std(resRel_std==0) = 1;
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
    res1 = (x - padgrab(res_mean,nan,1,1,1+r+off,N+r+off))./padgrab(res_std,nan,1,1,1+r+off,N+r+off);
    
    %if ismember('none',trans_types),
    if bitand(1,trans_types),
      y(end+1,:) = res1; %#ok<*AGROW>
      feature_names{end+1} = {'stat','zscore_neighbors','trans','none','radius',r,'offset',off};
    end
    
    %if ismember('abs',trans_types),
    if bitand(2,trans_types),
      y(end+1,:) = abs(res1);
      feature_names{end+1} = {'stat','zscore_neighbors','trans','abs','radius',r,'offset',off};
    end
    
    %if ismember('flip',trans_types),
    if bitand(4,trans_types),
      res2 = res1;
      res2(x<0) = -res2(x<0);
      y(end+1,:) = res2;
      feature_names{end+1} = {'stat','zscore_neighbors','trans','flip','radius',r,'offset',off};
    end
    
    %if ismember('relative',trans_types),
    if bitand(8,trans_types),
      resRel1 = (modX - padgrab(resRel_mean,nan,1,1,1+r+off,N+r+off))./padgrab(resRel_std,nan,1,1,1+r+off,N+r+off);
      y(end+1,:) = resRel1; %#ok<*AGROW>
      feature_names{end+1} = {'stat','zscore_neighbors','trans','relative','radius',r,'offset',off};
    end
    
    if SANITY_CHECK,
      
      %if ismember('none',trans_types),
      if bitand(1,trans_types),
        fastY = res1; %#ok<*AGROW>
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          s = nanstd(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off),1);
          if s == 0, s = 1; end
          res_dumb(n_dumb) = (x(n_dumb) - nanmean(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off)))/s;
        end
        checkSanity(fastY,res_dumb,r,off,'zscore_neighbor','none');        
      end
    
      %if ismember('abs',trans_types),
      if bitand(2,trans_types),
        fastY = abs(res1);
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          s = nanstd(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off),1);
          if s == 0, s = 1; end
          res_dumb(n_dumb) = abs(x(n_dumb) - nanmean(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off)))/s;
        end
        checkSanity(fastY,res_dumb,r,off,'zscore_neighbor','abs');
      end
    
      %if ismember('flip',trans_types),   %%%  no relative here as above?
      if bitand(4,trans_types),
        res2 = res1;
        res2(x<0) = -res2(x<0);
        fastY = res2;
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          s = nanstd(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off),1);
          if s == 0,s = 1;end
          res_dumb(n_dumb) = (x(n_dumb) - nanmean(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off)))/s;
          if x(n_dumb) < 0,
            res_dumb(n_dumb) = -res_dumb(n_dumb);
          end
        end
        checkSanity(fastY,res_dumb,r,off,'zscore_neighbor','flip');
      end
      
    end
    
  end
end
