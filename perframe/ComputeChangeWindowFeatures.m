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
feature_types = {};

% change window radii (width = 2*radius + 1)
change_window_radii = 0;

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
  change_window_radii...
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'feature_types',feature_types,...
  'trans_types',trans_types,...
  'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE,...
  'cache',cache,...
  'change_window_radii',change_window_radii); %#ok<ASGLU>

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
