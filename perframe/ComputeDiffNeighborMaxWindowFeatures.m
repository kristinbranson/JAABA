function [y,feature_names,cache] = ComputeDiffNeighborMaxWindowFeatures(x,varargin)
  
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
  trans_types = {'none','abs','flip'};
end

%% select default windows from various ways of specifying windows

[windows,window_radii,windowi2radiusi,nradii] = ...
  SetWindowParameters(...
  windows,window_offsets,...
  window_radii,min_window_radius,...
  max_window_radius,nwindow_radii);

%% compute per-frame transformations 
[x_trans,IDX_ORIG,IDX_ABS,IDX_FLIP,ntrans] = ComputePerFrameTrans(x,trans_types);

%% main computation

for radiusi = 1:nradii,
  r = window_radii(radiusi);
  
  % doesn't make sense for r == 0
  if r == 0,
    continue;
  end
  
  w = 2*r+1;
  
  if DOCACHE && ismember(r,cache.max.radii),
    cache_i = find(r == cache.max.radii,1);
    res = cache.max.data{cache_i};
  else
    % pad with infs for boundary conditions
    x_pad = [-inf(ntrans,r),x_trans,-inf(ntrans,r)];
    % maximum: use imdilate
    se = strel(ones(1,w));
    res = imdilate(x_pad,se);
    
    if DOCACHE,
      cache.max.radii(end+1) = r;
      cache.max.data{end+1} = res;
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
    res1 = padgrab(res,nan,1,ntrans,1+r+off,N+r+off) - x_trans;
    y(end+1,:) = res1(IDX_ORIG,:); %#ok<*AGROW>
    feature_names{end+1} = {'stat','diff_neighbor_max','trans','none','radius',r,'offset',off};
    
    if SANITY_CHECK,
      
      res_dumb = nan(1,N);
      for n_dumb = 1:N,
        res_dumb(n_dumb) = max(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off)) - x(n_dumb);
      end
      
      if any(isnan(y(end,:)) ~= isnan(res_dumb)),
        fprintf('SANITY CHECK: diff_neighbor_max, trans = none, r = %d, off = %d, nan mismatch\n',r,off);
      else
        fprintf('SANITY CHECK: diff_neighbor_max, trans = none, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
      end
      
    end
    
    
    if IDX_ABS > 0,
      y(end+1,:) = res1(IDX_ABS,:);
      feature_names{end+1} = {'stat','diff_neighbor_max','trans','abs','radius',r,'offset',off};
      
      
      if SANITY_CHECK,
        
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = nanmax(abs(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off))) - abs(x(n_dumb));
        end
        
        if any(isnan(y(end,:)) ~= isnan(res_dumb)),
          fprintf('SANITY CHECK: diff_neighbor_max, trans = abs, r = %d, off = %d, nan mismatch\n',r,off);
        else
          fprintf('SANITY CHECK: diff_neighbor_max, trans = abs, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
        end
        
      end
      
    end
    % flip is redundant with abs if r = 0 && off = 0
    if IDX_FLIP > 0 && ~( (r == 0) && (off == 0) && (IDX_ABS > 0) ),
      y(end+1,:) = res1(IDX_ORIG,:);
      y(end,x<0) = res1(IDX_FLIP,x<0);
      feature_names{end+1} = {'stat','diff_neighbor_max','trans','flip','radius',r,'offset',off};
      
      if SANITY_CHECK,
        
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          m_dumb = 1;
          if sign(x(n_dumb)) < 0,
            m_dumb = -1;
          end
          res_dumb(n_dumb) = nanmax(m_dumb*padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off)) - abs(x(n_dumb));
        end
        
        if any(isnan(y(end,:)) ~= isnan(res_dumb)),
          fprintf('SANITY CHECK: diff_neighbor_max, trans = flip, r = %d, off = %d, nan mismatch\n',r,off);
        else
          fprintf('SANITY CHECK: diff_neighbor_max, trans = flip, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
        end
        
      end
      
    end
    
  end
  
  
  
end
