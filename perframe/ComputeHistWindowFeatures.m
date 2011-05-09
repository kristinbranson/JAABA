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
  hist_edges...
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'feature_types',feature_types,...
  'trans_types',trans_types,...
  'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE,...
  'cache',cache,...
  'hist_edges',hist_edges); %#ok<ASGLU>

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

if isempty(hist_edges),
  warning('hist_edges is empty, but hist features are supposed to be computed. Skipping these features.');
  return;
end

% find bin for each value -- notice the transpose
nbins = numel(hist_edges)-1;
[~,bin] = histc(x_trans',hist_edges);
bin(bin > nbins) = nbins;

for radiusi = 1:nradii,
  r = window_radii(radiusi);
  w = 2*r+1;
  
  res = zeros(nbins,N+2*r,ntrans);
  for j = 1:ntrans,
    % initialize first value:
    % t contains [max(1,t-r),min(N,t+r)]
    % t = 1 will contain 1:r+1
    res(bin(1,j),1,j) = 1;
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
    y(end+1:end+nbins,:) = res1(:,:,IDX_ORIG);
    for bini = 1:nbins,
      feature_names{end+1} = {'stat','hist','trans','none','radius',r,'offset',off,'bin',bini,'lim',hist_edges(bini:bini+1)}; %#ok<*AGROW>
    end
    
    if SANITY_CHECK,
      
      res_real = y(end-nbins+1:end,:);
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
      if any(isnan(res_real(:)) ~= isnan(res_dumb(:))),
        warning('SANITY CHECK: hist, trans = none, r = %d, off = %d, nan mismatch\n',r,off);
      else
        fprintf('SANITY CHECK: hist, trans = none, r = %d, off = %d, max error = %f\n',r,off,max(abs(res_real(:)-res_dumb(:))));
      end
      
    end
    
    
    if IDX_ABS > 0,
      y(end+1:end+nbins,:) = res1(:,:,IDX_ABS);
      for bini = 1:nbins,
        feature_names{end+1} = {'stat','hist','trans','abs','radius',r,'offset',off,'bin',bini,'lim',hist_edges(bini:bini+1)};
      end
      
      
      if SANITY_CHECK,
        
        res_real = y(end-nbins+1:end,:);
        res_dumb = nan(nbins,N);
        for n_dumb = 1:N,
          tmp = abs(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
          Z_dumb = nnz(~isnan(tmp));
          if all(isnan(tmp)),
            res_dumb(:,n_dumb) = nan;
          else
            tmp = histc(tmp,hist_edges);
            tmp = [tmp(1:end-2),tmp(end-1)+tmp(end)];
            res_dumb(:,n_dumb) = tmp/Z_dumb;
          end
        end
        if any(isnan(res_real(:)) ~= isnan(res_dumb(:))),
          fprintf('SANITY CHECK: hist, trans = abs, r = %d, off = %d, nan mismatch\n',r,off);
        else
          fprintf('SANITY CHECK: hist, trans = abs, r = %d, off = %d, max error = %f\n',r,off,max(abs(res_real(:)-res_dumb(:))));
        end
        
      end
      
    end
    if IDX_FLIP > 0 && ~( (r == 0) && (off == 0) && (IDX_ABS > 0) ),
      for bini = 1:nbins,
        y(end+1,:) = res1(bini,:,IDX_ORIG);
        y(end,x<0) = res1(bini,x<0,IDX_FLIP);
        feature_names{end+1} = {'stat','hist','trans','flip','radius',r,'offset',off,'bin',bini,'lim',hist_edges(bini:bini+1)};
      end
      
      if SANITY_CHECK,
        
        res_real = y(end-nbins+1:end,:);
        res_dumb = nan(nbins,N);
        for n_dumb = 1:N,
          m = 1;
          if x(n_dumb) < 0,
            m = -1;
          end
          tmp = m*padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off);
          Z_dumb = nnz(~isnan(tmp));
          if all(isnan(tmp)),
            res_dumb(:,n_dumb) = nan;
          else
            tmp = histc(tmp,hist_edges);
            tmp = [tmp(1:end-2),tmp(end-1)+tmp(end)];
            res_dumb(:,n_dumb) = tmp/Z_dumb;
          end
        end
        if any(isnan(res_real(:)) ~= isnan(res_dumb(:))),
          fprintf('SANITY CHECK: hist, trans = flip, r = %d, off = %d, nan mismatch\n',r,off);
        else
          fprintf('SANITY CHECK: hist, trans = flip, r = %d, off = %d, max error = %f\n',r,off,max(abs(res_real(:)-res_dumb(:))));
        end
        
      end
      
    end
  end
end
