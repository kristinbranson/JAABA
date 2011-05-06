function [y,feature_names] = ComputeWindowFeatures(x,varargin)

%% initialize

x = x(:)';
N = numel(x);
y = nan(0,N);
feature_names = {};

%% default parameters

% defaults parameters for windows

% specify window using cross product of all radii and all offsets
windows = []; 

window_offsets = [-1,0,1];

% specify radii using minimum, maximum radii, and number of window radii
window_radii = [];
min_window_radius = 0;
max_window_radius = 20;
nwindow_radii = 5;

% use all features, transformation types by default
feature_types = 'all';
trans_types = 'all';

% edges of histogram bins
hist_edges = [];

% percentiles to compute (other than min, max)
prctiles = [];

% change window radii (width = 2*radius + 1)
change_window_radii = 0;

% harmonic features
num_harmonic = 1;

% for debugging purposes
SANITY_CHECK = true;

%% parse parameters

[...
  windows,...
  window_radii,window_offsets,...
  min_window_radius,max_window_radius,nwindow_radii,...
  feature_types,trans_types,...
  hist_edges,...
  prctiles,...
  SANITY_CHECK,...
  change_window_radii,...
  harmonic_fractions
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'feature_types',feature_types,'trans_types',trans_types,...
  'hist_edges',hist_edges,...
  'prctiles',prctiles,...
  'sanitycheck',SANITY_CHECK,...
  'change_window_radii',change_window_radii,...
  'num_harmonic',num_harmonic);

%% whether we've specified to use all features, trans

useallfeatures = ischar(feature_types) && strcmpi(feature_types,'all');
usealltrans = ischar(trans_types) && strcmpi(trans_types,'all');

%% select windows
if isempty(windows),
  
  % use radii and offsets if windows not input
  if isempty(window_radii),
    
    % use min, max, n if radii not input
    window_radii = unique(round(logspace(log10(min_window_radius+1),log10(max_window_radius+1),nwindow_radii)))-1;
  end
  if isempty(window_radii),
    error('window_radii is empty.');
  end
  if isempty(window_offsets),
    error('window_offsets is empty.');
  end
  
  % take all combinations of radii and offsets
  [all_radii,all_offsets] = meshgrid(window_radii,window_offsets);
  all_radii = all_radii(:);
  all_offsets = all_offsets(:);
  
  % offsets are fractions of radii
  all_offsets = round(all_offsets.*(max(all_radii,1)));

  % window is the pair of radius, offset
  windows = [all_radii,all_offsets];

  % make sure these are unique -- rounding to nearest frame might make them
  % not unique
  windows = unique(windows,'rows');
  
  % window i for frame t is
  % [t-windows(i,1)+windows(i,2), t+windows(i,1)+windows(i,2)]
  
end

if isempty(windows),
  error('windows is empty.');
end
if size(windows,2) ~= 2,
  error('windows must be nwindows x 2.');
end

[window_radii,~,windowi2radiusi] = unique(windows(:,1));
nradii = numel(window_radii);

%% compute transformations
x_trans = x;
IDX_ORIG = 1;
IDX_ABS = 0;
IDX_FLIP = 0;
if usealltrans || ismember('abs',trans_types),
  x_trans(end+1,:) = abs(x);
  IDX_ABS = size(x_trans,1);
end
if usealltrans || ismember('signflip',trans_types),
  x_trans(end+1,:) = -x;
  x_isneg = sign(x) < 0;
  IDX_FLIP = size(x_trans,1);
end
ntrans = size(x_trans,1);

%% window mean

if useallfeatures || ismember('mean',feature_types),
  for radiusi = 1:nradii,
    r = window_radii(radiusi);
    w = 2*r+1;
    % average
    fil = ones(1,w);
    % full: res(t+r) corresponds to frame t
    res = imfilter(x,fil,'full',0);
    % normalize
    res(w:end-w+1) = res(w:end-w+1) / w;
    % boundary conditions
    res(1:w-1) = bsxfun(@rdivide,res(1:w-1),1:w-1);
    res(end-w+2:end) = bsxfun(@rdivide,res(end-w+2:end),w-1:-1:1);
    
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
      
      if IDX_ABS > 0,
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
      if IDX_FLIP > 0 && ~( (r == 0) && (off == 0) && (IDX_ABS > 0) ), 
        y(end+1,:) = res1;
        y(end,x_isneg) = -res1(x_isneg);
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
end

%% window minimum

if useallfeatures || ismember('min',feature_types),
  for radiusi = 1:nradii,
    r = window_radii(radiusi);
    % no need to do minimum, maximum for r == 0 if we've done mean
    if r == 0 && (useallfeatures || any(ismember({'mean','max'},feature_types))),
      continue;
    end
    w = 2*r+1;
    % pad with infs for boundary conditions
    x_pad = [inf(ntrans,r),x_trans,inf(ntrans,r)];
    % minimum: use imerode
    se = strel(ones(1,w));
    res = imerode(x_pad,se);
    % offset
    windowis = find(windowi2radiusi == radiusi);
    for windowi = windowis',
      off = windows(windowi,2);
      % frame t, radius r, offset off:
      % [t-r+off, t+r+off]
      % so for r = 0, off = 1, we want [t+1,t+1]
      % which corresponds to res(t+r+off)
      % so we want to grab for 1+r+off through N+r+off
      res1 = padgrab(res,nan,1,ntrans,1+r+off,N+r+off);
      y(end+1,:) = res1(IDX_ORIG,:);
      feature_names{end+1} = {'stat','min','trans','none','radius',r,'offset',off};
      
      if SANITY_CHECK,
        
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = nanmin(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
        end
        
        if any(isnan(y(end,:)) ~= isnan(res_dumb)),
          fprintf('SANITY CHECK: min, trans = none, r = %d, off = %d, nan mismatch\n',r,off);
        else
          fprintf('SANITY CHECK: min, trans = none, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
        end
        
      end
      
      
      if IDX_ABS > 0,
        y(end+1,:) = res1(IDX_ABS,:);
        feature_names{end+1} = {'stat','min','trans','abs','radius',r,'offset',off};
        
        
        if SANITY_CHECK,
          
          res_dumb = nan(1,N);
          for n_dumb = 1:N,
            res_dumb(n_dumb) = nanmin(abs(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off)));
          end
          
          if any(isnan(y(end,:)) ~= isnan(res_dumb)),
            fprintf('SANITY CHECK: min, trans = abs, r = %d, off = %d, nan mismatch\n',r,off);
          else
            fprintf('SANITY CHECK: min, trans = abs, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
          end
          
        end
        
      end
      % flip is redundant with abs if r = 0 && off = 0
      if IDX_FLIP > 0 && ~( (r == 0) && (off == 0) && (IDX_ABS > 0) ), 
        y(end+1,:) = res1(IDX_ORIG,:);
        y(end,x_isneg) = res1(IDX_FLIP,x_isneg);
        feature_names{end+1} = {'stat','min','trans','flip','radius',r,'offset',off};
        
        if SANITY_CHECK,
          
          res_dumb = nan(1,N);
          for n_dumb = 1:N,
            m_dumb = 1;
            if sign(x(n_dumb)) < 0,
              m_dumb = -1;
            end
            res_dumb(n_dumb) = nanmin(m_dumb*padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
          end
          
          if any(isnan(y(end,:)) ~= isnan(res_dumb)),
            fprintf('SANITY CHECK: min, trans = flip, r = %d, off = %d, nan mismatch\n',r,off);
          else
            fprintf('SANITY CHECK: min, trans = flip, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
          end
          
        end

        
      end
    end
  end
end

%% window maximum

if useallfeatures || ismember('max',feature_types),
  for radiusi = 1:nradii,
    r = window_radii(radiusi);
    % no need to do minimum, maximum for r == 0 if we've done mean
    if r == 0 && (useallfeatures || ismember('mean',feature_types)),
      continue;
    end
    w = 2*r+1;
    % pad with infs for boundary conditions
    x_pad = [-inf(ntrans,r),x_trans,-inf(ntrans,r)];
    % maximum: use imdilate
    se = strel(ones(1,w));
    res = imdilate(x_pad,se);
    % offset
    windowis = find(windowi2radiusi == radiusi);
    for windowi = windowis',
      off = windows(windowi,2);
      % frame t, radius r, offset off:
      % [t-r+off, t+r+off]
      % so for r = 0, off = 1, we want [t+1,t+1]
      % which corresponds to res(t+r+off)
      % so we want to grab for 1+r+off through N+r+off
      res1 = padgrab(res,nan,1,ntrans,1+r+off,N+r+off);
      y(end+1,:) = res1(IDX_ORIG,:);
      feature_names{end+1} = {'stat','max','trans','none','radius',r,'offset',off};
      
      
      if SANITY_CHECK,
        
        res_dumb = nan(1,N);
        for n_dumb = 1:N,
          res_dumb(n_dumb) = nanmax(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
        end
        
        if any(isnan(y(end,:)) ~= isnan(res_dumb)),
          fprintf('SANITY CHECK: max, trans = none, r = %d, off = %d, nan mismatch\n',r,off);
        else
          fprintf('SANITY CHECK: max, trans = none, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
        end
        
      end
      
      if IDX_ABS > 0,
        y(end+1,:) = res1(IDX_ABS,:);
        feature_names{end+1} = {'stat','max','trans','abs','radius',r,'offset',off};
        
        if SANITY_CHECK,
          
          res_dumb = nan(1,N);
          for n_dumb = 1:N,
            res_dumb(n_dumb) = nanmax(abs(padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off)));
          end
          
          if any(isnan(y(end,:)) ~= isnan(res_dumb)),
            fprintf('SANITY CHECK: max, trans = abs, r = %d, off = %d, nan mismatch\n',r,off);
          else
            fprintf('SANITY CHECK: max, trans = abs, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
          end
          
        end
        
        
      end
      if IDX_FLIP > 0 && ~( (r == 0) && (off == 0) && (IDX_ABS > 0) ),
        y(end+1,:) = res1(IDX_ORIG,:);
        y(end,x_isneg) = res1(IDX_FLIP,x_isneg);
        feature_names{end+1} = {'stat','max','trans','flip','radius',r,'offset',off};
        
        
        if SANITY_CHECK,
          
          res_dumb = nan(1,N);
          for n_dumb = 1:N,
            m_dumb = 1;
            if sign(x(n_dumb)) < 0,
              m_dumb = -1;
            end
            res_dumb(n_dumb) = nanmax(m_dumb*padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off));
          end
          
          if any(isnan(y(end,:)) ~= isnan(res_dumb)),
            fprintf('SANITY CHECK: max, trans = flip, r = %d, off = %d, nan mismatch\n',r,off);
          else
            fprintf('SANITY CHECK: max, trans = flip, r = %d, off = %d, max error = %f\n',r,off,max(abs(y(end,:)-res_dumb)));
          end
          
        end
        
      end
    end
  end
end

%% histogram bin fraction

if useallfeatures || ismember('hist',feature_types),
  if isempty(hist_edges),
    warning('hist_edges is empty, but hist features are supposed to be computed. Skipping these features.');
  else
    nbins = numel(hist_edges)-1;
    % find bin for each value -- notice the transpose
    [~,bin] = histc(x_trans',hist_edges);
    for radiusi = 1:nradii,
      r = window_radii(radiusi);
      res = zeros(nbins,N+2*r,ntrans);
      for j = 1:ntrans,
        % initialize first value:
        % t contains [max(1,t-r),min(N,t+r)]
        % t = 1 will contain 1:r+1
        res(bin(1,j),1,j) = 1;
        for n = 1-r+1:N+r,
          res(:,n+r,j) = res(:,n+r-1,j);
          if n > r+1,
            res(bin(n-r-1,j),n+r,j) = res(bin(n-r-1,j),n+r,j) - 1;
          end
          if n <= N-r,
            res(bin(n+r,j),n+r,j) = res(bin(n+r,j),n+r,j) + 1;
          end
        end
      end
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
          feature_names{end+1} = {'stat','hist','trans','none','radius',r,'offset',off,'bin',bini,'lim',hist_edges(bini:bini+1)};
        end
        
        if SANITY_CHECK,
        
          res_real = y(end-nbins+1:end,:);
          res_dumb = nan(nbins,N);
          for n_dumb = 1:N,
            tmp = padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off);
            if all(isnan(tmp)),
              res_dumb(:,n_dumb) = nan;
            else
              tmp = histc(tmp,hist_edges);
              res_dumb(:,n_dumb) = tmp(1:end-1);
            end
          end       
          if any(isnan(res_real(:)) ~= isnan(res_dumb(:))),
            fprintf('SANITY CHECK: hist, trans = none, r = %d, off = %d, nan mismatch\n',r,off);
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
              if all(isnan(tmp)),
                res_dumb(:,n_dumb) = nan;
              else
                tmp = histc(tmp,hist_edges);
                res_dumb(:,n_dumb) = tmp(1:end-1);
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
            y(end,x_isneg) = res1(bini,x_isneg,IDX_FLIP);
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
              if all(isnan(tmp)),
                res_dumb(:,n_dumb) = nan;
              else
                tmp = histc(tmp,hist_edges);
                res_dumb(:,n_dumb) = tmp(1:end-1);
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
  end
end

%% percentiles

if useallfeatures || ismember('prctile',feature_types),
  if ~any(prctiles),
    warning('prctiles to compute is empty');
  else
    nprctiles = numel(prctiles);
    for radiusi = 1:nradii,
      r = window_radii(radiusi);
      % no need to do prctiles for r == 0
      if r == 0,
        continue;
      end

      % no clever way to do this
      res = nan(nprctiles,N+2*r);
      for n = 1-r:N+r,
        res(:,n+r) = prctile(x(max(1,n-r):min(N,n+r)),prctiles);
      end
      
      % offset
      windowis = find(windowi2radiusi == radiusi);
      for windowi = windowis',
        off = windows(windowi,2);
        % frame t, radius r, offset off:
        % [t-r+off, t+r+off]
        % so for r = 0, off = 1, we want [t+1,t+1]
        % which corresponds to res(t+r+off)
        % so we want to grab for 1+r+off through N+r+off
        res1 = padgrab(res,nan,1,nprctiles,1+r+off,N+r+off);
        y(end+1:end+nprctiles,:) = res1;
        for prctilei = 1:nprctiles,
          feature_names{end+1} = {'stat','prctile','trans','none','radius',r,'offset',off,'prctile',prctiles(prctilei)};
        end
        
        if SANITY_CHECK,
        
          res_real = y(end-nprctiles+1:end,:);
          res_dumb = nan(nprctiles,N);
          for n_dumb = 1:N,
            tmp = padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off);
            if all(isnan(tmp)),
              res_dumb(:,n_dumb) = nan;
            else
              res_dumb(:,n_dumb) = prctile(tmp,prctiles);
            end
          end       
          if any(isnan(res_real(:)) ~= isnan(res_dumb(:))),
            fprintf('SANITY CHECK: prctile, trans = none, r = %d, off = %d, nan mismatch\n',r,off);
          else
            fprintf('SANITY CHECK: prctile, trans = none, r = %d, off = %d, max error = %f\n',r,off,max(abs(res_real(:)-res_dumb(:))));
          end
        
        end
        
        
        if IDX_ABS > 0,
          y(end+1:end+nprctiles,:) = abs(res1);
          for prctilei = 1:nprctiles,
            feature_names{end+1} = {'stat','prctile','trans','abs','radius',r,'offset',off,'prctile',prctiles(prctilei)};
          end
          
          if SANITY_CHECK,
            
            res_real = y(end-nprctiles+1:end,:);
            res_dumb = nan(nprctiles,N);
            for n_dumb = 1:N,
              tmp = padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off);
              if all(isnan(tmp)),
                res_dumb(:,n_dumb) = nan;
              else
                res_dumb(:,n_dumb) = abs(prctile(tmp,prctiles));
              end
            end
            if any(isnan(res_real(:)) ~= isnan(res_dumb(:))),
              fprintf('SANITY CHECK: prctile, trans = abs, r = %d, off = %d, nan mismatch\n',r,off);
            else
              fprintf('SANITY CHECK: prctile, trans = abs, r = %d, off = %d, max error = %f\n',r,off,max(abs(res_real(:)-res_dumb(:))));
            end
            
          end
          
        end
        if IDX_FLIP > 0,
          for prctilei = 1:nprctiles,
            y(end+1,:) = res1(prctilei,:);
            y(end,x_isneg) = -res1(prctilei,x_isneg);
            feature_names{end+1} = {'stat','prctile','trans','flip','radius',r,'offset',off,'prctile',prctiles(prctilei)};
          end
          
          
          if SANITY_CHECK,
            
            res_real = y(end-nprctiles+1:end,:);
            res_dumb = nan(nprctiles,N);
            for n_dumb = 1:N,
              tmp = padgrab(x,nan,1,1,n_dumb-r+off,n_dumb+r+off);
              if all(isnan(tmp)),
                res_dumb(:,n_dumb) = nan;
              else
                tmp = prctile(tmp,prctiles);
                if x(n_dumb) < 0,
                  tmp = -tmp;
                end
                res_dumb(:,n_dumb) = tmp;
              end
            end
            if any(isnan(res_real(:)) ~= isnan(res_dumb(:))),
              fprintf('SANITY CHECK: prctile, trans = flip, r = %d, off = %d, nan mismatch\n',r,off);
            else
              fprintf('SANITY CHECK: prctile, trans = flip, r = %d, off = %d, max error = %f\n',r,off,max(abs(res_real(:)-res_dumb(:))));
            end
            
          end
          
          
        end
      end
      
      
    end
  end
end


%% change between end and start

if useallfeatures || ismember('change',feature_types),
  for change_r_i = 1:numel(change_window_radii),

    % take the mean across all windows of radius change_r
    change_r = change_window_radii(change_r_i);
    change_w = 2*change_r+1;

    % average filter
    fil = ones(1,change_w);
    % full: res_mean(t+change_r) corresponds to frame t
    res_mean = imfilter(x,fil,'full',0);

    % normalize
    res_mean(change_w:end-change_w+1) = res_mean(change_w:end-change_w+1) / change_w;
    % boundary conditions
    res_mean(1:change_w-1) = bsxfun(@rdivide,res_mean(1:change_w-1),1:change_w-1);
    res_mean(end-change_w+2:end) = bsxfun(@rdivide,res_mean(end-change_w+2:end),change_w-1:-1:1);
    
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
        y(end+1,:) = res1;
        feature_names{end+1} = {'stat','change','trans','none','radius',r,'offset',off,'change_window_radius',change_r};
      
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
        
        if IDX_ABS > 0,
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
        
        if IDX_FLIP > 0,
          res2 = res1;
          res2(x_isneg) = -res2(x_isneg);
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
end

%% standard deviation

if useallfeatures || ismember('std',feature_types),

  for radiusi = 1:nradii,
    r = window_radii(radiusi);
    w = 2*r+1;
    % average
    fil = ones(1,w);
    % full: res(t+r) corresponds to frame t
    % res_mean doesn't really need to be computed if mean has been computed
    res_mean = imfilter(x,fil,'full',0);
    res = imfilter(x.^2,fil,'full',0);
    % normalize
    res_mean(w:end-w+1) = res_mean(w:end-w+1) / w;
    res(w:end-w+1) = res(w:end-w+1) / w;
    % boundary conditions
    res_mean(1:w-1) = bsxfun(@rdivide,res_mean(1:w-1),1:w-1);
    res(1:w-1) = bsxfun(@rdivide,res(1:w-1),1:w-1);
    res_mean(end-w+2:end) = bsxfun(@rdivide,res_mean(end-w+2:end),w-1:-1:1);
    res(end-w+2:end) = bsxfun(@rdivide,res(end-w+2:end),w-1:-1:1);

    % combine to get standard deviation
    res = sqrt(res - res_mean.^2);
    
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
      y(end+1,:) = res1; 
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
end

%% harmonic

if useallfeatures || ismember('mean',feature_types) && num_harmonic > 0,
  
  for radiusi = 1:nradii,
    r = window_radii(radiusi);
    w = 2*r+1;
    
    for num_harmonic_curr = 1:num_harmonic,

      % can't split into more segments than frames
      if num_harmonic_curr >= w,
        continue;
      end
      
      fil = cos(linspace(0,pi*num_harmonic_curr,w))/w*(num_harmonic_curr+1);
      % res is nan from 1:r, N-r+1:N
      res = myconv(x,fil,nan,'corr','same');
      % boundary conditions: smallw used for frame
      for smallr = 1:r-1,
        smallw = 2*smallr+1;
        fil = cos(linspace(0,pi*num_harmonic_curr,smallw))/smallw*(num_harmonic_curr+1);
        % start
        res(smallr+1) = sum(fil.*x(1:smallw));
        % end
        res(end-smallr) = sum(fil.*x(end-smallw+1:end));
      end

      % all offsets for this radius
      windowis = find(windowi2radiusi == radiusi);
      for windowi = windowis',
        off = windows(windowi,2);
        res1 = padgrab(res,nan,1,1,1+off,N+off);
        y(end+1,:) = res1; 
        feature_names{end+1} = {'stat','harmonic','trans','none','radius',r,'offset',off,'num_harmonic',num_harmonic_curr};

        if SANITY_CHECK,
        
          res_dumb = nan(1,N);
          for n_dumb = 1:N,
            r_dumb = min([r,n_dumb+off-1,N-n_dumb-off]);
            if r_dumb < 1,
              continue;
            end
            w_dumb = 2*r_dumb+1;
            fil_dumb = cos(linspace(0,pi*num_harmonic_curr,w_dumb))/w_dumb*(num_harmonic_curr+1);
            res_dumb(n_dumb) = sum(fil_dumb.*x(n_dumb+off-r_dumb:n_dumb+off+r_dumb));
          end
        
          if any(isnan(y(end,:)) ~= isnan(res_dumb)),
            fprintf('SANITY CHECK: harmonic, trans = none, r = %d, off = %d, num_harmonic = %d, nan mismatch\n',r,off,num_harmonic_curr);
          else
            fprintf('SANITY CHECK: harmonic, trans = none, r = %d, off = %d, num_harmonic_curr = %d, max error = %f\n',r,off,num_harmonic_curr,max(abs(y(end,:)-res_dumb)));
          end
        
        end
      
        if IDX_ABS > 0,
          y(end+1,:) = abs(res1);
          feature_names{end+1} = {'stat','harmonic','trans','abs','radius',r,'offset',off,'num_harmonic',num_harmonic_curr};
          
          if SANITY_CHECK,
            
            res_dumb = nan(1,N);
            for n_dumb = 1:N,
              r_dumb = min([r,n_dumb+off-1,N-n_dumb-off]);
              if r_dumb < 1,
                continue;
              end
              w_dumb = 2*r_dumb+1;
              fil_dumb = cos(linspace(0,pi*num_harmonic_curr,w_dumb))/w_dumb*(num_harmonic_curr+1);
              res_dumb(n_dumb) = abs(sum(fil_dumb.*x(n_dumb+off-r_dumb:n_dumb+off+r_dumb)));
            end
            
            if any(isnan(y(end,:)) ~= isnan(res_dumb)),
              fprintf('SANITY CHECK: harmonic, trans = abs, r = %d, off = %d, num_harmonic = %d, nan mismatch\n',r,off,num_harmonic_curr);
            else
              fprintf('SANITY CHECK: harmonic, trans = abs, r = %d, off = %d, num_harmonic_curr = %d, max error = %f\n',r,off,num_harmonic_curr,max(abs(y(end,:)-res_dumb)));
            end
            
          end
        end
        
      end
    end
    
  end
end