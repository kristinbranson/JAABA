% [y,feature_names,cache] = ComputePrctileWindowFeatures(x,...)
% 
% Computes the percentile window features y for the input one-dimensional
% per-frame time series data x. y(i,t) will correspond to percentile
% window feature i and frame t, and is the (transformation of) the min of
% all the per-frame data within the window defined by i at t. 
% 
% For window feature i, let p_i be the percentile to compute, r_i be the
% window radius, and off_i be the window offset. If the transformation type
% for feature i is 'none', the window feature y(i,t) is the p_i th
% percentile of the per-frame data within the window from t-r_i+off_i
% through t+r_i+off_i. If the transformation type for feature i is
% 'abs', then we compute the absolute value of the percentile of the
% per-frame data within the window. If the transformation type for feature
% i is 'flip', then we compute the sign of x(t) times the percentile of the
% per-frame data within the window.
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
% {'stat','prctile','trans','abs','radius',1,'offset',1,'prctile',50}
% cache: if DOCACHE is set to true, then computations useful in other
% window feature calculations may be saved in the output cache. 
%
% Optional inputs:
%
% 'prctiles': The 1 x NPRCTILES percentiles to compute. default value: []. 
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
% min computations are useful when computing the diff_neighbor_min.
% default value: true.
%
% 'cache': input cached computations that might be useful in this
% computation. 

function [y,feature_names,cache] = ComputePrctileWindowFeatures(x,varargin)

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
%feature_types = {};

% percentiles to compute (other than min, max)
prctiles = [];

%% parse parameters

[...
  windows,...
  window_radii,window_offsets,...
  min_window_radius,max_window_radius,nwindow_radii,...
  trans_types,...
  SANITY_CHECK,...
  DOCACHE,...
  cache,...
  prctiles...
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'trans_types',trans_types,...
  'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE,...
  'cache',cache,...
  'prctiles',prctiles); %#ok<ASGLU>
%   'feature_types',feature_types,...

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

if ~any(prctiles),
  warning('prctiles to compute is empty');
  return;
end

nprctiles = numel(prctiles);
for radiusi = 1:nradii,

  r = window_radii(radiusi);

  % no need to do prctiles for r == 0
  %if r == 0 && any(ismember({'min','mean','max'},feature_types)),
  %  continue;
  %end
  
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
    
    if ismember('none',trans_types),
      
      y(end+1:end+nprctiles,:) = res1;
      for prctilei = 1:nprctiles,
        feature_names{end+1} = {'stat','prctile','trans','none','radius',r,'offset',off,'prctile',prctiles(prctilei)}; %#ok<*AGROW>
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
      
    end
    
    if ismember('abs',trans_types),
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
    if ismember('flip',trans_types),
      for prctilei = 1:nprctiles,
        y(end+1,:) = res1(prctilei,:);
        y(end,x<0) = -res1(prctilei,x<0);
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
