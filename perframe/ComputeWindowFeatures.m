% [y,feature_names] = ComputeWindowFeatures(x,...)
% 
% Wrapper function that computes all the window features y for the input
% one-dimensional time series x. Window features that can be computed by
% this function include:
% * mean: (transformation of) the average of the per-frame data within the
% window. 
% * min: minimum of the (transformation of) the per-frame data within the
% window. 
% * max: maximum of the (transformation of) the per-frame data within the
% window. 
% * hist: fraction of the (transformation of) the per-frame data within the
% window within each bin specified by the input edges. 
% * prctile: (transformation of) the percentiles of the per-frame data
% within the window. 
% * change: (transformation of) the difference between the mean of the
% sub-window of per-frame data centered at the end of the window and the
% mean of the sub-window of per-frame data centered at the start of the
% window. 
% * std: standard deviation of the per-frame data within the window. 
% * harmonic: (transformation of) the correlation between the per-frame
% data in the window and cosine functions for n / 2 periods, for n =
% 1,...,num_harmonic. 
% * diff_neighbor_mean: (transformation of) the difference between the
% per-frame data at the current frame and the mean of the per-frame data in
% the window around it. 
% diff_neighbor_min: (transformation of) the difference between the
% per-frame data at the current frame and the minimum of the per-frame data
% in the window around it. 
% diff_neighbor_max: (transformation of) the difference between the
% per-frame data at the current frame and the maximum of the per-frame data
% in the window around it. 
% zscore_neighbors: (transformation of) the per-frame data at the current
% frame z-scored by the statistics of the per-frame data in the in the
% window around it.
% 
% Input: 
%
% x: 1 x NFRAMES array of per-frame data. 
% 
% Output: 
%
% y: NWINDOWFEATURES x NFRAMES matrix of window data, where y(i,t)
% corresponds to window feature i and frame t. 
% feature_names: 1 x NWINDOWFEATURES cell in which feature_names{i}
% describes the ith window feature computed. feature_names{i} is itself a
% cell that can be interpreted as pairs of a string description followed by
% a value, e.g. 
% {'stat','zscore_neighbors','trans','abs','radius',1,'offset',1}
% 
% Optional inputs:
% 
% DEFAULT WINDOW LOCATIONS:
% These window locations are used if window locations are not specified on
% a per-feature type basis. Default default values set by
% SetDefaultWindowParameters.
% Inputs interpreted by SetWindowParameters. 
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
% 'feature_types': list of the types of window features to compute, chosen
% from the list of window feature types listed above. 'feature_types' can also be 
% the string 'all', in which case all features are computed. 
%
% 'trans_types': default types of transformations to apply for each feature
% type, if not otherwise specified. Options include 'abs','flip', and
% 'none'. 'abs' corresponds to the absolute value, flip corresponds to
% flipping the sign of the feature if the per-frame feature at the frame t
% is negative, and 'none' corresponds to no transformation. default value:
% {'none','abs','flip'}. 
%
% PER-FEATURE TYPE PARAMETERS:
% parameters specific to each feature-type, passed to the
% Compute*WindowFeatures functions. Default for all of these is {}. 
% 'mean_params': cell array of parameters to be passed to
% ComputeMeanWindowFeatures. 
% 'min_params': cell array of parameters to be passed to
% ComputeMinWindowFeatures. 
% 'max_params': cell array of parameters to be passed to
% ComputeMaxWindowFeatures. 
% 'hist_params': cell array of parameters to be passed to
% ComputeHistWindowFeatures. 
% 'prctile_params': cell array of parameters to be passed to
% ComputePrctileWindowFeatures. 
% 'change_params': cell array of parameters to be passed to
% ComputeChangeWindowFeatures. 
% 'std_params': cell array of parameters to be passed to
% ComputeStdWindowFeatures. 
% 'harmonic_params': cell array of parameters to be passed to
% ComputeHarmonicWindowFeatures. 
% 'diff_neighbor_mean_params': cell array of parameters to be passed to
% ComputeDiffNeighborMeanWindowFeatures. 
% 'diff_neighbor_min_params': cell array of parameters to be passed to
% ComputeDiffNeighborMinWindowFeatures. 
% 'diff_neighbor_max_params': cell array of parameters to be passed to
% ComputeDiffNeighborMaxWindowFeatures. 
% 'zscore_neighbors_params': cell array of parameters to be passed to
% ComputeZScoreNeighborsWindowFeatures. 
% 
% 'sanitycheck': whether to compute all the features a second time in the
% obvious way to make sure that the optimized computations are correct.
% default value: false. 
%
% 'docache': whether to cache computations that might be useful, e.g. the
% mean computations are useful when computing the diff_neighbor_mean.
% default value: true.
%
% 't0': first frame of per-frame data to compute window features for. if
% empty, then window features are computed for all frames. default value:
% []. 
% 't1': last frame of per-frame data to compute window features for. if
% empty, then window features are computed for all frames. default value:
% []. 
%
% TODO: Check for redundant features
%
function [y,feature_names] = ComputeWindowFeatures(x,varargin)

%% initialize

x = x(:)';
N = numel(x);
feature_names = {};
feature_types_computed = {};

%% default parameters

% defaults parameters for windows, if not specified per window feature
[default_windows,default_window_offsets,...
  default_window_radii,default_min_window_radius,...
  default_max_window_radius,default_nwindow_radii] = ...
  SetDefaultWindowParameters();

% use all features, transformation types by default
feature_types = 'all';
default_trans_types = 'all';

% per-feature parameters
mean_params = {};
min_params = {};
max_params = {};
hist_params = {};
prctile_params = {};
change_params = {};
std_params = {};
harmonic_params = {};
diff_neighbor_mean_params = {};
diff_neighbor_min_params = {};
diff_neighbor_max_params = {};
zscore_neighbors_params = {};

% for debugging purposes
SANITY_CHECK = false;

% whether to cache results
DOCACHE = true;

%% parse parameters

%   change_window_radii,...
%   num_harmonic,...

[...
  default_windows,...
  default_window_radii,default_window_offsets,...
  default_min_window_radius,default_max_window_radius,default_nwindow_radii,...
  feature_types,...
  default_trans_types,...
  mean_params,...
  min_params,...
  max_params,...
  hist_params,...
  prctile_params,...
  change_params,...
  std_params,...
  harmonic_params,...
  diff_neighbor_mean_params,...
  diff_neighbor_min_params,...
  diff_neighbor_max_params,...
  zscore_neighbors_params,...
  SANITY_CHECK,...
  DOCACHE,...
  t0,t1...
  ] = myparse(varargin,...
  'windows',default_windows,...
  'window_radii',default_window_radii,'window_offsets',default_window_offsets,...
  'min_window_radius',default_min_window_radius,'max_window_radius',default_max_window_radius,'nwindow_radii',default_nwindow_radii,...
  'feature_types',feature_types,...
  'trans_types',default_trans_types,...
  'mean_params',mean_params,...
  'min_params',min_params,...
  'max_params',max_params,...
  'hist_params',hist_params,...
  'prctile_params',prctile_params,...
  'change_params',change_params,...
  'std_params',std_params,...
  'harmonic_params',harmonic_params,...
  'diff_neighbor_mean_params',diff_neighbor_mean_params,...
  'diff_neighbor_min_params',diff_neighbor_min_params,...
  'diff_neighbor_max_params',diff_neighbor_max_params,...
  'zscore_neighbors_params',zscore_neighbors_params,...
  'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE,...
  't0',[],'t1',[]);

%% whether we've specified to use all features

useallfeatures = ischar(feature_types) && strcmpi(feature_types,'all');

%% select default windows from various ways of specifying windows

default_windows = SetWindowParameters(...
  default_windows,default_window_offsets,...
  default_window_radii,default_min_window_radius,...
  default_max_window_radius,default_nwindow_radii);

%% set all parameters, putting in default windows, trans as necessary

mean_params = SetDefaultParams(mean_params,default_windows,default_trans_types);
min_params = SetDefaultParams(min_params,default_windows,default_trans_types);
max_params = SetDefaultParams(max_params,default_windows,default_trans_types);
hist_params = SetDefaultParams(hist_params,default_windows,default_trans_types);
prctile_params = SetDefaultParams(prctile_params,default_windows,default_trans_types);
change_params = SetDefaultParams(change_params,default_windows,default_trans_types);
std_params = SetDefaultParams(std_params,default_windows,default_trans_types);
harmonic_params = SetDefaultParams(harmonic_params,default_windows,default_trans_types);
diff_neighbor_mean_params = SetDefaultParams(diff_neighbor_mean_params,default_windows,default_trans_types);
diff_neighbor_min_params = SetDefaultParams(diff_neighbor_min_params,default_windows,default_trans_types);
diff_neighbor_max_params = SetDefaultParams(diff_neighbor_max_params,default_windows,default_trans_types);
zscore_neighbors_params = SetDefaultParams(zscore_neighbors_params,default_windows,default_trans_types);

common_params = {'sanitycheck',SANITY_CHECK,...
  'docache',DOCACHE};

%% we'll need to compute means, stds mins, maxs for multiple features
cache = InitializeCache();

%% if t0 or t1 are input, then only compute in a boundary around t0 through t1
% compute the maximum dependency radius and only compute up to there. 
% there's more efficient ways to do this, but this is the simplest
% modification to the code

MAXDEPENDENCYRADIUS = ComputeMaxDependencyRadius({'windows',default_windows},...
  mean_params,min_params,max_params,hist_params,prctile_params,change_params,std_params,...
  harmonic_params,diff_neighbor_mean_params,diff_neighbor_min_params,...
  diff_neighbor_max_params,zscore_neighbors_params);
  
if isempty(t0),
  t0 = 1;
  OFF0 = 1;
else
  OFF0 = max(1,t0-MAXDEPENDENCYRADIUS);
end

if isempty(t1),
  t1 = N;
  OFF1 = N;
else
  OFF1 = min(N,t1+MAXDEPENDENCYRADIUS);
end
x = x(OFF0:OFF1);
y = zeros(0,OFF1-OFF0+1);

%% window mean

if useallfeatures || ismember('mean',feature_types),
  [y_new,feature_names_new,cache] = ...
    ComputeMeanWindowFeatures(x,'cache',cache,...
    common_params{:},mean_params{:});
  %'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'mean';
end

%% window minimum

if useallfeatures || ismember('min',feature_types),
  [y_new,feature_names_new,cache] = ...
    ComputeMinWindowFeatures(x,'cache',cache,...
    common_params{:},min_params{:});
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'min';
end

%% window maximum

if useallfeatures || ismember('max',feature_types),
  [y_new,feature_names_new,cache] = ...
    ComputeMaxWindowFeatures(x,'cache',cache,...
    common_params{:},max_params{:});
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'max';
end

%% histogram bin fraction

if useallfeatures || ismember('hist',feature_types),
  [y_new,feature_names_new,cache] = ...
    ComputeHistWindowFeatures(x,'cache',cache,...
    common_params{:},hist_params{:});
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'hist';
end

%% percentiles

if useallfeatures || ismember('prctile',feature_types),
  [y_new,feature_names_new,cache] = ...
    ComputePrctileWindowFeatures(x,'cache',cache,...
    common_params{:},prctile_params{:});
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'prctile';
end

%% change between end and start

if useallfeatures || ismember('change',feature_types),
  
  [y_new,feature_names_new,cache] = ...
    ComputeChangeWindowFeatures(x,'cache',cache,...
    common_params{:},change_params{:});
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'change';

end

%% standard deviation

if useallfeatures || ismember('std',feature_types),

  [y_new,feature_names_new,cache] = ...
    ComputeStdWindowFeatures(x,'cache',cache,...
    common_params{:},std_params{:});
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'std';

end

%% harmonic

if useallfeatures || ismember('harmonic',feature_types),
  
  [y_new,feature_names_new,cache] = ...
    ComputeHarmonicWindowFeatures(x,'cache',cache,...
    common_params{:},harmonic_params{:});
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'harmonic';

end

%% difference from neighbor mean

if useallfeatures || ismember('diff_neighbor_mean',feature_types),

  [y_new,feature_names_new,cache] = ...
    ComputeDiffNeighborMeanWindowFeatures(x,'cache',cache,...
    common_params{:},diff_neighbor_mean_params{:});
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'diff_neighbor_mean';

end

%% difference from neighbor min

if useallfeatures || ismember('diff_neighbor_min',feature_types),

  [y_new,feature_names_new,cache] = ...
    ComputeDiffNeighborMinWindowFeatures(x,'cache',cache,...
    common_params{:},diff_neighbor_min_params{:});
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'diff_neighbor_min';
  
end

%% difference from neighbor max

if useallfeatures || ismember('diff_neighbor_max',feature_types),

  [y_new,feature_names_new,cache] = ...
    ComputeDiffNeighborMaxWindowFeatures(x,'cache',cache,...
    common_params{:},diff_neighbor_max_params{:});
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'diff_neighbor_max';

end

%% z-score wrt neighbors

if useallfeatures || ismember('zscore_neighbors',feature_types),

  [y_new,feature_names_new,cache] = ...
    ComputeZScoreNeighborsWindowFeatures(x,'cache',cache,...
    common_params{:},zscore_neighbors_params{:}); %#ok<NASGU>
  % 'feature_types',feature_types_computed,...
  y(end+1:end+size(y_new,1),:) = y_new;
  feature_names(end+1:end+numel(feature_names_new)) = feature_names_new;
  feature_types_computed{end+1} = 'zscore_neighbors'; %#ok<NASGU>

end

%% get the part of y that we actually want to return
pad0 = t0-OFF0;
pad1 = OFF1-t1;
y = y(:,1+pad0:end-pad1);