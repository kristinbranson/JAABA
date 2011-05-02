function [y,feature_names] = ComputeWindowFeatures(x,varargin)

x = x(:)';
N = numel(x);
y = nan(0,N);
feature_names = {};

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

[...
  windows,...
  window_radii,window_offsets,...
  min_window_radius,max_window_radius,nwindow_radii,...
  feature_types,trans_types
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii,...
  'feature_types',feature_types,'trans_types',trans_types);

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
nwindows = size(windows,1);

%% compute transformations
x_trans = nan(1,N);
IDX_ABS = 0;
IDX_FLIP = 0;
if usealltrans || ismember('abs',trans_types),
  x_trans(end+1,:) = abs(x);
  IDX_ABS = size(x_trans,2);
end
if usealltrans || ismember('signflip',trans_types),
  x_trans(end+1,:) = -x;
  IDX_FLIP = size(x_trans,2);
end

%% window mean

if useallfeatures || ismember('mean',feature_types),
  for windowi = 1:nwindows,
    
  end
end