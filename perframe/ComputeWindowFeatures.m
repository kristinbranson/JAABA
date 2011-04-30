function y = ComputeWindowFeatures(x,varargin)

% defaults parameters for windows

% specify window using cross product of all radii and all offsets
windows = []; 

% specify radii using minimum, maximum radii, and number of window radii
window_radii = [];
min_window_radius = 0;
max_window_radius = 20;
nwindow_radii = 5;

window_offsets = [-1,0,1];

[...
  windows,...
  window_radii,window_offsets,...
  min_window_radius,max_window_radius,nwindow_radii,...
  ] = myparse(varargin,...
  'windows',windows,...
  'window_radii',window_radii,'window_offsets',window_offsets,...
  'min_window_radius',min_window_radius,'max_window_radius',max_window_radius,'nwindow_radii',nwindow_radii);

% select windows
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
  
  % compute start and end dts
  dt0s = -all_radii + all_offsets;
  dt1s = all_radii + all_offsets;
  windows = [dt0s,dt1s];
  % make sure these are unique -- rounding to nearest frame might make them
  % not unique
  windows = unique(windows,'rows');
  
  % window i for frame t is t+t0(i):t+t1(i)
  
end

if isempty(windows),
  error('windows is empty.');
end
if size(windows,2) ~= 2,
  error('windows must be nwindows x 2.');
end

