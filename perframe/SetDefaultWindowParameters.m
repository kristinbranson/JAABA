function [default_windows,default_window_offsets,...
  default_window_radii,default_min_window_radius,...
  default_max_window_radius,default_nwindow_radii] = ...
  SetDefaultWindowParameters()

% specify window using cross product of all radii and all offsets
default_windows = []; 

default_window_offsets = [-1,0,1];

% specify radii using minimum, maximum radii, and number of window radii
default_window_radii = [];
default_min_window_radius = 0;
default_max_window_radius = 20;
default_nwindow_radii = 5;