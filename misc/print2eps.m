%PRINT2EPS  Prints figures to eps with improved line styles
%
% Examples:
%   print2eps filename
%   print2eps(filename, fig_handle)
%   print2eps(filename, fig_handle, options)
%
% This function saves a figure as an eps file, and improves the line style,
% making dashed lines more like those on screen and giving grid lines their
% own dotted style.
%
%IN:
%   filename - string containing the name (optionally including full or
%              relative path) of the file the figure is to be saved as. A
%              ".eps" extension is added if not there already. If a path is
%              not specified, the figure is saved in the current directory.
%   fig_handle - The handle of the figure to be saved. Default: gcf.
%   options - Additional parameter strings to be passed to print.

% Copyright (C) Oliver Woodford 2008-2009

% The idea of editing the EPS file to change line styles comes from Jiro
% Doke's FIXPSLINESTYLE (fex id: 17928)
% The idea of changing dash length with line width came from comments on
% fex id: 5743, but the implementation is mine :)

function print2eps(name, fig, varargin)
options = {'-depsc2'};
if nargin < 2
    fig = gcf;
elseif nargin > 2
    options = [options varargin];
end
% Construct the filename
if numel(name) < 5 || ~strcmpi(name(end-3:end), '.eps')
    name = [name '.eps']; % Add the missing extension
end
% Set paper size
old_mode = get(fig, 'PaperPositionMode');
set(fig, 'PaperPositionMode', 'auto');
% Print to eps file
print(fig, options{:}, name);
% Reset paper size
set(fig, 'PaperPositionMode', old_mode);
% Fix the line styles
try
    fix_lines(name);
catch
    warning('fix_lines failed. This is usually because the figure contains a large number of patch objects. Consider exporting to a bitmap format in this case.');
end
return
