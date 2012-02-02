function h_group = impolygon(h_parent, x, y)

n = length(x);
if length(y) ~= n,
  error('x and y must be the same length');
end;

if ~ishandle(h_parent)
  error('Images:impolygon:invalidHandle', ...
    'HPARENT must be a valid graphics handle.');
end
if strcmp(get(h_parent, 'type'), 'axes')
  h_axes = h_parent;
else
  h_axes = ancestor(h_parent, 'axes');
  if isempty(h_axes)
    error('Images:impolygon:noAxesAncestor', ...
      'HPARENT must be a descendent of an axes object.');
  end
end

try
  h_group = hggroup('ButtonDownFcn', @startDrag, 'Parent', h_parent, ...
    'HitTest', 'on');
catch
  error('Images:impolygon:failureToParent', ...
    'HPARENT must be able to have an hggroup object as a child.');
end

draw_api = wingedPolygon(h_group,n);


function draw_API = wingedPolygon(h_group,n)
%wingedRect Creates renderer for winged rectangle symbol.
%   DRAW_API = wingedRect(H_GROUP) creates a DRAW_API for use in association
%   with IMRECT that draws rectangles with wings that show only if the rectangle
%   is very small. DRAW_API is a structure of function handles that are used by
%   IMRECT to draw the rectangle and update its properties.
%
%       DRAW_API.setColor
%       DRAW_API.updateView
%       DRAW_API.getBoundingBox

%   Copyright 2005-2006 The MathWorks, Inc.
%   $Revision $  $Date: 2006/06/15 20:10:16 $

% initialize variables needing function scope
bounding_box = [];
h_axes = iptancestor(h_group,'axes');

% The line objects should have a width of one screen pixel.
% line_width = ceil(getPointsPerScreenPixel);
line_width = getPointsPerScreenPixel();
h_line = zeros(1,n);
for i = 1:n,
  h_line = line('Color', 'w', ...
    'LineStyle', '-', ...
    'LineWidth', line_width, ...
    'HitTest', 'off', ...
    'Parent', h_group);
end;

h_patch = patch('FaceColor', 'none', 'EdgeColor', 'none', ...
    'HitTest', 'off', ...
    'Parent', h_group);

draw_API.setColor         = @setColor;
draw_API.updateView       = @updateView;
draw_API.getBoundingBox   = @getBoundingBox;

    %----------------------------
    function pos = getBoundingBox
        pos = bounding_box;
    end

    %----------------------------
    function updateView(position)

        if ~ishandle(h_group)
            return;
        end

        [dx_per_screen_pixel, dy_per_screen_pixel] = getAxesScale(h_axes);

        min_decorated_rect_size = 30;
        x_left = position(1) / dx_per_screen_pixel;
        x_right = (position(1) + position(3)) / dx_per_screen_pixel;
        x_wing_size = max(ceil((min_decorated_rect_size - ...
            (x_right - x_left)) / 2), 0);

        y_bottom = position(2) / dy_per_screen_pixel;
        y_top = (position(2) + position(4)) / dy_per_screen_pixel;
        y_wing_size = max(ceil((min_decorated_rect_size - ...
            (y_top - y_bottom)) / 2), 0);

        x1 = x_left - x_wing_size;
        x2 = x_left;
        x3 = (x_left + x_right) / 2;
        x4 = x_right;
        x5 = x_right + x_wing_size;

        y1 = y_bottom - y_wing_size;
        y2 = y_bottom;
        y3 = (y_bottom + y_top) / 2;
        y4 = y_top;
        y5 = y_top + y_wing_size;

        % (x,y) is a polygon that strokes the middle line.  Here it is in
        % screen pixel units.
        x = [x1 x2 x2 x3 x3 x3 x4 x4 x5 x4 x4 x3 x3 x3 x2 x2 x1];
        y = [y3 y3 y2 y2 y1 y2 y2 y3 y3 y3 y4 y4 y5 y4 y4 y3 y3];

        % Convert the (x,y) polygon back to data units.
        [x,y] = pixel2DataUnits(h_axes,x,y);

        xx1 = x1 - 1;
        xx2 = x2 - 1;
        xx3 = x3 - 1;
        xx4 = x3 + 1;
        xx5 = x4 + 1;
        xx6 = x5 + 1;

        yy1 = y1 - 1;
        yy2 = y2 - 1;
        yy3 = y3 - 1;
        yy4 = y3 + 1;
        yy5 = y4 + 1;
        yy6 = y5 + 1;

        % (xx,yy) is a polygon that strokes the outer line.  Here it is in
        % screen pixel units.
        xx = [xx1 xx2 xx2 xx3 xx3 xx4 xx4 xx5 xx5 xx6 xx6 xx5 xx5 ...
            xx4 xx4 xx3 xx3 xx2 xx2 xx1 xx1];
        yy = [yy3 yy3 yy2 yy2 yy1 yy1 yy2 yy2 yy3 yy3 yy4 yy4 yy5 ...
            yy5 yy6 yy6 yy5 yy5 yy4 yy4 yy3];

        % Convert the (xx,yy) polygon back to data units.
        [xx,yy] = pixel2DataUnits(h_axes,xx,yy);

        % Set the outer position to include the entire extent of the drawn
        % rectangle, including decorations.
        bounding_box = findBoundingBox(xx,yy);

        setXYDataIfChanged([h_top_line,h_bottom_line,h_patch], x, y);
    end

    %------------------
    function setColor(c)
        if ishandle(h_top_line)
            set(h_top_line, 'Color', c);
        end
    end
end

%----------------------------------
function setXYDataIfChanged(h, x, y)
% Set XData and YData of HG object h to x and y if they are different.  h
% must be a valid HG handle to an object having XData and YData properties.
% No validation is performed.

if ~isequal(get(h, 'XData'), x) || ~isequal(get(h, 'YData'), y)
    set(h, 'XData', x, 'YData', y);
end
end

