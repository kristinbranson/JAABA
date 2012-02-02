function truesize(varargin)
%TRUESIZE Adjust display size of image.
%   TRUESIZE(FIG,[MROWS NCOLS]) adjusts the display size of an
%   image. FIG is a figure containing a single image or a single 
%   image with a colorbar. [MROWS NCOLS] is a 1-by-2 vector that
%   specifies the requested screen area (in pixels) that the
%   image should occupy.
%
%   TRUESIZE(FIG) uses the image height and width for 
%   [MROWS NCOLS]. This results in the display having one screen
%   pixel for each image pixel.
%
%   If you omit the figure argument, TRUESIZE works on the
%   current figure.
%
%   Example
%   --------
%       % Fit image to figure window.
%       imshow(checkerboard,'InitialMagnification','fit') 
%
%       % Resize image and figure to show image at its
%       % 80-by-80 pixel size.
%       truesize
%
%   See also IMSHOW.

%   Copyright 1993-2007 The MathWorks, Inc.
%   $Revision: 1.1.8.7 $  $Date: 2007/11/09 20:22:33 $

[axHandle, imHandle, imSize, ...
    oneImageAndOneAxes, oneImageAndOneColorbar] = ParseInputs(varargin{:});

figHandle = ancestor(axHandle, 'figure');
if strcmp(get(figHandle, 'WindowStyle'), 'docked')
    wid = 'Images:truesize:dockedFigure';
    warning(wid, 'TRUESIZE cannot adjust the size of a docked figure.');
    return
end

% Find out if image fills the figure.
axesPos = getpixelposition(axHandle);
figPos = getpixelposition(figHandle);
isBorderTight = isequal(axesPos(3:4), figPos(3:4));

if oneImageAndOneAxes
    if isempty(imSize)
        % Figure contains one image with or without a colorbar.
        % Syntax for initSize(imageHandle, screenPerImageImagePixel,
        % isBorderTight)
        initSize(imHandle, 1, isBorderTight);
    else
        Resize1(axHandle, imHandle, imSize);
    end
elseif oneImageAndOneColorbar && isempty(imSize)
    % Figure contains image and a colorbar and did not specify imSize
    % Syntax for initSize(imageHandle, screenPerImageImagePixel,
    % isBorderTight)
    initSize(imHandle, 1, isBorderTight);
else
    Resize2(axHandle, imHandle, imSize);
end


%--------------------------------------------
% Subfunction Resize1
%--------------------------------------------
function Resize1(axHandle, imHandle, imSize)
% Resize figure containing a single axes
% object with a single image.

if (isempty(imSize))
    % How big is the image?
    imageWidth = size(get(imHandle, 'CData'), 2);
    imageHeight = size(get(imHandle, 'CData'), 1);
else
    imageWidth = imSize(2);
    imageHeight = imSize(1);
end

if (imageWidth * imageHeight == 0)
    % Don't try to handle the degenerate case.
    return;
end

axUnits = get(axHandle, 'Units');
set(axHandle, 'Units', 'pixels');
axPos = get(axHandle, 'Position');

figHandle = ancestor(axHandle, 'figure');
figUnits = get(figHandle, 'Units');
rootUnits = get(0, 'Units');
set(figHandle, 'Units', 'pixels');
set(0, 'Units', 'pixels');

figLeftBorder = 10;  % assume left figure decorations are 10 pixels
figRightBorder = 10;
figBottomBorder = 10;
figTopBorder = 100;

minFigWidth = 128; % don't try to display a figure smaller than this.
minFigHeight = 128;

% What are the gutter sizes?
figPos = get(figHandle, 'Position');
gutterLeft = max(axPos(1) - 1, 0);

% What are the screen dimensions
screenSize = get(0, 'ScreenSize');
screenWidth = screenSize(3);
screenHeight = screenSize(4);
if ((screenWidth <= 1) || (screenHeight <= 1))
    screenWidth = Inf;
    screenHeight = Inf;
end

scale = 100;
done = 0;
defAxesPos = get(0,'DefaultAxesPosition');
nonzeroGutters = (gutterLeft > 0);
while (~done)
    if (nonzeroGutters)
        gutterWidth = round((1 - defAxesPos(3)) * imageWidth / defAxesPos(3));
        gutterHeight = round((1 - defAxesPos(4)) * imageHeight / defAxesPos(4));
        newFigWidth = imageWidth + gutterWidth;
        newFigHeight = imageHeight + gutterHeight;
    else
        newFigWidth = imageWidth;
        newFigHeight = imageHeight;
    end
    if (((newFigWidth + figLeftBorder + figRightBorder) > screenWidth) || ...
                ((newFigHeight + figBottomBorder + figTopBorder) > screenHeight))
        scale = 3 * scale / 4;
        imageWidth = round(imageWidth * scale / 100);
        imageHeight = round(imageHeight * scale / 100);
    else
        done = 1;
    end
end

newFigWidth = max(newFigWidth, minFigWidth);
newFigHeight = max(newFigHeight, minFigHeight);

figPos(1) = max(1, figPos(1) - floor((newFigWidth - figPos(3))/2));
figPos(2) = max(1, figPos(2) - floor((newFigHeight - figPos(4))/2));
figPos(3) = newFigWidth;
figPos(4) = newFigHeight;

% Translate figure position if necessary
deltaX = (screenSize(3) - figRightBorder) - (figPos(1) + figPos(3));
if (deltaX < 0)
    figPos(1) = figPos(1) + deltaX;
end
deltaY = (screenSize(4) - figTopBorder) - (figPos(2) + figPos(4));
if (deltaY < 0)
    figPos(2) = figPos(2) + deltaY;
end

% Figure out where to place the axes object in the
% resized figure.  Make sure axes width and height are at least 1.
gutterWidth = figPos(3) - imageWidth;
gutterHeight = figPos(4) - imageHeight;
gutterLeft = floor(gutterWidth/2);
gutterBottom = floor(gutterHeight/2);

axPos(1) = gutterLeft + 1;
axPos(2) = gutterBottom + 1;
axPos(3) = max(imageWidth,1);
axPos(4) = max(imageHeight,1);

set(figHandle, 'Position', figPos)
set(axHandle, 'Position', axPos);

% Restore the units
drawnow;  % necessary to work around HG bug   -SLE
set(figHandle, 'Units', figUnits);
set(axHandle, 'Units', axUnits);
set(0, 'Units', rootUnits);

% Set the nextplot property of the figure so that the
% axes object gets deleted and replaced for the next plot.
% That way, the new plot gets drawn in the default position.
set(figHandle, 'NextPlot', 'replacechildren');

% Warn if the display is not truesize.
if (scale < 100)
    wid = sprintf('Images:%s:imageTooBigForScreen',mfilename);
    warning(wid, ...
        'Image is too big to fit on screen; displaying at %d%% scale.', ...
        floor(scale));
end

%--------------------------------------------
% Subfunction Resize2
%--------------------------------------------
% Resize figure containing multiple axes or
% other objects.  Basically we're going to
% compute a global figure scaling factor
% that will bring the target image into
% truesize mode.  This works reasonably well
% for subplot-type figures as long as all
% the images have the same size.  This is
% basically the guts of truesize.m from IPT
% version 1.
function Resize2(axHandle, imHandle, imSize)

if (isempty(imSize))
    imSize = size(get(imHandle, 'CData'));
end

if (prod(imSize) == 0)
    % Don't try to handle the degenerate case.
    return;
end

axUnits = get(axHandle, 'Units');
set(axHandle, 'Units', 'pixels');
axPosition = get(axHandle, 'Position');
% Do we need to do anything?
if norm(axPosition(3:4) - imSize([2 1])) < sqrt(eps)
    set(axHandle, 'Units', axUnits)
    return;
end

figHandle = ancestor(axHandle, 'figure');
figUnits = get(figHandle, 'Units');
rootUnits = get(0,'Units');
set(axHandle, 'Units', 'normalized');
axPosition = get(axHandle, 'Position');
set([figHandle 0], 'Units', 'pixels');
figPosition = get(figHandle, 'Position');
screenSize = get(0,'ScreenSize');

% What should the new figure size be?
dx = ceil(imSize(2)/axPosition(3) - figPosition(3));
dy = ceil(imSize(1)/axPosition(4) - figPosition(4));
newFigWidth = figPosition(3) + dx; 
newFigHeight = figPosition(4) + dy;

% Is the new figure size too big or too small?
figLeftBorder = 10;
figRightBorder = 10;
figBottomBorder = 10;
figTopBorder = 100;
minFigWidth = 128;
minFigHeight = 128;

if ((newFigWidth + figLeftBorder + figRightBorder) > screenSize(3))
    scaleX = (screenSize(3) - figLeftBorder - figRightBorder) / newFigWidth;
else
    scaleX = 1;
end

if ((newFigHeight + figBottomBorder + figTopBorder) > screenSize(4))
    scaleY = (screenSize(4) - figBottomBorder - figTopBorder) / newFigHeight;
else
    scaleY = 1;
end

if (newFigWidth < minFigWidth)
    scaleX = minFigWidth / newFigWidth;
end
if (newFigHeight < minFigHeight)
    scaleY = minFigHeight / newFigHeight;
end


% make sure scaling factor does not make newFigWidth and newFigHeight less
% than one.

if (min(scaleX, scaleY) < 1)
    % Yes, the new figure is too big for the screen.
    scale = min(scaleX, scaleY); 
    newFigWidth = floor(max(newFigWidth*scale,1)); 
    newFigHeight = floor(max(newFigHeight*scale,1)); 

elseif (max(scaleX, scaleY) > 1)
    % Yes, the new figure is too small.
    scale = max(scaleX, scaleY);
    newFigWidth = floor(max(newFigWidth*scale,1)); 
    newFigHeight = floor(max(newFigHeight*scale,1));
else
    scale = 1;
end

figPosition(3) = newFigWidth; 
figPosition(4) = newFigHeight;

% Translate figure position if necessary
deltaX = (screenSize(3) - figRightBorder) - (figPosition(1) + figPosition(3));
if (deltaX < 0)
    figPosition(1) = figPosition(1) + deltaX;
end
deltaY = (screenSize(4) - figTopBorder) - (figPosition(2) + figPosition(4));
if (deltaY < 0)
    figPosition(2) = figPosition(2) + deltaY;
end

% Change axes position to get exactly one pixel per image pixel.
% That is, as long as scale = 1.
dx = scale*imSize(2)/figPosition(3) - axPosition(3);
dy = scale*imSize(1)/figPosition(4) - axPosition(4);
axPosition = axPosition + [-dx/2 -dy/2 dx dy];
axPosition(3:4) = max(axPosition(3:4),eps);

% OK, make the changes
set(axHandle, 'Position', axPosition);
set(figHandle, 'Position', figPosition);

% Restore the original units
set(axHandle, 'Units', axUnits);
set(figHandle, 'Units', figUnits);
set(0, 'Units', rootUnits);

% Warn if the display is not truesize

if (scale < 1)
    wid = sprintf('Images:%s:imageTooBigForScreen',mfilename);
    warning(wid, ...
        'Image is too big to fit on screen; displaying at %d%% scale.', ...
        round(100*scale));
elseif (scale > 1)
    wid = sprintf('Images:%s:imageTooSmallToScale',mfilename);        
    warning(wid, ...
        'Image is too small for truesize figure scaling; \ndisplaying at %d%% scale.', ...
        round(100*scale));
end

%--------------------------------------------
% Subfunction ParseInputs
%--------------------------------------------
function [axHandle,imHandle,imSize,oneImage1Axes,oneImage1Colorbar] = ...
        ParseInputs(varargin)

imSize = [];
  
if nargin >= 2
    %TRUESIZE(FIG,[M N],...)
    imSize = varargin{2};
    if (~isequal(size(imSize), [1 2]))
        eid = sprintf('Images:%s:reqsizeMustBe1by2',mfilename);    
        error(eid, 'REQSIZE must be a 1-by-2 vector.');
    end
    figHandle = varargin{1};
    iptcheckhandle(figHandle, {'figure'}, mfilename, 'FIG', 1);
    axHandle = get(figHandle, 'CurrentAxes');
elseif nargin == 1
    %TRUESIZE([M N])
    %TRUESIZE(FIG)
    if isequal(size(varargin{1}), [1 2])
        figHandle = get(0, 'CurrentFigure');
        imSize = varargin{1};
    else
        figHandle = varargin{1};
        iptcheckhandle(figHandle, {'figure'}, mfilename, 'FIG', 1);
    end
    axHandle = get(figHandle, 'CurrentAxes');
else
    %TRUESIZE
    figHandle = get(0, 'CurrentFigure');
    axHandle = get(figHandle, 'CurrentAxes');
end

if (isempty(axHandle))
    eid = sprintf('Images:%s:currentFigureMissingAxes',mfilename);
    error(eid, 'Current figure has no axes.');
end

% Find all the images and texturemapped surfaces
% in the current figure.  These are the candidates.
h = [findobj(figHandle, 'Type', 'image'); ...
    findobj(figHandle, 'Type', 'surface', 'FaceColor', 'texturemap')];

% If there's a colorbar, ignore it.
colorbarHandle = findobj(figHandle, 'type', 'image', ...
        'Tag', 'TMW_COLORBAR');
if ~isempty(colorbarHandle)
      h(ismember(h,colorbarHandle)) = [];
end

if isempty(h)
    eid = sprintf('Images:%s:noImagesOrSurfaces',mfilename);    
    error(eid, 'No images or texturemapped surfaces in the figure.');
end

% Start with the first object on the list as the
% initial candidate.  If it's not in the current
% axes, look for another one that is.

imHandle = h(1);
if (ancestor(imHandle,'axes') ~= axHandle)
    for k = 2:length(h)
        if (ancestor(h(k),'axes') == axHandle)
            imHandle = h(k);
            break;
        end
    end
end

% Determine if there are more than one axes and or a colorbar in 
% the figure. You need this information to determine appropriate
% truesize method.
axesKids = findobj(figHandle, 'type', 'axes');
noUicontrol = isempty( ...
    findobj(figHandle, 'Type', 'uicontrol', 'Visible', 'on'));
oneImage = length(h) == 1;
numAxes = length(axesKids);
oneImage1Axes = oneImage && numAxes == 1 && ...
    isequal(ancestor(imHandle, 'axes'), axesKids) && noUicontrol;
oneImage1Colorbar = oneImage && length(colorbarHandle) == 1 && ...
    numAxes == 2 && noUicontrol;