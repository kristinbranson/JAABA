%PRINT2ARRAY  Exports a figure to an image array
%
% Examples:
%   A = print2array
%   A = print2array(figure_handle)
%   A = print2array(figure_handle, resolution)
%   A = print2array(figure_handle, resolution, renderer)
%
% This function outputs a bitmap image of the given figure, at the desired
% resolution.
%
% If renderer is '-painters' then ghostcript needs to be installed. This
% can be downloaded from: http://www.ghostscript.com
%
% IN:
%   figure_handle - The handle of the figure to be exported. Default: gcf.
%   resolution - Resolution of the output, as a factor of screen
%                resolution. Default: 1.
%   renderer - string containing the renderer paramater to be passed to
%              print. Default: '-opengl'.
%
% OUT:
%   A - MxNx3 uint8 image of the figure.

% Copyright (C) Oliver Woodford 2008-2009

function A = print2array(fig, res, renderer)
% Generate default input arguments, if needed
if nargin < 2
    res = 1;
    if nargin < 1
        fig = gcf;
    end
end
% Warn if output is large
old_mode = get(fig, 'Units');
set(fig, 'Units', 'pixels');
px = get(fig, 'Position');
set(fig, 'Units', old_mode);
px = prod(px(3:4)*res)/1e6;
if px > 30
    % 30M pixels or larger!
    warning('MATLAB:LargeImage', 'print2array generating a %.1fM pixel image. This could be slow and might also cause memory problems.', px);
end
% Generate temporary file name
tmp_nam = [tempname '.tif'];
if nargin > 2 && strcmp(renderer, '-painters')
    % Print to eps file
    tmp_eps = [tempname '.eps'];
    print2eps(tmp_eps, fig, renderer, '-loose');
    try
        % Export to tiff using ghostscript
        ghostscript(['-dEPSCrop -q -dNOPAUSE -dBATCH -r' num2str(ceil(get(0, 'ScreenPixelsPerInch')*res)) ' -sDEVICE=tiff24nc -sOutputFile="' tmp_nam '" "' tmp_eps '"']);
    catch
        % Delete the intermediate file
        delete(tmp_eps);
        rethrow(lasterror);
    end
    % Delete the intermediate file
    delete(tmp_eps);
    % Read in the generated bitmap
    A = imread(tmp_nam);
    % Retrieve the background colour
    bcol = get(fig, 'Color');
    % Set border pixels to the correct colour
    if ~isequal(bcol, 'none') && ~isequal(bcol, [1 1 1])
        for l = 1:size(A, 2)
            if ~all(reshape(A(:,l,:) == 255, [], 1))
                break;
            end
        end
        for r = size(A, 2):-1:l
            if ~all(reshape(A(:,r,:) == 255, [], 1))
                break;
            end
        end
        for t = 1:size(A, 1)
            if ~all(reshape(A(t,:,:) == 255, [], 1))
                break;
            end
        end
        for b = size(A, 1):-1:t
            if ~all(reshape(A(b,:,:) == 255, [], 1))
                break;
            end
        end
        bcol = round(bcol * 255);
        for c = 1:size(A, 3)
            A(:,[1:l-1, r+1:end],c) = bcol(c);
            A([1:t-1, b+1:end],:,c) = bcol(c);
        end
    end
else
    if nargin < 3
        renderer = '-opengl';
    end
    % Set paper size
    old_mode = get(fig, 'PaperPositionMode');
    set(fig, 'PaperPositionMode', 'auto');
    % Print to tiff file
    print(fig, renderer, ['-r' num2str(ceil(get(0, 'ScreenPixelsPerInch')*res))], '-dtiff', tmp_nam);
    % Reset paper size
    set(fig, 'PaperPositionMode', old_mode);
    % Read in the printed file
    A = imread(tmp_nam);
end
% Delete the temporary file
delete(tmp_nam);
return
