function O = imwarp(I, u, v, nopad)
%IMWARP   Warp image with flow field
%   O = IMWARP(I, U, V[, NOPAD]) warps the image I with optical flow field
%   U and V.  The flow field is to be given in the coordinate system of
%   O, i.e. the operation warps I toward O.  If the optional argument
%   NOPAD is given and true, the undefined pixels (source pixel outside
%   the boundary) are given by the nearest boundary pixels, otherwise
%   they are set to NaN.
%
%   Author:  Stefan Roth, Department of Computer Science, TU Darmstadt
%   Contact: sroth@cs.tu-darmstadt.de
%   $Date: 2007-03-27 14:09:11 -0400 (Tue, 27 Mar 2007) $
%   $Revision: 252 $

% Copyright 2004-2007, Brown University, Providence, RI. USA
% Copyright 2007-2010 TU Darmstadt, Darmstadt, Germany.
% 
%                          All Rights Reserved
% 
% All commercial use of this software, whether direct or indirect, is
% strictly prohibited including, without limitation, incorporation into in
% a commercial product, use in a commercial service, or production of other
% artifacts for commercial purposes.     
%
% Permission to use, copy, modify, and distribute this software and its
% documentation for research purposes is hereby granted without fee,
% provided that the above copyright notice appears in all copies and that
% both that copyright notice and this permission notice appear in
% supporting documentation, and that the name of the author and Brown
% University not be used in advertising or publicity pertaining to
% distribution of the software without specific, written prior permission.        
%
% For commercial uses contact the Technology Venture Office of Brown University
% 
% THE AUTHOR AND BROWN UNIVERSITY DISCLAIM ALL WARRANTIES WITH REGARD TO
% THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY AND
% FITNESS FOR ANY PARTICULAR PURPOSE.  IN NO EVENT SHALL THE AUTHOR OR
% BROWN UNIVERSITY BE LIABLE FOR ANY SPECIAL, INDIRECT OR CONSEQUENTIAL
% DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
% PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
% ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
% THIS SOFTWARE.        
  
  % Image size
  sx = size(I, 2); 
  sy = size(I, 1); 

  % Image size w/ padding
  spx = sx + 2;
  spy = sy + 2;
  
  
  if (nargin > 3 & nopad)  
    % Warped image coordinates
    [X, Y] = meshgrid(1:sx, 1:sy);
    XI = reshape(X + u, 1, sx * sy);
    YI = reshape(Y + v, 1, sx * sy);
    
    % Bound coordinates to valid region
    XI = max(1, min(sx - 1E-6, XI));
    YI = max(1, min(sy - 1E-6, YI));
    
    % Perform linear interpolation (faster than interp2)
    fXI = floor(XI);
    cXI = ceil(XI);
    fYI = floor(YI);
    cYI = ceil(YI);
    
    alpha_x = XI - fXI;
    alpha_y = YI - fYI;
    
    O = (1 - alpha_x) .* (1 - alpha_y) .* I(fYI + sy * (fXI - 1)) + ...
        alpha_x .* (1 - alpha_y) .* I(fYI + sy * (cXI - 1)) + ...
        (1 - alpha_x) .* alpha_y .* I(cYI + sy * (fXI - 1)) + ...
        alpha_x .* alpha_y .* I(cYI + sy * (cXI - 1));
  
  else
    % Pad image with NaNs
    Z = [NaN(1, sx+2); NaN(sy, 1), I, NaN(sy, 1); NaN(1, sx+2)];
    
    % Warped image coordinates in padded image
    [X, Y] = meshgrid(2:sx+1, 2:sy+1);
    XI = reshape(X + u, 1, sx * sy);
    YI = reshape(Y + v, 1, sx * sy);
    
    % Bound coordinates to valid region
    XI = max(1, min(spx - 1E-6, XI));
    YI = max(1, min(spy - 1E-6, YI));
    
    % Perform linear interpolation (faster than interp2)
    fXI = floor(XI);
    cXI = ceil(XI);
    fYI = floor(YI);
    cYI = ceil(YI);
    
    alpha_x = XI - fXI;
    alpha_y = YI - fYI;
    
    O = (1 - alpha_x) .* (1 - alpha_y) .* Z(fYI + spy * (fXI - 1)) + ...
        alpha_x .* (1 - alpha_y) .* Z(fYI + spy * (cXI - 1)) + ...
        (1 - alpha_x) .* alpha_y .* Z(cYI + spy * (fXI - 1)) + ...
        alpha_x .* alpha_y .* Z(cYI + spy * (cXI - 1));
  end
  
  O = reshape(O, sy, sx);
