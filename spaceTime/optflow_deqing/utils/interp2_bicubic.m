function [ZI, ZXI, ZYI, C, alpha_x, alpha_y, fXI, fYI] = interp2_bicubic(Z, XI, YI, Dxfilter)
% function varargout = interp2_bicubic(Z, XI, YI)
%   Author:  Stefan Roth, Department of Computer Science, TU Darmstadt
%   Contact: sroth@cs.tu-darmstadt.de
%   $Date$
%   $Revision$

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
  % Implementation according to Numerical Recipes
  
  % modified by dqsun
  if nargin < 4
      Dxfilter = [-0.5 0 0.5];
  end;
  Dyfilter = Dxfilter';
  Dxyfilter = conv2(Dxfilter, Dyfilter, 'full');
  
  input_size = size(XI);
  
  % Reshape input coordinates into a vector
  XI = reshape(XI, 1, prod(input_size));
  YI = reshape(YI, 1, prod(input_size));
  
  % Bound coordinates to valid region
  sx = size(Z, 2);
  sy = size(Z, 1);
  
% $$$   vXI = max(1, min(sx - 1E-6, XI));
% $$$   vYI = max(1, min(sy - 1E-6, YI));
  
  % Neighbor coordinates
% $$$   fXI = floor(vXI);
% $$$   cXI = fXI + 1;
% $$$   fYI = floor(vYI);
% $$$   cYI = fYI + 1;
  fXI = floor(XI);
  cXI = fXI + 1;
  fYI = floor(YI);
  cYI = fYI + 1;
 
%   indx = (fXI>sx) | (fXI<1) | (cXI>sx) | (cXI<1) | (fYI>sy) | (fYI<1) | (cYI>sy) | (cYI<1);
  indx = (fXI<1) | (cXI>sx) | (fYI<1) | (cYI>sy);
  
  fXI = max(1, min(sx, fXI));
  cXI = max(1, min(sx, cXI));
  fYI = max(1, min(sy, fYI));
  cYI = max(1, min(sy, cYI));

  
  % Image at 4 neighbors
  Z00 = Z(fYI + sy * (fXI - 1));
  Z01 = Z(cYI + sy * (fXI - 1));
  Z10 = Z(fYI + sy * (cXI - 1));
  Z11 = Z(cYI + sy * (cXI - 1));
  
  % x-derivative at 4 neighbors
%   DX = imfilter(Z, [-0.5, 0, 0.5], 'symmetric', 'corr');
  DX = imfilter(Z, Dxfilter, 'symmetric', 'corr'); % modified by dqsun
  DX00 = DX(fYI + sy * (fXI - 1));
  DX01 = DX(cYI + sy * (fXI - 1));
  DX10 = DX(fYI + sy * (cXI - 1));
  DX11 = DX(cYI + sy * (cXI - 1));

  % y-derivative at 4 neighbors
%   DY = imfilter(Z, [-0.5, 0, 0.5]', 'symmetric', 'corr'); 
  DY = imfilter(Z, Dyfilter, 'symmetric', 'corr'); % modified by dqsun
  DY00 = DY(fYI + sy * (fXI - 1));
  DY01 = DY(cYI + sy * (fXI - 1));
  DY10 = DY(fYI + sy * (cXI - 1));
  DY11 = DY(cYI + sy * (cXI - 1));

  % xy-derivative at 4 neighbors
%   DXY = imfilter(Z, 0.25 * [1, 0, -1; 0 0 0; -1 0 1], 'symmetric', 'corr');
  DXY = imfilter(Z, Dxyfilter, 'symmetric', 'corr'); % modified by dqsun
  DXY00 = DXY(fYI + sy * (fXI - 1));
  DXY01 = DXY(cYI + sy * (fXI - 1));
  DXY10 = DXY(fYI + sy * (cXI - 1));
  DXY11 = DXY(cYI + sy * (cXI - 1));

  
  W = [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
        0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0;
       -3,  0,  0,  3,  0,  0,  0,  0, -2,  0,  0, -1,  0,  0,  0,  0;  
        2,  0,  0, -2,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,  0;
        0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0;
        0,  0,  0,  0, -3,  0,  0,  3,  0,  0,  0,  0, -2,  0,  0, -1;
        0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0,  0,  1,  0,  0,  1;
       -3,  3,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
        0,  0,  0,  0,  0,  0,  0,  0, -3,  3,  0,  0, -2, -1,  0,  0;
        9, -9,  9, -9,  6,  3, -3, -6,  6, -6, -3,  3,  4,  2,  1,  2;
       -6,  6, -6,  6, -4, -2,  2,  4, -3,  3,  3, -3, -2, -1, -1, -2;
        2, -2,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0;
        0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  1,  1,  0,  0;
       -6,  6, -6,  6, -3, -3,  3,  3, -4,  4,  2, -2, -2, -2, -1, -1;
        4, -4,  4, -4,  2,  2, -2, -2,  2, -2, -2,  2,  1,  1,  1,  1];

  
  V = [Z00(:)'; Z10(:)'; Z11(:)'; Z01(:)'; ...
       DX00(:)'; DX10(:)'; DX11(:)'; DX01(:)'; ...
       DY00(:)'; DY10(:)'; DY11(:)'; DY01(:)'; ...
       DXY00(:)'; DXY10(:)'; DXY11(:)'; DXY01(:)'];
  
  C = W * V;
  
  alpha_x = reshape(XI - fXI, input_size);
  alpha_y = reshape(YI - fYI, input_size);
  
  % Clip out-of-boundary pixels to boundary
  % Modified by Deqing Sun (7-10-2008)
  alpha_x(indx) = 0;
  alpha_y(indx) = 0;

  
  fXI = reshape(fXI, input_size);
  fYI = reshape(fYI, input_size);
  
  % Interpolation

  ZI = zeros(input_size);
  ZXI = zeros(input_size);
  ZYI = zeros(input_size);
  
  idx = 1;
  for i = 0:3
    for j = 0:3
      ZI = ZI + reshape(C(idx, :), input_size) .* alpha_x.^i .* ...
           alpha_y.^j;
      if (i > 0) & (nargout>=2)
        ZXI = ZXI + i * reshape(C(idx, :), input_size) .* alpha_x.^(i-1) .* ...
              alpha_y.^j;
      end
      if (j > 0) & (nargout>=3)
        ZYI = ZYI + j * reshape(C(idx, :), input_size) .* alpha_x.^i .* ...
              alpha_y.^(j - 1);
      end
      idx = idx + 1;
    end
  end  

ZI(indx) = NaN;
