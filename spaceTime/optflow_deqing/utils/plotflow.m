function [] = plotflow(f, kind)
%PLOTFLOW   Plot flow field
%   PLOTFLOW(F[, KIND]) plots an optical flow field F.  The optional
%   argument KIND specifies the kind of flow field plot:
%    - 'quiver', 'vector': needle-type plot (default)
%    - 'rgb': color plot (blue encodes U, green encodes V)
%    - 'hsv': color plot (hue encodes angle, value encodes velocity)
%    - 'bw': grayscale plot (U on the left, V on the right) 
%    - 'bwscale': also print out the flow scaling
%    - 'mag': grayscale plot of flow magnitude
%    - 'magscale': also print out the flow scaling
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
  
  if (nargin < 2)
    kind = 'vector';
  end
  
  switch (kind)
   case {'quiver', 'vector'}
    s = size(f);
%     step = max(s / 120);
%     step = max(s / 60);
    step = max(s / 40);
%     step = max(s / 20);
    
    [X, Y] = meshgrid(1:step:s(2), s(1):-step:1);
    u = interp2(f(:, :, 1), X, Y);
    v = interp2(f(:, :, 2), X, Y);
    
%     quiver(X, -Y, u, -v, 0.7);
    quiver(X, -Y, u, -v, 1, 'k', 'LineWidth', 1);
    axis image;
    
   case 'rgb'
    b = f(:, :, 1);
    b = b - min(b(:));
    b = b / max(b(:));
    
    g = f(:, :, 2);
    g = g - min(g(:));
    g = g / max(g(:));

    r = zeros(size(b));
    
    [ignore, rad] = cart2pol(f(:, :, 1), f(:, :, 2));
    
    nanidx  = isnan(f(:, :, 1)) & isnan(f(:, :, 2));
    zeroidx = (rad < 0.1);
    
    r(nanidx) = 0;
    g(nanidx) = 0;
    b(nanidx) = 0;
    r(zeroidx) = 1;
    g(zeroidx) = 1;
    b(zeroidx) = 1;
    
    im = cat(3, r, g, b);
    image(im);
    
   case 'hsv'
    [theta, rho] = cart2pol(f(:, :, 1), f(:, :, 2));
    
    theta = (theta + pi) / (2*pi);
    rho   = rho / max(rho(:));
    
    im = cat(3, theta, ones(size(theta)), rho);
    image(hsv2rgb(im));
    
   case {'bw', 'bwscale'}
    f1 = f(:, :, 1);
    f2 = f(:, :, 2);
    im = [f1, f2];
    
    m = max(abs(im(:)));
    imagesc(im, [-m m]);
    colormap gray(256);
    axis image
    axis off
    if (strcmp(kind, 'bwscale'))
      title(['[' num2str(min(im(:))) '; ' num2str(max(im(:))) ']'])
      colorbar
    end
    
   case {'mag', 'magscale'}
    m = sqrt(sum(f.^2, 3));
    imagesc(m);
    colormap gray(256);
    axis image
    axis off
    if (strcmp(kind, 'magscale'))
      title(['[' num2str(max(m(:))) ']'])
      colorbar
    end
    
   otherwise
    error('Invalid plot type');
    
  end