function y = ximresize(x, nsz, method)
%XIMRESIZE   Fast version of IMRESIZE that also supports adjoint operators
%   XIMRESIZE(X, SZ[, METHOD]) resizes image X to the new size SZ.  The
%   optional METHOD argument specifies the method to be used:
%    - 'nearest': resize using nearest neighbor interpolation  
%    - 'nearest_adj': resize using the adjoint operator of nearest
%                     neighbor interpolation
%    - 'bilinear': resize using bilinear interpolation  
%    - 'bilinear_adj': resize using the adjoint operator of bilinear
%                      interpolation  
%  
%   See also IMRESIZE.
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

  
  if (ndims(x) > 2)
    error('>= 3D not implemented');
  end
  
  if (nargin < 3)
    method = 'nearest';
  end
  
  % Resize factor
  sz = size(x);
  if (length(nsz) == 1)
    nsz = round(nsz * sz);
  end
  
  switch (method)
   case 'nearest'
    
    rows = 1 + (sz(1)/nsz(1)) * (0:nsz(1)-1)';
    cols = 1 + (sz(2)/nsz(2)) * (0:nsz(2)-1);
  
    rows = floor(rows); rows(rows < 1) = 1; rows(rows > sz(1)) = sz(1);
    cols = floor(cols); cols(cols < 1) = 1; cols(cols > sz(2)) = sz(2);
    
    ind = xrepmat(rows, 1, length(cols)) + ...
          xrepmat(sz(1) * (cols-1), length(rows), 1);
    
    y = x(ind);
    
    
   case 'nearest_adj'
    
    rows = 1 + (nsz(1)/sz(1)) * (0:sz(1)-1)';
    cols = 1 + (nsz(2)/sz(2)) * (0:sz(2)-1);
  
    rows = floor(rows); rows(rows < 1) = 1; rows(rows > nsz(1)) = nsz(1);
    cols = floor(cols); cols(cols < 1) = 1; cols(cols > nsz(2)) = nsz(2);
    
    ind = xrepmat(rows, 1, length(cols)) + ...
          xrepmat(nsz(1) * (cols-1), length(rows), 1);

    A = sparse(ind(:), 1:numel(ind), ones(numel(ind), 1), prod(nsz), prod(sz));
    y = reshape(A * x(:), nsz);

    
   case 'bilinear'
    rows = 1 + ((sz(1) - 1)/(nsz(1) - 1)) * (0:nsz(1)-1)';
    cols = 1 + ((sz(2) - 1)/(nsz(2) - 1)) * (0:nsz(2)-1);

    brows = floor(rows); brows(brows < 1) = 1; brows(brows > sz(1)) = sz(1);
    bcols = floor(cols); bcols(bcols < 1) = 1; bcols(bcols > sz(2)) = sz(2);
    trows = ceil(rows); trows(trows < 1) = 1; trows(trows > sz(1)) = sz(1);
    tcols = ceil(cols); tcols(tcols < 1) = 1; tcols(tcols > sz(2)) = sz(2);
    
    ralpha = rows - brows;
    calpha = cols - bcols;

    % bottom, bottom
    ind = xrepmat(brows, 1, length(cols)) + ...
          xrepmat(sz(1) * (bcols-1), length(rows), 1);
    alpha = xrepmat(1 - ralpha, 1, length(cols)) .* ...
            xrepmat(1 - calpha, length(rows), 1);
    
    y = alpha .* x(ind);

    % top, bottom
    ind = xrepmat(trows, 1, length(cols)) + ...
          xrepmat(sz(1) * (bcols-1), length(rows), 1);
    alpha = xrepmat(ralpha, 1, length(cols)) .* ...
            xrepmat(1 - calpha, length(rows), 1);
    
    y = y + alpha .* x(ind);

    % bottom, top
    ind = xrepmat(brows, 1, length(cols)) + ...
          xrepmat(sz(1) * (tcols-1), length(rows), 1);
    alpha = xrepmat(1 - ralpha, 1, length(cols)) .* ...
            xrepmat(calpha, length(rows), 1);
    
    y = y + alpha .* x(ind);

    % top, top
    ind = xrepmat(trows, 1, length(cols)) + ...
          xrepmat(sz(1) * (tcols-1), length(rows), 1);
    alpha = xrepmat(ralpha, 1, length(cols)) .* ...
            xrepmat(calpha, length(rows), 1);
    
    y = y + alpha .* x(ind);

    
   case 'bilinear_adj' % Adjoint bilinear resize operator
    
    rows = 1 + ((nsz(1) - 1)/(sz(1) - 1)) * (0:sz(1)-1)';
    cols = 1 + ((nsz(2) - 1)/(sz(2) - 1)) * (0:sz(2)-1);

    brows = floor(rows); brows(brows < 1) = 1; brows(brows > nsz(1)) = nsz(1);
    bcols = floor(cols); bcols(bcols < 1) = 1; bcols(bcols > nsz(2)) = nsz(2);
    trows = ceil(rows); trows(trows < 1) = 1; trows(trows > nsz(1)) = nsz(1);
    tcols = ceil(cols); tcols(tcols < 1) = 1; tcols(tcols > nsz(2)) = nsz(2);
    
    ralpha = rows - brows;
    calpha = cols - bcols;

    % bottom, bottom
    ind = xrepmat(brows, 1, length(cols)) + ...
          xrepmat(nsz(1) * (bcols-1), length(rows), 1);
    alpha = xrepmat(1 - ralpha, 1, length(cols)) .* ...
            xrepmat(1 - calpha, length(rows), 1);
    
    A = sparse(ind(:), 1:numel(ind), ones(numel(ind), 1));
    tmp = alpha .* x;
    y = reshape(A * tmp(:), nsz);
    

    % top, bottom
    ind = xrepmat(trows, 1, length(cols)) + ...
          xrepmat(nsz(1) * (bcols-1), length(rows), 1);
    alpha = xrepmat(ralpha, 1, length(cols)) .* ...
            xrepmat(1 - calpha, length(rows), 1);
    
    A = sparse(ind(:), 1:numel(ind), ones(numel(ind), 1));
    tmp = alpha .* x;
    y = y + reshape(A * tmp(:), nsz);


    % bottom, top
    ind = xrepmat(brows, 1, length(cols)) + ...
          xrepmat(nsz(1) * (tcols-1), length(rows), 1);
    alpha = xrepmat(1 - ralpha, 1, length(cols)) .* ...
            xrepmat(calpha, length(rows), 1);
    
    A = sparse(ind(:), 1:numel(ind), ones(numel(ind), 1));
    tmp = alpha .* x;
    y = y + reshape(A * tmp(:), nsz);


    % top, top
    ind = xrepmat(trows, 1, length(cols)) + ...
          xrepmat(nsz(1) * (tcols-1), length(rows), 1);
    alpha = xrepmat(ralpha, 1, length(cols)) .* ...
            xrepmat(calpha, length(rows), 1);
    
    A = sparse(ind(:), 1:numel(ind), ones(numel(ind), 1));
    tmp = alpha .* x;
    y = y + reshape(A * tmp(:), nsz);

    
   otherwise
    error('Method not implemented');
  end