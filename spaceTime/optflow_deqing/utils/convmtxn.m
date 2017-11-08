function M = convmtxn(F, sz, valid)
%CONVMTXN   N-D convolution matrix
%   M = CONVMTXN(F, SZ[, VALID]) returns the convolution matrix for the
%   matrix F.  SZ gives the size of the array that the convolution should
%   be applied to.  The returned matrix M is sparse.
%   If X is of size SZ, then reshape(M * X(:), SZ + size(F) - 1) is the
%   same as convn(X, F).
%   The optional parameter VALID controls which rows of M (corresponding
%   to the pixels of X) should be set to zero to suppress the convolution
%   result.  VALID must have size SZ+size(F)-1.
%  
%   See also CONVN, CONVMTX2, MAKE_CONVN_MTX.
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

  
  ndims  = length(sz);
  blksz  = prod(size(F));
  nblks  = prod(sz);
  nelems = blksz * nblks;

  
  % Build index array for all possible image positions
  tmp = zeros(size(F) + sz - 1);
  for d = 1:ndims
    sub{d} = 1:sz(d);
  end
  tmp(sub{:}) = 1;  
  imgpos = find(tmp(:));

  % Build index array for all possible filter positions
  tmp = zeros(size(F) + sz - 1);
  for d = 1:ndims
    sub{d} = 1:size(F, d);
  end
  tmp(sub{:}) = 1;  
  fltpos = find(tmp(:));
  
  
  % The loop code below is replaced with the vectorized version below  
% $$$   rows = zeros(nelems, 1);
% $$$   cols = zeros(nelems, 1);
% $$$   vals = zeros(nelems, 1);
% $$$   
% $$$   for i = 1:prod(sz)
% $$$     % For every possible image position (cols) insert filter values
% $$$     % into the appropriate output position (rows).
% $$$   
% $$$     j = 1 + (i-1) * blksz;
% $$$     k = i * blksz;
% $$$     rows(j:k) = imgpos(i) + fltpos - 1;
% $$$     cols(j:k) = i;
% $$$     vals(j:k) = F(:);
% $$$   end

  % Vectorized version of loop code above.
  rows = reshape(repmat(imgpos', blksz, 1), nelems, 1) + ...
         repmat(fltpos - 1, nblks, 1);
  cols = reshape(repmat(1:nblks, blksz, 1), nelems, 1);
  vals = repmat(F(:), nblks, 1);
%   rows = reshape(xrepmat(imgpos', blksz, 1), nelems, 1) + ...
%          xrepmat(fltpos - 1, nblks, 1);
%   cols = reshape(xrepmat(1:nblks, blksz, 1), nelems, 1);
%   vals = xrepmat(F(:), nblks, 1);

  % Pick out valid rows
  if (nargin > 2)
    valid_idx = valid(rows);

    rows = rows(valid_idx);
    cols = cols(valid_idx);
    vals = vals(valid_idx);
  end
  
  % Build sparse output array
  M = sparse(rows, cols, vals, prod(size(F) + sz - 1), nblks);
