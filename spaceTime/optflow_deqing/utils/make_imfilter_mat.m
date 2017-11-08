function M = make_imfilter_mat(F, sz, bndry, shape)
%MAKE_IMFILTER_MAT   Flexible N-D image filtering matrix
%   M = MAKE_IMFILTER_MAT(F, SZ[, BNDRY[, SHAPE]]) returns the
%   image filtering matrix for the matrix F.  SZ gives the size of the
%   array that the filtering (convolution) should be applied to.  The
%   returned matrix M is sparse. 
%   If X is of size SZ and SHAPE is 'full', then reshape(M * X(:), SZ +
%   size(F) - 1) is the same as imfilter(X, F, BNDRY, 'conv').
%
%   The optional parameter BNDRY controls the boundary handling:
%     * 'replicate':  Input array values outside the bounds of the array
%                     are assumed to equal the nearest array border
%                     value.
%     * 'circular':   Input array values outside the bounds of the array
%                     wrap around to the opposite bound.
%     * '0':          (default) Input array values outside of the bounds
%                     of the array are assumed to be zero.  
%
%   The optional parameter SHAPE controls which parts of the filtering
%   result to return in M:
%     * 'full':      returns the full N-D convolution, i.e. the 
%                    behavior is identical to imfilter(X, F, bndry,
%                    'full', 'conv'). 
%     * 'same':      (default) returns the central part of the
%                    convolution that is the same size as X, i.e. the
%                    behavior is identical to imfilter(X, F, bndry,
%                    'same', 'conv').
%     * 'valid':     returns only the part of the result that can be
%                    computed without assuming zero-padded arrays,
%                    i.e. the behavior is identical to imfilter(X, F,
%                    bndry, 'valid', 'conv').
%
%   See also IMFILTER, IMFILTERMTX, MAKE_CONVN_MAT.
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


  % Default to 'same' matrix
  if (nargin < 4)
    shape = 'same';
  end
  
  % Default to zero boundaries
  if (nargin < 3)
    bndry = '0';
  end
  
  ndims = length(sz);
  
  % Border sizes for 'same'
  Fsize_2 = round((size(F) - 1) / 2);

  switch(shape)
    case 'same'
      % Mark valid and invalid pixels (i.e. the ones within and outside
      % of the part to be returned)
      valid = true(sz+size(F)-1);
      
      for d = 1:ndims
        for e = 1:ndims
          sub{e} = ':';
        end
        
        sub{d} = 1:Fsize_2(d);      
        valid(sub{:}) = false;
        sub{d} = size(valid, d)-Fsize_2(d)+1:size(valid, d);      
        valid(sub{:}) = false;
      end
      
      % Image filter matrix with appropriate boundary handling
      M = imfiltermtx(F, sz, bndry);
      % Suppress rows of M outside of the valid area
      M = M(valid, :);
      
    case 'valid'
      % Mark valid and invalid pixels (i.e. the ones within and outside
      % of the part to be returned)
      valid = true(sz+size(F)-1);
      
      for d = 1:ndims
        for e = 1:ndims
          sub{e} = ':';
        end
        
        sub{d} = 1:2*Fsize_2(d);      
        valid(sub{:}) = false;
        sub{d} = size(valid, d)-2*Fsize_2(d)+1:size(valid, d);      
        valid(sub{:}) = false;
      end
      
      % Image filter matrix with appropriate boundary handling
      M = imfiltermtx(F, sz, bndry);
      % Suppress rows of M outside of the valid area
      M = M(valid, :);
      
    otherwise
      % Full convolution; return everything
      M = imfiltermtx(F, sz, bndry);
      
  end
