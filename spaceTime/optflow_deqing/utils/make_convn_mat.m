function M = make_convn_mat(F, sz, shape, pad)
%MAKE_CONVN_MAT   Flexible N-D convolution matrix
%   M = MAKE_CONVN_MAT(F, SZ[, SHAPE[, PAD]]) returns the convolution
%   matrix for the matrix F.  SZ gives the size of the array that the
%   convolution should be applied to.  The returned matrix M is sparse. 
%   If X is of size SZ and SHAPE is 'full', then reshape(M * X(:), SZ +
%   size(F) - 1) is the same as convn(X, F).
%
%   The optional parameter SHAPE controls which parts of the convolution
%   result to return in M:
%     * 'full':      (default) returns the full N-D convolution, i.e. the 
%                    behavior is identical to convn(X, F, 'full').
%     * 'same':      returns the central part of the convolution that 
%                    is the same size as X, i.e. the behavior is
%                    identical to convn(X, F, 'same').
%     * 'sameswap':  identical to 'same' except that rounding needed for
%                    even-sized filters is performed the opposite way.
%     * 'valid':     returns only the part of the result that can be
%                    computed without assuming zero-padded arrays,
%                    i.e. the behavior is identical to convn(X, F,
%                    'valid').
%
%   The optional parameter PAD controls whether to return an M that pads
%   the result of the convolution with zeros (default no padding):
%     * 'full':      returns the result padded to the output size of a
%                    full N-D convolution (works with SHAPE set to 'same'
%                    and 'valid'). 
%     * 'same':      returns the result padded to the central part of the
%                    convolution that is the same size as X (works with
%                    SHAPE set to 'valid').
%     * 'sameswap':  identical to 'same' except that rounding needed for
%                    even-sized filters is performed the opposite way.
%  
%   See also CONVMTXN.
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


  
  % Default to full matrix
  if (nargin < 3)
    shape = 'full';
  end
  
  ndims = length(sz);

  % Border sizes for 'same' and 'sameswap'
  Fsize_lo_2 = ceil((size(F) - 1) / 2);
  Fsize_hi_2 = floor((size(F) - 1) / 2);
  
  % Border sizes for 'valid'
  Fsize = size(F) - 1;

  switch(shape)
    case 'same'
      % Mark valid and invalid pixels (i.e. the ones within and outside
      % of the part to be returned)
      valid = true(sz+size(F)-1);
      
      for d = 1:ndims
        for e = 1:ndims
          sub{e} = ':';
        end
        
        sub{d} = 1:Fsize_lo_2(d);      
        valid(sub{:}) = false;
        sub{d} = size(valid, d)-Fsize_hi_2(d)+1:size(valid, d);      
        valid(sub{:}) = false;
      end

      if (nargin > 3 && strcmp(pad, 'full'))
        % If we're padding to 'full' set the coefficients on the border
        % to zero
        M = convmtxn(F, sz, valid);
      else
        % If we're *not* padding, then suppress the rows of M that
        % correspond to the border
        M = convmtxn(F, sz);
        M = M(valid, :);
      end
    
   case 'sameswap'
     % Mark valid and invalid pixels (i.e. the ones within and outside
     % of the part to be returned), but round the other way
     valid = true(sz+size(F)-1);
     
     for d = 1:ndims
       for e = 1:ndims
         sub{e} = ':';
       end
       
       sub{d} = 1:Fsize_hi_2(d);      
       valid(sub{:}) = false;
       sub{d} = size(valid, d)-Fsize_lo_2(d)+1:size(valid, d);      
       valid(sub{:}) = false;
     end
     
     if (nargin > 3 && strcmp(pad, 'full'))
       % If we're padding to 'full' set the coefficients on the border
       % to zero
       M = convmtxn(F, sz, valid);
     else
       % If we're *not* padding, then suppress the rows of M that
       % correspond to the border
       M = convmtxn(F, sz);
       M = M(valid, :);
     end
    
   case 'valid'
     % Mark valid and invalid pixels (i.e. the ones within and outside
     % of the part to be returned)
     valid = true(sz+size(F)-1);
     
     for d = 1:ndims
       for e = 1:ndims
         sub{e} = ':';
       end
       
       sub{d} = 1:Fsize(d);      
       valid(sub{:}) = false;
       sub{d} = size(valid, d)-Fsize(d)+1:size(valid, d);      
       valid(sub{:}) = false;
     end
     
     if (nargin > 3)
       % If we're padding, then figure out the area to be padded       

       switch (pad)
         case 'same'
           % Mark valid and invalid pixels (i.e. the ones within and outside
           % of the part to be padded)
           pad_valid = true(sz+size(F)-1);
           
           for d = 1:ndims
             for e = 1:ndims
               sub{e} = ':';
             end
             
             sub{d} = 1:Fsize_lo_2(d);      
             pad_valid(sub{:}) = false;
             sub{d} = size(valid, d)-Fsize_hi_2(d)+1:size(valid, d);      
             pad_valid(sub{:}) = false;
           end
           
           % Set coefficients on the border to zero
           M = convmtxn(F, sz, valid);
           
           % Suppress rows of M outside of the padded area
           M = M(pad_valid, :);
           
         case 'sameswap'
           % Mark valid and invalid pixels (i.e. the ones within and outside
           % of the part to be padded), but round the other way
           pad_valid = true(sz+size(F)-1);
           
           for d = 1:ndims
             for e = 1:ndims
               sub{e} = ':';
             end
             
             sub{d} = 1:Fsize_hi_2(d);      
             pad_valid(sub{:}) = false;
             sub{d} = size(valid, d)-Fsize_lo_2(d)+1:size(valid, d);      
             pad_valid(sub{:}) = false;
           end
           
           % Set coefficients on the border to zero
           M = convmtxn(F, sz, valid);
           
           % Suppress rows of M outside of the padded area
           M = M(pad_valid, :);
           
         otherwise
           % Padding to 'full'; only set coefficients on the border to zero
           M = convmtxn(F, sz, valid);           
       end
     else
       % No padding; suppress all rows on the border
       M = convmtxn(F, sz);
       M = M(valid, :);
     end
     
    otherwise
      % Full convolution; return everything
      M = convmtxn(F, sz);
      
  end
