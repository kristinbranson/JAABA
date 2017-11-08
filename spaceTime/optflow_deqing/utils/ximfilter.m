function B = ximfilter(A, H, varargin)

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

  bndry  = 0;
  output = 'same';
  type   = 'corr';
  
  % Parse options
  for i = 1:length(varargin)
    if (isnumeric(varargin{i}))
      bndry = varargin{i};
    else
      switch(varargin{i})
        case {'symmetric', 'replicate', 'circular'}
          bndry = varargin{i};
        case {'same', 'full'}
          output = varargin{i};
        case {'corr', 'conv'}
          type = varargin{i};
      end
    end
  end
  
  
  % Pad input A depending on boundary handling
  if (strcmp(output, 'same'))
    padl = floor((size(H) - 1) / 2);
    padh = size(H) - 1 - padl;
  else
    padl = size(H) - 1;
    padh = size(H) - 1;
  end 
  
  if (isnumeric(bndry))
    tmp = xrepmat(bndry(1), size(A) + padl + padh);
    tmp(padl(1)+1:end-padh(1), padl(2)+1:end-padh(2)) = A;
    A = tmp;
    clear tmp;
  else
    switch (bndry)
      case 'symmetric'
        sz = size(A);
        idx1 = [padl(1):-1:1, 1:sz(1), sz(1):-1:sz(1)-padh(1)+1];
        idx2 = [padl(2):-1:1, 1:sz(2), sz(2):-1:sz(2)-padh(2)+1];
        
      case 'replicate'
        sz = size(A);
        idx1 = [xrepmat(1, 1, padl(1)), 1:sz(1), xrepmat(sz(1), 1, padh(1))];
        idx2 = [xrepmat(1, 1, padl(2)), 1:sz(2), xrepmat(sz(2), 1, padh(2))];
      
      case 'circular'
        sz = size(A);
        idx1 = [sz(1)-padl(1)+1:sz(1), 1:sz(1), 1:padh(1)];
        idx2 = [sz(2)-padl(2)+1:sz(2), 1:sz(2), 1:padh(2)];
        
      otherwise
        error('Invalid boundary handling option!');
    end
   
    A = A(idx1, idx2);
  end
  
  % Mirror filter, if we use correlation
  if (strcmp(type, 'corr'))
    H = reshape(H(end:-1:1), size(H));
  end
  
  B = conv2(A, H, 'same');
  switch (output)
    case 'same'
      B = B(padl(1)+1:end-padh(1), padl(2)+1:end-padh(2));
    case 'full'
      padl = floor((size(H) - 1) / 2);
      padh = size(H) - 1 - padl;
      B = B(padl(1)+1:end-padh(1), padl(2)+1:end-padh(2));
    otherwise
      error('Invalid return shape!');
  end