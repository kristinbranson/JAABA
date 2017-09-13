function out = resample_flow(uv, sz, method)
% function out = resample_flow(uv, factor, method)
%RESAMPLE_FLOW   Resample flow field
%   OUT = RESAMPLE_FLOW(IN, FACTOR[, METHOD]) resamples (resizes) the flow
%   field IN using a factor of FACTOR.  The optional argument METHOD
%   specifies the interpolation method ('bilinear' (default) or
%   'bicubic'). 
%  
%   This is a private member function of the class 'clg_2d_optical_flow'. 
%
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
  
  % Make bilinear the default method
  if (nargin < 3)
    method = 'bilinear';
%     method = 'bicubic';
  end
%   
%   % Resize u and v 
%   tmp_u = ximresize(uv(:, :, 1), factor, method);
%   tmp_v = ximresize(uv(:, :, 2), factor, method);
%   out = cat(3, tmp_u, tmp_v)*factor;
%   
  ratio = sz(1) / size(uv,1);
  u     = imresize(uv(:,:,1), sz, method)*ratio;
  v     = imresize(uv(:,:,2), sz, method)*ratio;
  out   = cat(3, u, v); 



  
