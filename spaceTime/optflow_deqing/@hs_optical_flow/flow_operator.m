function [A, b, params, iterative] = flow_operator(this, uv)
%FLOW_OPERATOR   Linear flow operator (equation) for flow estimation
%   [A, b] = FLOW_OPERATOR(THIS, UV, INIT)  
%   returns a linear flow operator (equation) of the form A * x = b.  The
%   flow equation is linearized around UV with the initialization INIT
%   (e.g. from a previous pyramid level).  
%
%   [A, b, PARAMS, ITER] = FLOW_OPERATOR(...) returns optional parameters
%   PARAMS that are to be passed into a linear equation solver and a flag
%   ITER that indicates whether solving for the flow requires multiple
%   iterations of linearizing.
%  
%   This is a member function of the class 'hs_optical_flow'. 
%   $Revision: $
%
% Copyright 2007-2010, Brown University, Providence, RI. USA
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

  % Compute spatial and temporal partial derivatives
  [It Ix Iy] = partial_deriv(this.images, uv, this.interpolation_method, this.deriv_filter);
  
  sz        = [size(Ix, 1) size(Ix, 2)];
  npixels   = prod(sz);

  % Spatial/prior term
  L     = [0 1 0; 1 -4 1; 0 1 0];     % Laplacian operator
  F     = make_imfilter_mat(L, sz, 'replicate', 'same');

  % Replicate operator for u and v
  M     = [F, sparse(npixels, npixels);
           sparse(npixels, npixels), F];

       
  % For color processing
  Ix2 = mean(Ix.^2, 3);
  Iy2 = mean(Iy.^2, 3);
  Ixy = mean(Ix.*Iy, 3);
  Itx = mean(It.*Ix, 3);
  Ity = mean(It.*Iy, 3);

  duu   = spdiags(Ix2(:), 0, npixels, npixels);
  dvv   = spdiags(Iy2(:), 0, npixels, npixels);
  duv   = spdiags(Ixy(:), 0, npixels, npixels);

  % Compute the operator
  A     = [duu duv; duv dvv]/this.sigmaD2 - this.lambda*M/this.sigmaS2;  
  b     = this.lambda *M*uv(:)/this.sigmaS2 - [Itx(:); Ity(:)]/this.sigmaD2;
  
  % No auxiliary parameters
  params    = [];
  iterative = true;
 
  % If the non-linear weights are non-uniform, we have to iterate
%   if (max(pp_s(:)) - min(pp_s(:)) < 1E-6 && ...
%       max(pp_d(:)) - min(pp_d(:)) < 1E-6)
%     iterative = false;
%   else
%     iterative = true;
%   end