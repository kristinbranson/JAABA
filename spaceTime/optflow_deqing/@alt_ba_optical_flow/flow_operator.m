function [A, b, params, iterative] = flow_operator(this, uv, duv, It, Ix, Iy)

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
%   This is a member function of the class 'alt_ba_optical_flow'. 
%
% Authors: Deqing Sun, Department of Computer Science, Brown University
% Contact: dqsun@cs.brown.edu
% $Date: 2009 $
% $Revision: $
%
% Copyright 2009-2010, Brown University, Providence, RI. USA
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

  sz        = [size(Ix,1) size(Ix,2)];
  npixels   = prod(sz);

  % spatial term
  S = this.spatial_filters;

  FU = sparse(npixels, npixels);
  FV = sparse(npixels, npixels);
  for i = 1:length(S)

      FMi = make_convn_mat(S{i}, sz, 'valid', 'sameswap');
      Fi  = FMi';

      % Use flow increment to update the nonlinearity
      u_        = FMi*reshape(uv(:, :, 1)+duv(:, :, 1), [npixels 1]);      
      v_        = FMi*reshape(uv(:, :, 2)+duv(:, :, 2), [npixels 1]);      

      if isa(this.rho_spatial_u{i}, 'robust_function')          
          pp_su     = deriv_over_x(this.rho_spatial_u{i}, u_);
          pp_sv     = deriv_over_x(this.rho_spatial_v{i}, v_);          
      elseif isa(this.rho_spatial_u{i}, 'gsm_density')
          pp_su     = -evaluate_log_grad_over_x(this.rho_spatial_u{i}, u_')';
          pp_sv     = -evaluate_log_grad_over_x(this.rho_spatial_v{i}, v_')';          
      else
          error('evaluate_log_posterior: unknown rho function!');
      end;
      
      FU        = FU+ Fi*spdiags(pp_su, 0, npixels, npixels)*FMi;
      FV        = FV+ Fi*spdiags(pp_sv, 0, npixels, npixels)*FMi;
      
  end;

  M = [-FU, sparse(npixels, npixels);
      sparse(npixels, npixels), -FV];

  
  % data term   
  Ix2 = Ix.^2;
  Iy2 = Iy.^2;
  Ixy = Ix.*Iy;
  Itx = It.*Ix;
  Ity = It.*Iy;
  
  % Perform linearization - note the change in It
  It = It + Ix.*repmat(duv(:,:,1), [1 1 size(It,3)]) ...
          + Iy.*repmat(duv(:,:,2), [1 1 size(It,3)]);
  
      
%   figure; imshow(It, []); %temporarily added to observe BCE
  
  if isa(this.rho_data, 'robust_function')      
      pp_d  = deriv_over_x(this.rho_data, It(:));      
  elseif isa(this.rho_data, 'gsm_density')
      pp_d = -evaluate_log_grad_over_x(this.rho_data, It(:)')';
  else
      error('flow_operator: unknown rho function!');
  end;  
  
  % Display estimated data outlier
  % figure;  imshow(reshape(pp_d, sz), []);
  
  
  tmp = pp_d.*Ix2(:);
  duu = spdiags(tmp, 0, npixels, npixels);
  tmp = pp_d.*Iy2(:);
  dvv = spdiags(tmp, 0, npixels, npixels);
  tmp = pp_d.*Ixy(:);
  dduv = spdiags(tmp, 0, npixels, npixels);

  A = [duu dduv; dduv dvv] - this.lambda*M;

  % right hand side
  b =  this.lambda * M * uv(:) - [pp_d.*Itx(:); pp_d.*Ity(:)];

  % No auxiliary parameters
  params    = [];
  
  % If the non-linear weights are non-uniform, do more linearization
  if (max(pp_su(:)) - min(pp_su(:)) < 1E-6 && ...
     max(pp_sv(:)) - min(pp_sv(:)) < 1E-6 && ...
      max(pp_d(:)) - min(pp_d(:)) < 1E-6)
    iterative = false;
  else
    iterative = true;
  end