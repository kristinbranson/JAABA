function duv = compute_flow_increment_cgd(this,uv);

%DUV = compute_flow_increment_cgd(THIS, UV) computes the flow increment
%   to the current flow estimate uv by minimizing a quadratic function
%   using the cgd minimizer
%   E(du, dv) = \sum (It +Ix.du + Iy.*dv).^2 + \sum |Dx (u+du)|^2 +     
%               |Dy (u+du)|^2 + |Dx (v+dv)|^2 + |Dy (v+dv)|^2
%               Dx and Dy are horizontal and vertical derivative operator
%
% It's slower and less accurate than backslash 
%
%   This is a member function of the class 'hs_optical_flow'. 
%
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-11-30 $
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
  
  lambda = this.lambda*this.sigmaD2/this.sigmaS2;
  
  duv = zeros(size(uv));  
  duv = duv(:);
  [duv  E]= minimize(duv, @evaluate_flow_increment_energy, 1000, It, Ix, Iy, uv, lambda);  

  
  duv = reshape(duv, size(uv));

  
  function [E dE] = evaluate_flow_increment_energy(duv, It, Ix, Iy, uv, lambda)
  % evaluate the energy function and its derivative to the flow increment
  
  duv = reshape(duv, size(uv));
  
  % Energy function -- data term  
  tmp = It + Ix.*duv(:,:,1) + Iy.*duv(:,:,2);
  
  d   = sum(tmp(:).^2);
  
  % Energy function -- prior (spatial) term  
  S = {[1 -1], [1; -1]};
  p = 0;
  for i = 1:length(S)
      u_{i}  = conv2(uv(:,:,1)+duv(:,:,1), S{i}, 'valid');
      v_{i}  = conv2(uv(:,:,2)+duv(:,:,2), S{i}, 'valid');
      p   = p + sum(u_{i}(:).^2) + sum(v_{i}(:).^2);
  end;

  E = lambda*p + d;
      
  if nargout == 2
      % Energy gradient -- data term
      dE1 = [tmp(:).*Ix(:); tmp(:).*Iy(:)];
      
      % Energy gradient -- prior (spatial) term  
      dEu = zeros(size(uv,1), size(uv,2));
      dEv = dEu;
      for i = 1:length(S)
          du  = conv2(u_{i}, S{i}(end:-1:1), 'full');
          dv  = conv2(v_{i}, S{i}(end:-1:1), 'full');
          dEu = dEu + du;
          dEv = dEv + dv;
      end;
      
      dE = dE1 + lambda*[dEu(:); dEv(:)];
      dE = 2*dE;    % analytical formula
  end;
      
