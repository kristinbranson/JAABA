function guv = evaluate_log_posterior_grad(this, uv)
%EVALUATE_LOG_POSTERIOR computes the gradient of the log-posterior
%   (negative energy) wrt the flow fields UV 
%   Actually only proportional to the log posterior since the variance of neither the
%   spatial nor the data terms is considered
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

% Spatial term
S       = this.spatial_filters;
gu1     = zeros(size(uv,1), size(uv,2));
gv1     = gu1;

for i = 1:length(S)

    u_ = conv2(uv(:,:,1), S{i}, 'valid');
    v_ = conv2(uv(:,:,2), S{i}, 'valid');

    Si = reshape(S{i}(end:-1:1), size(S{i}));
    
    if isa(this.rho_spatial_u{i}, 'robust_function')        
        u_ = -reshape(deriv(this.rho_spatial_u{i}, u_(:)), size(u_));                
        v_ = -reshape(deriv(this.rho_spatial_v{i}, v_(:)), size(v_));                        
    elseif isa(this.rho_spatial_u{i}, 'gsm_density')
        u_ = reshape(evaluate_log_grad(this.rho_spatial_u{i}, u_(:)'), size(u_));                
        v_ = reshape(evaluate_log_grad(this.rho_spatial_v{i}, v_(:)'), size(v_));                        
    else
        error('evaluate_log_posterior: unknown rho function!');
    end;
    
    gu1 = gu1+conv2(u_, Si, 'full');
    gv1 = gv1+conv2(v_, Si, 'full');    
end;

gu2     = zeros(size(uv,1), size(uv,2));
gv2     = gu2;

% Data term
[It Ix Iy] = partial_deriv(this.images, uv, this.interpolation_method);    
    
if isa(this.rho_data, 'robust_function')
    temp   = -reshape(deriv(this.rho_data, It(:)), size(It));
    
elseif isa(this.rho_data, 'gsm_density')    
    temp   = reshape(evaluate_log_grad(this.rho_data, It(:)'), size(It));
else
    error('evaluate_log_posterior: unknown rho function!');
end;

gu2     = temp.*Ix;
gv2     = temp.*Iy;

guv = cat(3, gu2+this.lambda*gu1, gv2+this.lambda*gv1);
