function flow_operator_test(this)
%FLOW_OPERATOR_TEST compare the analytical flow operator with numerical
%   approximated derivatives w.r.t. flow fields 
%
%   This is a member function of the class 'ba_optical_flow'. 
%
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007- $
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


this.display = false;
imsz        = [20 20];
this.images = randn([imsz 2]);
uv          = max(-3, min(3, randn([imsz 2]))); 
this.interpolation_method = 'bi-cubic'; 
this.blend   = 1;


%% -- Black & Anandan
% Analytical results
[It Ix Iy] = partial_deriv(this.images, uv, this.interpolation_method, this.deriv_filter, this.blend);
[A, b] = flow_operator(this, uv, zeros(size(uv)), It, Ix, Iy);
Duv    = b;

% Numerical approximation
delta = 1E-6;

Duv2 = zeros(size(Duv));

for i = 1:prod(imsz)*2
    
    uvp = uv; uvm = uv;
    uvp(i) = uvp(i) + delta;
    uvm(i) = uvm(i) - delta;

    L1 = evaluate_log_posterior(this, uvp);
    L2 = evaluate_log_posterior(this, uvm);

    Duv2(i) = (L1-L2) / (2*delta);
end;

err     = (Duv2-Duv);
nerr    = err/max(max(abs(Duv)), max(abs(Duv2)) );
disp('maximum absolute error and maximum absolute relative error -- Black & Anandan');
[max(abs(err(:)))  max(abs(nerr(:)))]    % should be around or less to 1e-6


if max(abs(err(:))) > 1e-3 | max(abs(nerr(:))) > 1e-3
    figure;             % set a breakpoint here
    tmp1 = reshape(err(1:end/2), imsz);
    tmp2 = reshape(err(1:end/2), imsz);
    imshow([tmp1 tmp2], []);
    pause; close;
end;