function I = denoise_robust2(In, lambda, method, num_iter)
% function I = denoise_robust(In, lambda, gt)
%
% Denoise with a robust function by gradient descent
% (I-In)^2 + lambda*rho(f*I), where f is filter
%
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2009$
%   $Revision $
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


% Convert input image to double.
im = double(In);

% PDE (partial differential equation) initial condition.
diff_im = im;

filters = {[1 -1], [1; -1], [1 0; 0 -1], [0 1; -1 0], ...
           [1 0 0; 0 0 0; 0 0 -1], [ 1 0 ; 0 0; 0 -1], [1; 0; -1], [0 1; 0 0; -1 0], ...
           [0 0 1; 0 0 0; -1 0 0], [1 0 0; 0 0 -1], [ 1 0 -1], [0 0 -1; 1 0 0]};


if nargin < 3
%     method = 'charbonnier';
    method = 'lorentzian';
end;
if  nargin <4
    num_iter = 50;
end;


switch method
    case 'geman_mcclure'
        rho = robust_function('geman_mcclure', 0.3);
    case 'lorentzian'
        rho = robust_function('lorentzian', 0.1);
    case 'charbonnier'
        rho = robust_function('charbonnier', 1e-3);
end;

% need minimize from
% http://www.kyb.tuebingen.mpg.de/bs/people/carl/code/minimize/
diff_im = minimize(diff_im(:), @evaluate_denoise_energy, num_iter, im, filters, lambda, rho);

I = reshape(diff_im, size(im));


function [e g] = evaluate_denoise_energy(I, im, filters, lambda, rho)
% evaluate energy and gradient

% energy - data 
I = reshape(I, size(im));
e = sum((I(:)-im(:)).^2);

for iF = 1:length(filters)
    I_  = conv2(I, filters{iF}, 'valid');
    e   = e+ lambda* sum(evaluate(rho, I_(:)));
end;
if nargout == 2
    g = zeros(size(im));
    for iF = 1:length(filters)
        I_  = conv2(I, filters{iF}, 'valid');
        rI_ = reshape(deriv(rho, I_(:)), size(I_));
        
        % To check
        invF = reshape(filters{iF}(end:-1:1), size(filters{iF}) );
        
        g  = g + conv2(rI_, invF, 'full');
    end;
    
    g  = 2*(I -im) + lambda*g;
    g  = g(:);    
end;
