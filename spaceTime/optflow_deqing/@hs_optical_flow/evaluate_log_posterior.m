function L = evaluate_log_posterior(this, uv)
%EVALUATE_LOG_POSTERIOR computes the log-posterior (negative energy) of the
%   flow fields UV 
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

% spatial term

S = {[1 -1], [1; -1]};
p = 0;
for i = 1:length(S)
    u_  = conv2(uv(:,:,1), S{i}, 'valid');
    v_  = conv2(uv(:,:,2), S{i}, 'valid');
    p   = p - sum(u_(:).^2) - sum(v_(:).^2);
end;

% data term
It  = partial_deriv(this.images, uv, this.interpolation_method, this.deriv_filter);
l   = -sum(It(:).^2);

L = this.lambda/this.sigmaS2*p + l/this.sigmaD2;

if this.display
    fprintf('spatial\t%3.2e\tdata\t%3.2e\n', this.lambda*p/this.sigmaS2, l/this.sigmaD2);
end;
