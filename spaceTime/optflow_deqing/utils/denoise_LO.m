function uo = denoise_LO(un, median_filter_size, lambda, niters)
%% Denoising using the Li & Osher median formula

% Y. Li and Osher "A New Median Formula with Applications to PDE Based
% Denoising"
% lambda is on the data term
% lambda |u-un|^2 + \sum_{j\in N_i} |u_i-u_j| 

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



if nargin < 4    
    niters = 1;
end;

mfsize = median_filter_size(1);
hfsize = floor(mfsize/2);

% padding
n   = (mfsize*mfsize-1)/2;
tmp = -n:n;
tmp = repmat(tmp', [1 numel(un(:))])/lambda;
tmp = repmat(un(:)', [2*n+1, 1]) + tmp;

uo  = un;
for i = 1:niters
    
    % replicate
    u  = padarray(uo, hfsize*[1 1], 'symmetric', 'both');        
    u2 = im2col(u, mfsize*[1 1]);   
    u2 = [u2(1:floor(mfsize*mfsize/2),:); u2(ceil(mfsize*mfsize/2) + 1: end, :)];     
    
    uo = reshape(median([u2; tmp]), [size(un,1) size(un,2)]);
end;

