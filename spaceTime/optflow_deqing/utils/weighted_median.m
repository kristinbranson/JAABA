function uo = weighted_median(w, u)
%%
% compute the weighted median 
%   uo = \min_u \sum w(i)|uo - u(i)| 
% using the formula (3.13) in
%   Y. Li and Osher "A New Median Formula with Applications to PDE Based
%   Denoising"
% applied to every corresponding columns of w and u

% Authors: Deqing Sun, Department of Computer Science, Brown University 
% Contact: dqsun@cs.brown.edu
% $Date: $
% $Revision: $
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

%

[H W] = size(u);
[sort_u ir] = sort(u);
ic       = repmat(1:W, [H, 1]);
ind      = sub2ind([H W], ir, ic);

w        = w(ind);

% number of negative infinites for each column
k = ones(1, W);
pp = -sum(w);
for i = H:-1:1;
    pc = pp + 2*w(i, :);
    indx = (pc > 0) & (pp < 0);
    k(indx) = H-i+1;
    pp = pc;
end;

k   = H-k+1;
ind = sub2ind([H W], k, 1:W);

uo  = sort_u(ind);
