function L = evalaute_new_log_posterior(this, uv, uvhat)

% Test code
% model = load_of_method_ECCV_March('alt-ba');
% model.images = [];
% model.lambda = 5;
% model.lambda2 = 1e1;
% model.weightRatio = 10;
% model.itersLO     =10;

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


% compute auxiliary flow field uv_bar,
% compute joint energy of E(uv) + lambda2 |uv-uvhat|^2 + E_area(uvhat)

if nargin < 3
    uvhat = uv;
end;

L1 = evaluate_log_posterior(this, uv);


E3 = sum((uv(:)-uvhat(:)).^2); 

fsz = this.median_filter_size(1);
hsz = round(fsz/2);
mask = zeros(fsz,fsz);
mask(hsz, hsz) = 1;

filters = {[1 -1], [1; -1], [1 0; 0 -1], [0 1; -1 0], ...
    [1 0 0; 0 0 0; 0 0 -1], [ 1 0 ; 0 0; 0 -1], [1; 0; -1], [0 1; 0 0; -1 0], ...
    [0 0 1; 0 0 0; -1 0 0], [1 0 0; 0 0 -1], [ 1 0 -1], [0 0 -1; 1 0 0]};

E4 = 0;

method = 'charbonnier'; %'geman_mcclure'; %
pen = robust_function(method, 1e-3); % 6.3

for iF = 1:length(filters)
    
    tmp = conv2(uvhat(:,:,1), filters{iF});
    
    %E4  = E4 + sum(evaluate(pen, tmp(:))); % charbonnier
    E4  = E4 + sum(abs(tmp(:)));            % L1
    
    tmp = conv2(uvhat(:,:,2), filters{iF});
    %E4  = E4 + sum(evaluate(pen, tmp(:))); % charbonnier
    E4  = E4 + sum(abs(tmp(:)));            % L1
end;

% [L1 this.lambda2*E3  this.lambda2*E4/this.weightRatio]/1e5
L = L1 - this.lambda2*E3 -this.lambda3* E4;
