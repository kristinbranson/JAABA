function L = evalaute_wa_log_posterior(this, uv, uvhat, occ, bfhsz, mfsz)
% compute joint energy of E(uv) + lambda2 |uv-uvhat|^2 +
% E_{weight_area}(uvhat)
%
% This is a member function of the class 'classic_nl_optical_flow'. 
%
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

if nargin < 3
    uvhat = uv;
end;

L1 = evaluate_log_posterior(this, uv);

E3 = sum((uv(:)-uvhat(:)).^2); 

E4 = 0;

sigma_x = 7;   %  spatial distance (7)
sigma_i = 5; %3.5; % intensity distance (10) 5 better than 10 and 20

dilate_sz = 5*[1 1];  % dilation window size for flow edge region [5 5]

if nargin < 3
    occ = ones(size(im));
end;

if nargin < 4
    bfhsz = 10; % half window size
end;

uvo = uv;

e1 = edge(uv(:,:,1), 'sobel');
e2 = edge(uv(:,:,2), 'sobel');
e  = e1|e2;
mask = imdilate(e, ones(dilate_sz) );


[indx_row, indx_col] = find(mask ==1);
 
% [H W] = size(im);

pad_u  = padarray(uv(:,:,1), bfhsz*[1 1], 'symmetric', 'both');        
pad_v  = padarray(uv(:,:,2), bfhsz*[1 1], 'symmetric', 'both');        
pad_im = padarray(im, bfhsz*[1 1], 'symmetric', 'both');        
pad_occ= padarray(occ, bfhsz*[1 1], 'symmetric', 'both');        

% Divide into several groups for memory reasons ~70,000 causes out of memory

Indx_Row = indx_row;
Indx_Col = indx_col;
N        = length(Indx_Row); % number of elements to process
n        = 4e4;              % number of elements per batch
nB       = ceil(N/n);

for ib = 1:nB;
    istart = (ib-1)*n + 1;
    iend   = min(ib*n, N);
    indx_row = Indx_Row(istart:iend);
    indx_col = Indx_Col(istart:iend);
    
    neighbors_u = zeros((bfhsz*2+1)^2, length(indx_row));
    neighbors_v = neighbors_u;
    weights     = neighbors_u;
    
    for i = 1:length(indx_row)
        
        % crop window
        r1 = indx_row(i);
        r2 = indx_row(i) + 2*bfhsz;
        c1 = indx_col(i);
        c2 = indx_col(i) + 2*bfhsz;
        
        rc = indx_row(i) + bfhsz; % row  center
        cc = indx_col(i) + bfhsz; % column center
        ic = pad_im(rc, cc);       % intensity of the center pixel
        
        tmp_u = pad_u(r1:r2, c1:c2);
        tmp_v = pad_v(r1:r2, c1:c2);
        tmp_i = pad_im(r1:r2, c1:c2);
        
        % spatial weight
        [C R] = meshgrid(c1:c2, r1:r2);
        w = exp( -((C-cc).^2+(R-rc).^2)/2/sigma_x^2 );
        % Uncomment below: no spatial weight for test
        % w = ones(size(w));
        
        % intensity weight; comment below to disable the term
        w = w.* exp(- (tmp_i - ic).^2/2/sigma_i^2);
        
        % occluded weight;  comment below to disable the term
        w = w.*pad_occ(r1:r2, c1:c2);
        
        % Normalize
        
        w = w/sum(w(:));
        
        neighbors_u(:,i) = tmp_u(:);
        neighbors_v(:,i) = tmp_v(:);
        weights(:,i)     = w(:);
                
    end;
    
    % solve weighted median filtering
    indx = sub2ind(size(im), indx_row, indx_col);
    uo   = uv(:,:,1);
    tmp  = uo(indx);
    tmp  = tmp(:)'; 
    
    [H W] = size(neighbors_u);
    
    E4 = E4 + sum(sum(abs(neighbors_u - repmat(tmp, [H, 1]) ) ) ); 
    
    vo   = uv(:,:,2);
    tmp  = vo(indx);
    tmp  = tmp(:)'; 
    E4 = E4 + sum(sum(abs(neighbors_v - repmat(tmp, [H, 1]) ) ) ); 
end;

% energy of area term in nonboundary regions

% replicate

mfsize = median_filter_size(1);
hfsize = floor(mfsize/2);

u = padarray(uv(:,:,1), hfsize*[1 1], 'symmetric', 'both');
v = padarray(uv(:,:,2), hfsize*[1 1], 'symmetric', 'both');

tmpu = sum( abs(u - repmat(u(round(mfsize*mfsize/2), :), [size(u,1) 1])) );
tmpu(indx) = 0;
E4 = E4 + sum(tmpu);

tmpv = sum( abs(v - repmat(v(round(mfsize*mfsize/2), :), [size(v,1) 1])) );
tmpv(indx) = 0;
E4 = E4 + sum(tmpv);

% [L1 this.lambda2*E3  this.lambda2*E4/this.weightRatio]/1e5
L = L1 - this.lambda2*E3 -this.lambda3* E4;
