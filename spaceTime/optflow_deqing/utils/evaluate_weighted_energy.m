function E = evaluate_weighted_energy(uv, im, occ, bfhsz, mfsz, sigma_i)

% edge region: weighted median filtering, the weights are determined by
%  spatial distance, intensity distance, occlusion state
% smooth region: 

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


sigma_x = 7;   %  spatial distance (7)

dilate_sz = 5*[1 1];  % dilation window size for flow edge region [5 5]

sz = size(im);
sz = sz(1:2);

if nargin < 3
    occ = ones(sz);
end;

if nargin < 4
    bfhsz = 10; % half window size
end;

%mfsz = [7 7]; % for test

if nargin < 5
    uvo = uv; 
else
%     uvo(:,:,1) = medfilt2(uv(:,:,1), mfsz, 'symmetric');
%     uvo(:,:,2) = medfilt2(uv(:,:,2), mfsz, 'symmetric');
end;

if nargin < 6
    sigma_i = 5; %3.5; % intensity distance (10) 5 better than 10 and 20
end;

% % % WMF first then MF
% % uvo = uv;

e1 = edge(uv(:,:,1), 'sobel');
e2 = edge(uv(:,:,2), 'sobel');
e  = e1|e2;
mask = imdilate(e, ones(dilate_sz) );

% nearest 4 neighbors smaller than [3 3] 
% tmp = [ 0 1 0; 1 1 1; 0 1 0];
% mask = imdilate(e, tmp );

% mask = e; % no dilation

% Update non-edge regions first
% uv(repmat(mask, [1 1 2])) = uvo(repmat(mask, [1 1 2]));

% Uncomment below to apply WMF to all regions
% mask = ones(size(im));

% Select boundary regions
[indx_row, indx_col] = find(mask ==1); % 
 
% [H W] = size(im);
pad_u  = padarray(uv(:,:,1), bfhsz*[1 1], 'symmetric', 'both');        
pad_v  = padarray(uv(:,:,2), bfhsz*[1 1], 'symmetric', 'both');        
pad_im = padarray(im, bfhsz*[1 1], 'symmetric', 'both');        
pad_occ= padarray(occ, bfhsz*[1 1], 'symmetric', 'both');        

% fprintf('%d\n', length(indx_row));

% Divide into several groups for memory reasons ~70,000 causes out of memory

Indx_Row = indx_row;
Indx_Col = indx_col;
N        = length(Indx_Row); % number of elements to process
n        = 4e4;              % number of elements per batch
nB       = ceil(N/n);

E = 0;

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
        ic = pad_im(rc, cc, :);       % intensity of the center pixel
        
        tmp_u = pad_u(r1:r2, c1:c2);
        tmp_v = pad_v(r1:r2, c1:c2);
        tmp_i = pad_im(r1:r2, c1:c2, :);
        
        % spatial weight
        [C R] = meshgrid(c1:c2, r1:r2);
        w = exp( -((C-cc).^2+(R-rc).^2)/2/sigma_x^2 );
        % Uncomment below: no spatial weight for test
        % w = ones(size(w));
        
        % intensity weight; comment below to disable the term
        w = w.* mean( exp(- (tmp_i - repmat(ic, [size(tmp_i,1) size(tmp_i, 2) 1])).^2/2/sigma_i^2), 3);
        
        % occluded weight;  comment below to disable the term
        w = w.*pad_occ(r1:r2, c1:c2);
        
        % Normalize
        w = w/sum(w(:));
        
%         neighbors_u(:,i) = tmp_u(:);
%         neighbors_v(:,i) = tmp_v(:);
%         weights(:,i)     = w(:);       
        
        E = E + sum( abs(tmp_u(:) - pad_u(rc, cc)).* w(:));
        E = E + sum( abs(tmp_v(:) - pad_v(rc, cc)).* w(:));
        
    end;
    
    
end;

% for 5x5 equal weight

% Select non-boundary regions
[indx_row, indx_col] = find(mask ~=1);
bfhsz = floor(mfsz(1)/2);

pad_u  = padarray(uv(:,:,1), bfhsz*[1 1], 'symmetric', 'both');        
pad_v  = padarray(uv(:,:,2), bfhsz*[1 1], 'symmetric', 'both');        
pad_im = padarray(im, bfhsz*[1 1], 'symmetric', 'both');        
pad_occ= padarray(occ, bfhsz*[1 1], 'symmetric', 'both');        

% fprintf('%d\n', length(indx_row));

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
        ic = pad_im(rc, cc, :);       % intensity of the center pixel
        
        tmp_u = pad_u(r1:r2, c1:c2);
        tmp_v = pad_v(r1:r2, c1:c2);
        tmp_i = pad_im(r1:r2, c1:c2, :);
        
        % spatial weight
        [C R] = meshgrid(c1:c2, r1:r2);
        w = exp( -((C-cc).^2+(R-rc).^2)/2/sigma_x^2 );
        % Uncomment below: no spatial weight for test
        w = ones(size(w));
        
%         % intensity weight; comment below to disable the term
%         w = w.* mean( exp(- (tmp_i - repmat(ic, [size(tmp_i,1) size(tmp_i, 2) 1])).^2/2/sigma_i^2), 3);
%         
%         % occluded weight;  comment below to disable the term
%         w = w.*pad_occ(r1:r2, c1:c2);
        
        % Normalize
        w = w/sum(w(:));
        
%         neighbors_u(:,i) = tmp_u(:);
%         neighbors_v(:,i) = tmp_v(:);
%         weights(:,i)     = w(:);       
        
        E = E + sum( abs(tmp_u(:) - pad_u(rc, cc)).* w(:));
        E = E + sum( abs(tmp_v(:) - pad_v(rc, cc)).* w(:));
        
    end;
    
    
end;
