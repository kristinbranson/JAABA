function uvo = denoise_color_weighted_medfilt2(uv, im, occ, bfhsz, mfsz, sigma_i, fullVersion)

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

if nargin < 5
    uvo = uv; 
else
    uvo(:,:,1) = medfilt2(uv(:,:,1), mfsz, 'symmetric');
    uvo(:,:,2) = medfilt2(uv(:,:,2), mfsz, 'symmetric');
end;

if nargin < 6
    sigma_i = 5; 
end;

if nargin < 7
    fullVersion = false;
end;

% % % WMF first then MF
% % uvo = uv;

e1 = edge(uv(:,:,1), 'sobel');
e2 = edge(uv(:,:,2), 'sobel');
e  = e1|e2;
mask = imdilate(e, ones(dilate_sz) );

% below to apply WMF to all regions
if fullVersion
    mask = ones(size(mask));
end;

[indx_row, indx_col] = find(mask ==1);
 
pad_u  = padarray(uv(:,:,1), bfhsz*[1 1], 'symmetric', 'both');        
pad_v  = padarray(uv(:,:,2), bfhsz*[1 1], 'symmetric', 'both');        
pad_im = padarray(im, bfhsz*[1 1], 'symmetric', 'both');        
pad_occ= padarray(occ, bfhsz*[1 1], 'symmetric', 'both');        


[H W] = size(pad_u);

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

    [C R] = meshgrid(-bfhsz:bfhsz, -bfhsz:bfhsz);
    nindx = R + C*H;    
    cindx = indx_row +bfhsz  + (indx_col+bfhsz-1)*H;
    
    pad_indx = repmat(nindx(:), [1 length(indx_row)]) + ...
               repmat(cindx(:)', [(bfhsz*2+1)^2, 1] );
           
    % spatial weight
    tmp = exp(- (C.^2 + R.^2) /2/sigma_x^2 );
    weights = repmat(tmp(:), [1 length(indx_row)]);    
    
    % %Uncomment below: no spatial weight for test
    % weights = ones(size(weights));    
    
    % Intensity weight
    tmp_w = zeros(size(weights));
    
    for i = 1:size(pad_im,3)
        tmp = pad_im(:,:,i);
        tmp_w = tmp_w + (tmp(pad_indx) - repmat(tmp(cindx(:))', [(bfhsz*2+1)^2, 1])).^2;
    end;    
    tmp_w = tmp_w/size(pad_im,3);
    
    weights = weights.* exp(-tmp_w/2/sigma_i^2);

    % Occlusion weight    
    weights = weights.*pad_occ(pad_indx);
    
    % Normalize
    weights = weights./repmat(sum(weights, 1), [(bfhsz*2+1)^2, 1]);
    
    neighbors_u = pad_u(pad_indx);
    neighbors_v = pad_v(pad_indx);

    % solve weighted median filtering
    indx = sub2ind(sz, indx_row, indx_col);
    uo   = uvo(:,:,1);
    u    = weighted_median(weights, neighbors_u);
%     u    = weighted_median_iter(weights, neighbors_u);
    uo(indx) = u;
    vo   = uvo(:,:,2);
    v    = weighted_median(weights, neighbors_v);
%     v    = weighted_median_iter(weights, neighbors_v);
    vo(indx) = v;
    uvo = cat(3, uo, vo);
    
end;
