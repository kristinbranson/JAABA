function [It, Ix, Iy] = partial_deriv(images, uv_prev, interpolation_method, deriv_filter, b)

%PARTIAL_DERIV   Spatio-temporal derivatives
%   P = PARTIAL_DERIV(IMAGES, INIT) computes the spatio-temporal derivatives
%
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-11-30 $
%   $Revision $
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

if nargin == 2
    interpolation_method = 'cubic';
end;

if nargin == 4
    h = deriv_filter;
else
    h = [1 -8 0 8 -1]/12; % used in Wedel etal "improved TV L1"
end;

if nargin < 5
    b = 0.5;    %blending ratio
end;

H   = size(images, 1);
W   = size(images, 2);

[x,y]   = meshgrid(1:W,1:H);
x2      = x + uv_prev(:,:,1);        
y2      = y + uv_prev(:,:,2);  

% Record out of boundary pixels
B = (x2>W) | (x2<1) | (y2>H) | (y2<1);

if size(images, 4) ~= 1
    B = repmat(B, [1 1 size(images, 3)]);
end;

if size(images, 4) == 1
    img1 = images(:,:,1);
    img2 = images(:,:,2);
else
    img1 = images(:,:,:,1);
    img2 = images(:,:,:,2);
end;

if strcmp(interpolation_method, 'bi-cubic')

    %h = [-0.5, 0, 0.5]; % consistent with bi-cubic interpolation        
    
    % Bicubic interpolation
    if nargout == 1
        
        if size(images, 4) == 1
            % gray-level
            warpIm = interp2_bicubic(images(:,:,2),x2,y2, h);
        else
            % color
            warpIm  = zeros(size(images(:,:,:,1)));
            for j = 1:size(images,3)
                warpIm(:,:,j) = interp2_bicubic(images(:,:,j, 2),x2,y2, h);
            end;
        end;
        
    elseif nargout == 3
        
        if size(images, 4) == 1
            % gray-level
            [warpIm Ix Iy] = interp2_bicubic(images(:,:,2),x2,y2, h);
        else
            % color
            warpIm  = zeros(size(images(:,:,:,1)));
            Ix      = warpIm;
            Iy      = warpIm;
            for j = 1:size(images,3)
                [warpIm(:,:,j) Ix(:,:,j) Iy(:,:,j)] = interp2_bicubic(images(:,:,j,2),x2,y2, h);
            end;
        end;
    
    else
        error('partial_deriv: number of output wrong!');
    end;

    indx        = isnan(warpIm);
    if size(images, 4) == 1
        It          = warpIm - images(:,:,1);
    else
        It          = warpIm - images(:,:,:,1);
    end;
    
    % Disable those out-of-boundary pixels in warping
    It(indx)    = 0;
    if nargout == 3        
        
        % Temporal average
        I1x = imfilter(img1, h,  'corr', 'symmetric', 'same');  %
        I1y = imfilter(img1, h', 'corr', 'symmetric', 'same');
        
        % b = 0.6;    % recommended in Wedal etal "improved TV-L1" 2008
        % b = 0.5;
%         b = 1;
        
        
        Ix  = b*Ix+(1-b)*I1x;
        Iy  = b*Iy+(1-b)*I1y;

        Ix(indx) = 0;
        Iy(indx) = 0;       
    end;
    
elseif strcmp(interpolation_method, 'bi-linear') || strcmp(interpolation_method, 'cubic') 

    if strcmp(interpolation_method, 'bi-linear')
        method = 'linear';
    else
        method = 'cubic';
    end;
    
    % Matlab built-in 
    if size(images, 4) == 1
        % Gray-level
        warpIm  = interp2(x,y,images(:,:,2),x2,y2,method);
        It      = warpIm - images(:,:,1);      
        
    else
        % Color
        warpIm  = zeros(size(images(:,:,:,1)));
        for j = 1:size(images,3)
            warpIm(:,:,j) = interp2(x,y,images(:,:,j,2),x2,y2,method, NaN);
        end;
        It      = warpIm - images(:,:,:,1);
        
    end;

        % spline-based bicubic interpolation code with incorrect derivatives (below)
%         if size(images, 4) == 1
%             % gray-level
%             warpIm = interp2_bicubic(images(:,:,2),x2,y2, h);
%             It      = warpIm - images(:,:,1);      
%         else
%             % color
%             warpIm  = zeros(size(images(:,:,:,1)));
%             for j = 1:size(images,3)
%                 warpIm(:,:,j) = interp2_bicubic(images(:,:,j, 2),x2,y2, h);
%             end;
%             It      = warpIm - images(:,:,:,1);
%         end;
        
    if nargout == 3
         
        % First compute derivative then warp
        I2x = imfilter(img2, h,  'corr', 'symmetric', 'same');
        I2y = imfilter(img2, h', 'corr', 'symmetric', 'same');
        
        if size(images, 4) == 1
            % Gray-level
            Ix  = interp2(x,y,I2x,x2,y2,method);
            Iy  = interp2(x,y,I2y,x2,y2,method);

%             % spline-based bicubic interpolation code with incorrect derivatives            
%             Ix  = interp2_bicubic(I2x,x2,y2, h);
%             Iy  = interp2_bicubic(I2y,x2,y2, h);
            
        else
            % Color
            Ix  = zeros(size(images(:,:,:,1)));
            Iy  = Ix;
            
            for j = 1:size(images,3)
                Ix(:,:,j)  = interp2(x,y,I2x(:,:,j),x2,y2,method);
                Iy(:,:,j)  = interp2(x,y,I2y(:,:,j),x2,y2,method);
            end;
            
        end;
        
        % warp then derivative        
%         if size(images, 4) == 1
%             tmp           = images(:,:,1);
%         else
%             tmp           = images(:,:,:,1);
%         end;
%         warpIm(B) = tmp(B);                
%         Ix        = imfilter(warpIm, h,  'corr', 'symmetric', 'same');  %
%         Iy        = imfilter(warpIm, h', 'corr', 'symmetric', 'same');
        % end of warp then derivative
        

        % Temporal average
        I1x = imfilter(img1, h,  'corr', 'symmetric', 'same');  %
        I1y = imfilter(img1, h', 'corr', 'symmetric', 'same');
        
        % b = 0.6;    % recommended in Wedal etal "improved TV-L1" 2008
        % b = 0.5;
        
        Ix  = b*Ix+(1-b)*I1x;
        Iy  = b*Iy+(1-b)*I1y;        
       
    end;
    
    % Disable those out-of-boundary pixels in warping
    It(B)   = 0;
    if nargout == 3
        Ix(B) = 0;
        Iy(B) = 0;
    end;
    
else
    error('partial_deriv: unknown interpolation method!');
end;