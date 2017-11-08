function uv = compute_flow(this, init, gt)
%
%COMPUTE_FLOW   Compute flow field
%   UV = COMPUTE_FLOW(THIS[, INIT]) computes the flow field UV with
%   algorithm THIS and the optional initialization INIT.
%  
%   This is a member function of the class 'hs_optical_flow'. 
%
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-10-30 $
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

  % Frame size
  sz = [size(this.images, 1), size(this.images, 2)];

  % If we have no initialization argument, initialize with all zeros
  if (nargin < 2)
    uv = zeros([sz, 2]);
  else
    uv = init;
  end

  % Perform structure texture decomposition to get the texture component
  if this.texture == true;
      % Texture constancy
      
      if size(this.images, 4) ==1
          images  = structure_texture_decomposition_rof( this.images); %, 1/8, 100, this.alp);
      else
          
          for i = 1:size(this.images,4)
              images(:,:,i,:)  = structure_texture_decomposition_rof( squeeze(this.images(:,:,i,:)));
          end;
          
      end;

  else
      images  = scale_image(this.images, 0, 255);      
  end;

  % Automatic determine pyramid level
  this.pyramid_levels  =  1 + floor( log(min(size(images, 1), size(images,2))/16) / log(this.pyramid_spacing) );  
  
  % Construct image pyramid, using setting in Bruhn et al in  "Lucas/Kanade.." (IJCV2005') page 218

  factor            = sqrt(2);  % sqrt(3)
  smooth_sigma      = sqrt(this.pyramid_spacing)/factor;   % or sqrt(3) recommended by Manuel Werlberger   
  f                 = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);    

  pyramid_images    = compute_image_pyramid(images, f, this.pyramid_levels, 1/this.pyramid_spacing);
   
  % Iterate through all pyramid levels starting at the top
  all_flow = {};
  for l = this.pyramid_levels:-1:1
         
    if this.display
        disp(['Pyramid level: ', num2str(l)])    
    end
    
    % Scale flow to the current pyramid level
    uv    =  resample_flow(uv, [size(pyramid_images{l}, 1) size(pyramid_images{l}, 2)]);
   
    % Generate copy of algorithm with single pyramid level and the appropriate subsampling
    small = this;
    small.pyramid_levels = 1;
    small.images         = pyramid_images{l};
    
    % Run flow method on subsampled images
    uv = compute_flow_base(small, uv);    
    all_flow{end+1} = uv;
  end
        
  if this.display
      fprintf('energy of solution \t%3.3e\n', -evaluate_log_posterior(small, uv));
  end;
      
  if ~isempty(this.median_filter_size)
      uv(:,:,1) = medfilt2(uv(:,:,1), this.median_filter_size, 'symmetric');
      uv(:,:,2) = medfilt2(uv(:,:,2), this.median_filter_size, 'symmetric');
  end;