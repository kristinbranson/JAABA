function uv = compute_flow(this, init, gt)
%
%COMPUTE_FLOW   Compute flow field
%   UV = COMPUTE_FLOW(THIS[, INIT]) computes the flow field UV with
%   algorithm THIS and the optional initialization INIT.
%  
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


  % Frame size
  sz = [size(this.images, 1), size(this.images, 2)];

  % If we have no initialization argument, initialize with all zeros
  if (nargin < 2)
    uv = zeros([sz, 2]);
  else
    uv = init;
  end
  
  uvhat = uv;  
   
  % Perform structure texture decomposition to get the texture component
  if this.texture == true;
      images  = structure_texture_decomposition_rof( this.images); 
  else
      images  = scale_image(this.images, 0, 255);      
  end;
    
  % Automatic determine pyramid level
  this.pyramid_levels  =  1 + floor( log(min(size(images, 1), size(images,2))/16) / log(this.pyramid_spacing) );  
  
  % Construct image pyramid, using filter setting in Bruhn et al in "Lucas/Kanade.." (IJCV2005') page 218
  
  % For gnc stage 1    
  factor            = sqrt(2);  % sqrt(3)
  smooth_sigma      = sqrt(this.pyramid_spacing)/factor;   % or sqrt(3) recommended by Manuel Werlberger 
  f                 = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);        
  pyramid_images    = compute_image_pyramid(images, f, this.pyramid_levels, 1/this.pyramid_spacing);

  % For segmentation purpose
  org_pyramid_images = compute_image_pyramid(this.images, f, this.pyramid_levels, 1/this.pyramid_spacing);
  
  % For gnc stage 2 to gnc_iters  
  smooth_sigma      = sqrt(this.gnc_pyramid_spacing)/factor;
  f                 = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);        
  gnc_pyramid_images= compute_image_pyramid(images, f, this.gnc_pyramid_levels, 1/this.gnc_pyramid_spacing); 

  % For segmentation purpose
  org_gnc_pyramid_images = compute_image_pyramid(this.images, f, this.gnc_pyramid_levels, 1/this.gnc_pyramid_spacing); 
  
  tic


  for ignc = 1:this.gnc_iters     
      
      
      if ignc == this.gnc_iters
          this.replacment = false;
      else
          this.replacment = true;
      end;

      if this.display
          disp(['GNC stage: ', num2str(ignc)])
      end
          
      if ignc == 1
          pyramid_levels = this.pyramid_levels;
          pyramid_spacing = this.pyramid_spacing;
      else
          pyramid_levels = this.gnc_pyramid_levels;
          pyramid_spacing = this.gnc_pyramid_spacing;
      end;      
     
      % Iterate through all pyramid levels starting at the top
      for l = pyramid_levels:-1:1

          if this.display
              disp(['-Pyramid level: ', num2str(l)])
          end

          % Generate copy of algorithm with single pyramid level and the
          % appropriate subsampling
          small = this;
          if ignc == 1
              nsz               = [size(pyramid_images{l}, 1) size(pyramid_images{l}, 2)];
              small.images      = pyramid_images{l};                            
%               small.max_linear  = 1;             % number of linearization performed per warping, 1 OK for quadratic formulation
              im1   = org_pyramid_images{l}(:,:,1);
          else
              small.images         = gnc_pyramid_images{l};
              nsz   = [size(gnc_pyramid_images{l}, 1) size(gnc_pyramid_images{l}, 2)];
              im1   = org_gnc_pyramid_images{l}(:,:,1);
          end;
          
          % Turn off the quadratic constaint term in the very beginning
          if l == pyramid_levels && ignc == 1
              small.qterm = false;
          end;
          
          % Rescale the flow field
          uv        = resample_flow(uv, nsz);
          uvhat     = resample_flow(uvhat, nsz);
              
%           [tmp small.seg] = msseg(im1,7,6.5, round(50/pyramid_spacing.^(l-1)) );
%           figure; imagesc(small.seg);
%           small.seg = im1;
          small.seg = small.images(:,:,1);
          
          % Run flow method on subsampled images
          [uv uvhat] = compute_flow_base(small, uv, uvhat);
      end
      
      
      % Update GNC parameters (linearly)
      if this.gnc_iters > 1
          new_alpha  = 1 - ignc / (this.gnc_iters-1);
          this.alpha = min(this.alpha, new_alpha);
          this.alpha = max(0, this.alpha);          
      end;

      if true; %this.display

          %fprintf('log-posterior\t%3.3e\n', evaluate_log_posterior(this, uv));
          fprintf('\nlog-posterior uv uvhat \t%3.3e\n', evalaute_new_log_posterior(small, uv, uvhat));
          fprintf('log-posterior uv \t%3.3e\n', evalaute_new_log_posterior(small, uv));
          fprintf('log-posterior uvhat \t%3.3e\n', evalaute_new_log_posterior(small, uvhat));
          
          if nargin == 3
              [aae stdae aepe] = flowAngErr(gt(:,:,1), gt(:,:,2), uv(:,:,1), uv(:,:,2), 0);        % ignore 0 boundary pixels
              fprintf('AAE %3.3f STD %3.3f average end point error %3.3f \n', aae, stdae, aepe);
              
              [aae stdae aepe] = flowAngErr(gt(:,:,1), gt(:,:,2), uvhat(:,:,1), uvhat(:,:,2), 0);        % ignore 0 boundary pixels
              fprintf(' uvhat AAE %3.3f STD %3.3f average end point error %3.3f \n', aae, stdae, aepe);
          end;

          %figure; imshow(flowToColor(uv), []);          
          %figure; imshow(flowToColor(uvhat), []);   title('axu');
      end;

      toc/60
  end;
  

  uv = uvhat;