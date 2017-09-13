function this = ba_optical_flow(varargin)
%
%BA_OPTICAL_FLOW   Optical flow computation with Black & Anandan method
%       Black, M. J. and Anandan, P. The robust estimation of multiple
%       motions: Parametric and piecewise-smooth flow fields,
%       CVIU, 63(1), pp. 75-104, Jan. 1996.
%       http://www.cs.brown.edu/people/black/Papers/cviu.63.1.1996.pdf
%       
%   BA_OPTICAL_FLOW([IMGS]) constructs a BA optical flow object
%   with the optional image sequence IMGS ([n x m x 2] array). 
%   BA_OPTICAL_FLOW(O) constructs BA optical flow object by copying O.
%  
%   This is a member function of the class 'ba_optical_flow'. 
%
% Authors: Deqing Sun, Department of Computer Science, Brown University
%          Stefan Roth, Department of Computer Science, TU Darmstadt
% Contact: dqsun@cs.brown.edu, sroth@cs.tu-darmstadt.de
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

error(nargchk(0, 1, length(varargin)));
  
  switch (length(varargin))
    case 0
        
      this.images          = [];              
      this.lambda          = 1;
      this.lambda_q        = 1;    % Quadratic formulation of the objective function
      
      this.sor_max_iters   = 1e4;       % 100 seems sufficient

      this.limit_update    = true;      % limit the flow incrment to be less than 1 per linearization step
      this.display         = false;            
      
      this.solver          = 'backslash';   % 'sor' 'pcg' for machines with limited moemory       
      this.deriv_filter    = [1 -8 0 8 -1]/12; % 5-point 7 point [-1 9 -45 0 45 -9 1]/60; 
      this.blend           = 0.5;           % temporal blending ratio
      
      this.texture         = false;     % use texture component as input
      this.fc              = false;     % use filter constancy

      this.median_filter_size   = []; %[5 5];      
     
      this.interpolation_method = 'cubic';  % 'bi-cubic', 'cubic', 'bi-linear'
            
      % For Graduated Non-Convexity (GNC) optimization
      this.gnc_iters       = 3;
      this.alpha           = 1;             % change linearly from 1 to 0 through the GNC stages

      this.max_iters       = 10;            % number of warping per pyramid level
      this.max_linear      = 1;             % maximum number of linearization performed per warping, 1 OK for HS
      
      % For GNC stage 1
      this.pyramid_levels  = 4;           
      this.pyramid_spacing = 2;

      % For GNC stage 2 to last
      this.gnc_pyramid_levels     = 2;
      this.gnc_pyramid_spacing    = 1.25;                
      
      method = 'lorentzian'; %'geman_mcclure'; %charbonnier
      this.spatial_filters = {[1 -1], [1; -1]};  
      for i = 1:length(this.spatial_filters);
          this.rho_spatial_u{i}   = robust_function(method, 0.03);
          this.rho_spatial_v{i}   = robust_function(method, 0.03);      
      end;
      this.rho_data        = robust_function(method, 1.5); % 6.3             
      
      this.alp             = 0.95;  % for ROF texture decomposition      
      
      this.color_images     = [];
      this.auto_level       = true;
      
      this = class(this, 'ba_optical_flow');         
      
    case 1
      if isa(varargin{1}, 'ba_optical_flow')
        this = varargin{1};        
      else    
          this = ba_optical_flow;
          this.images = varargin{1};  
      end
      
    otherwise
      error('Incompatible arguments!');
      
  end