function ope = load_of_method(method)
%%
%LOAD_OF_METHOD loads the models for Optical Flow Estimation
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

median_filter_size = [5 5]; 

switch method

    case 'classic+nl-fast'
        ope = load_of_method('classic+nl');        
        ope.max_iters       = 3;
        ope.gnc_iters       = 2;
        ope.display  = true;
        
    case 'classic+nl'        
        ope = classic_nl_optical_flow;
        
        ope.texture  = true;
        ope.median_filter_size   = median_filter_size;
        
        ope.alp = 0.95;
        ope.area_hsz = 7;
        ope.sigma_i  = 7;
        ope.color_images = ones(1,1,3);
        
        ope.lambda = 3;
        ope.lambda_q =3;
        %ope.display  = true;        
        
    case 'classic+nl-full'
        ope = load_of_method('classic+nl');
        ope.fullVersion      = true;
       
    case {'hs-brightness'}
        % brightness constancy, with median filtering
        ope = hs_optical_flow;
        ope.median_filter_size   = median_filter_size;
        ope.lambda   = 10;
        ope.lambda_q = 10;        
        
    case {'hs'}
        % ROF texture constancy, with median filtering
        ope = hs_optical_flow;
        ope.median_filter_size   = median_filter_size;
        ope.texture              = true;        

        ope.lambda   = 40;
        ope.lambda_q = 40;
        ope.display  = true;       

    case {'ba-brightness'}        
        % B&A brightness constancy model
        
        ope = ba_optical_flow;
        ope.median_filter_size   = median_filter_size;
        
        method = 'lorentzian'; %'geman_mcclure'; 
        ope.spatial_filters = {[1 -1], [1; -1]};
        for i = 1:length(ope.spatial_filters);
            ope.rho_spatial_u{i}   = robust_function(method, 0.1); 
            ope.rho_spatial_v{i}   = robust_function(method, 0.1);
        end;
        ope.rho_data        = robust_function(method, 3.5);
        
        ope.lambda = 0.045;
        ope.lambda_q = 0.045 ;        
       
    case {'classic-l', 'ba'}
        
        ope = load_of_method('ba-brightness');               
        ope.median_filter_size   = median_filter_size;
        ope.texture              = true;        
       
        method = 'lorentzian'; %'geman_mcclure'; 
        ope.spatial_filters = {[1 -1], [1; -1]};
        for i = 1:length(ope.spatial_filters);
            ope.rho_spatial_u{i}   = robust_function(method, 0.03); 
            ope.rho_spatial_v{i}   = robust_function(method, 0.03);
        end;
        ope.rho_data        = robust_function(method, 1.5);        
      
        ope.lambda = 0.06;
        ope.lambda_q = 0.06 ;
        
    case 'classic-c-a'
        ope = alt_ba_optical_flow;            
        ope.median_filter_size   = median_filter_size;
        ope.texture              = true;
        method = 'charbonnier'; 
        ope.spatial_filters = {[1 -1], [1; -1]};
        for i = 1:length(ope.spatial_filters);
            ope.rho_spatial_u{i}   = robust_function(method, 1e-3);
            ope.rho_spatial_v{i}   = robust_function(method, 1e-3);
        end;
        ope.rho_data        = robust_function(method, 1e-3);
        ope.display         = false;
        
        ope.lambda2         = 1e2;    
        ope.lambda3         = 1;      
        
        ope.weightRatio     = ope.lambda2/ope.lambda3;     % lambda2/lambda3
        ope.itersLO         = 5;     % # Li & Osher median formula iteration
        
        ope.lambda   = 5;
        ope.lambda_q = 5;
        
    case {'classic-c-brightness'}
        %charbonnier penalty function
        ope = ba_optical_flow;
        ope.median_filter_size   = median_filter_size;
        ope.texture              = false;        
        method = 'charbonnier'; 
        ope.spatial_filters = {[1 -1], [1; -1]};
        for i = 1:length(ope.spatial_filters);
            ope.rho_spatial_u{i}   = robust_function(method, 1e-3); 
            ope.rho_spatial_v{i}   = robust_function(method, 1e-3);
        end;
        ope.rho_data        = robust_function(method, 1e-3); 
        
        ope.lambda   = 3;
        ope.lambda_q = 3;
        
    case {'classic-c'}
        ope = load_of_method('classic-c-brightness');
        ope.texture            = true;        
        ope.lambda   = 5;
        ope.lambda_q = 5;
        
    case {'classic++'}
        ope = ba_optical_flow;
        ope.median_filter_size = median_filter_size;
        ope.texture            = true;
        ope.interpolation_method = 'bi-cubic';        
        
        method = 'generalized_charbonnier'; 
        ope.spatial_filters = {[1 -1], [1; -1]};
        a   = 0.45;
        sig = 1e-3;
        for i = 1:length(ope.spatial_filters);
            ope.rho_spatial_u{i}   = robust_function(method, sig, a);
            ope.rho_spatial_v{i}   = robust_function(method, sig, a);
        end;
        ope.rho_data        = robust_function(method, sig, a);                
        
        ope.lambda = 3;
        ope.lambda_q =3;
        %ope.display = true;
        
    otherwise
        error('unknown optical flow estimation method!');
end;
