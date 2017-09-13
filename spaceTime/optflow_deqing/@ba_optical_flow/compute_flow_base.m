function uv = compute_flow_base(this, uv)
%%
%COMPUTE_FLOW_BASE   Base function for computing flow field
%   UV = COMPUTE_FLOW_BASE(THIS, INIT) computes the flow field UV with
%   algorithm THIS and the initialization INIT.
%  
%   This is a member function of the class 'ba_optical_flow'. 
%
% Authors: Deqing Sun, Department of Computer Science, Brown University
%          Stefan Roth, Department of Computer Science, TU Darmstadt
% Contact: dqsun@cs.brown.edu, sroth@cs.tu-darmstadt.de
% $Date: 2007-10-28$
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

  % Construct quadratic formulation
  qua_this          = this;
  qua_this.lambda   = this.lambda_q;  
  
  a = this.rho_spatial_u{1}.param;
  b = this.rho_data.param;    
  ta = b(1)/a(1);
  
  % Spatial
  if isa(this.rho_spatial_u{1}, 'robust_function')      
      for i = 1:length(this.rho_spatial_u)          
          qua_this.rho_spatial_u{i}   = robust_function('quadratic', 1);
          qua_this.rho_spatial_v{i}   = robust_function('quadratic', 1);
      end;     
  else
      error('evaluate_log_posterior: unknown rho function!');
  end;      
  
  % Data
  if isa(qua_this.rho_data, 'robust_function')      
      qua_this.rho_data        = robust_function('quadratic', ta);
  else
      error('evaluate_log_posterior: unknown rho function!');
  end;
  
 
  % Iterate flow computation
  for i = 1:this.max_iters
      
    duv = zeros(size(uv));   
        
    % Compute spatial and temporal partial derivatives
    [It Ix Iy] = partial_deriv(this.images, uv, this.interpolation_method, this.deriv_filter, this.blend);

    for j = 1:this.max_linear       
        
        % Every linearization step updates the nonlinearity using the
        % previous flow increments
        
        % Compute linear flow operator
        if this.alpha == 1
            [A, b, parm, iterative] = ...
                flow_operator(qua_this, uv, duv, It, Ix, Iy);        
            
        elseif this.alpha > 0
            [A, b] = ...
                flow_operator(qua_this, uv, duv, It, Ix, Iy);        
            [A1, b1, parm, iterative] = ...
                flow_operator(this, uv, duv, It, Ix, Iy);        
            A = this.alpha * A + (1-this.alpha) * A1;
            b = this.alpha * b + (1-this.alpha) * b1;

        elseif this.alpha == 0
            [A, b, parm, iterative] = ...
                flow_operator(this, uv, duv, It, Ix, Iy);        

        else
            error('flow_operator@ba_optical_flow: wrong gnc parameter alpha %3.2e', this.alpha);
        end;

        % Invoke the selected linear equation solver
        switch (lower(this.solver))
            case 'backslash'
                x = reshape(A \ b, size(uv));
            case 'sor'
                [x, flag, res, n] = sor(A', b, 1.9, this.sor_max_iters, 1E-2, uv(:));
                x = reshape(x, size(uv));
                fprintf('%d %d %d  ', flag, res, n);
            case 'bicgstab'
                [x,flag] = reshape(bicgstab(A, b, 1E-3, 200, [], [], uv(:)), size(uv)); %, parm
            case 'pcg'
                [x flag] = pcg(A,b, [], 10); 
                x        = reshape(x, size(uv));
            otherwise
                error('Invalid solver!')
        end

        % If limiting the incremental flow to [-1, 1] is requested, do so
        if (this.limit_update)
            x(x > 1)  = 1;
            x(x < -1) = -1;
        end
        
        % Print status information
        if this.display
            disp(['--Iteration: ', num2str(i), '   ', num2str(j), '   (', ...
                num2str(norm(x(:)-duv(:))), ')'])
        end;

        % Terminate iteration early if flow doesn't change substantially
%         if norm(x(:)-duv(:)) < 1E-3
%             break;
%         end
        
        duv = x;               
        
        uv0 = uv;
        
        uv  = uv+duv;        
        if ~isempty(this.median_filter_size)
            % Perform median filter
            uv(:,:,1) = medfilt2(uv(:,:,1), this.median_filter_size, 'symmetric');
            uv(:,:,2) = medfilt2(uv(:,:,2), this.median_filter_size, 'symmetric');
        end;
        duv = uv - uv0;
        
        uv  = uv0;

    end;   

    % Update flow fileds
    uv = uv + duv;
    
  end  
