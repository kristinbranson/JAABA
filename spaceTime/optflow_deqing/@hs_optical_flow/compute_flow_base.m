function uv = compute_flow_base(this, uv)

%COMPUTE_FLOW_BASE   Base function for computing flow field
%   UV = COMPUTE_FLOW_BASE(THIS, INIT) computes the flow field UV with
%   algorithm THIS and the initialization INIT.
%  
%   This is a member function of the class 'hs_optical_flow'. 
%
%   Author: Deqing Sun, Department of Computer Science, Brown University
%   Contact: dqsun@cs.brown.edu
%   $Date: 2007-11-30 $
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

  % Iterate flow computation
  for i = 1:this.max_warping_iters             
   
    % Compute linear flow operator
    [A, b, parm, iterative] = ...
        flow_operator(this, uv);
    
    % Invoke the selected linear equation solver
    switch (lower(this.solver))
      case 'backslash'
        x = reshape(A \ b, size(uv));
      case 'sor'
        [x, flag, res, n] = sor(A', b, 1.9, this.sor_max_iters, 1E-2, uv(:));
        x = reshape(x, size(uv));
        fprintf('%d %d %d  ', flag, res, n);
      case 'bicgstab'
        x = reshape(bicgstab(A, b, 1E-3, 200, [], [], uv(:), parm), size(uv));
      case 'pcg'
          [x flag] = pcg(A,b, [], 100);  %100           
      otherwise
        error('Invalid solver!')
    end
    
    % Print status information
    if this.display
        disp(['--Iteration: ', num2str(i), '    (', ...
             num2str(norm(x(:))), ')'])
    end;
    
    % Terminate iteration early if flow doesn't change substantially
    if (length(this.lambda) == 1 && norm(x(:)) < 1E-3)
      break
    end
   
    % If limiting the incremental flow to [-1, 1] is requested, do so
    if (this.limit_update)
      x(x > 1)  = 1;
      x(x < -1) = -1;
    end
    
    uv = uv + x;
        
    % Perform median filtering to remove outliers
    if ~isempty(this.median_filter_size)
        for m = 1:this.mf_iter; % extensive MF filtering 
            uv(:,:,1) = medfilt2(uv(:,:,1), this.median_filter_size, 'symmetric');
            uv(:,:,2) = medfilt2(uv(:,:,2), this.median_filter_size, 'symmetric');
        end;
    end;
    
    
    % Terminate early if the flow_operator doesn't require multiple
    % interations 
%     if (~iterative)
%       break;
%     end
    
  end
