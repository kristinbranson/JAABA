function r = subsref(this, selector)
%SUBSREF   Access parameters of ALT_BA optical flow algorithm
%   R = SUBSREF(THIS, SEL) extracts and returns parameter of a HS
%   optical algorithm THIS using given selector SEL.
%   Available selectors:
%     - ".optical_flow":   Flow base object
%     - ".pyramid_levels": Levels of the image pyramid
%     - ".warping_mode":   'forward' for forward warping and 'backward'
%                          for backward warping
%     - ".solver":         Selector for the linear equation solver
%                          ('backslash' for the builtin solver and 'sor'
%                          for a MEX-based SOR solver)
%     - ".max_iters":      Maximum number of iterations at each pyramid
%                          level
%     - ".limit_update":   Flag whether to restrict the incremental flow
%                          at each level to [-1, 1]
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

  
  % Allow full access to all fields of the structure
  r = builtin('subsref', this, selector);