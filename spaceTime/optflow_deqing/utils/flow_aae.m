function [aae, sae] = flow_aae(f1, f2, mask)
%FLOW_AAE   Average angular error in optical flow computation
%   AAE = FLOW_AAE(F1, F2[, MASK]) computes the average angular error in
%   degrees between the flow fields F1 and F2.  The optional argument
%   MASK can specify which of the pixels should be taken into account.
%  
%   [AAE, SAE] = FLOW_AAE(F1, F2[, MASK]) also returns the standard
%   deviation of the angular error.
%
%   Author:  Stefan Roth, Department of Computer Science, TU Darmstadt
%   Contact: sroth@cs.tu-darmstadt.de
%   $Date: 2007-03-27 14:09:11 -0400 (Tue, 27 Mar 2007) $
%   $Revision: 252 $

% Copyright 2004-2007, Brown University, Providence, RI. USA
% Copyright 2007-2010 TU Darmstadt, Darmstadt, Germany.
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
  
  
  aae = acos((sum(f1 .* f2, 3) + 1) ./ ...
             sqrt((sum(f1 .^ 2, 3) + 1) .* (sum(f2 .^ 2, 3) + 1)));
  if (nargin > 2)
    aae = aae(mask);
  end
    
  sae = std(real(aae(:))) * (180 / pi);
  aae = mean(real(aae(:))) * (180 / pi);

