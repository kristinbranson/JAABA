function n = psnr(x, y)
%PSNR   Compute the peak signal to noise ratio between two images
%   PSNR(X, Y) computes the peak signal to noise ratio in decibels (dB)
%   between images X and Y.
%   The images are assumed to have 256 gray levels (0 .. 255).
%
%   Class support for inputs X and Y:
%        float:   double, single
%        integer: uint*, int*
%
%   See also SNR, IMMSE.
%  
%   Author:  Stefan Roth, Department of Computer Science, Brown University
%   Contact: roth@cs.brown.edu
%   $Date: 2005-06-08 17:10:38 -0400 (Wed, 08 Jun 2005) $
%   $Revision: 48 $
%
% Copyright 2005-2007, Brown University, Providence, RI. USA
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
  
  diff = double(x) - double(y);
  n = 20 * log10(255 / std(diff(:)));

