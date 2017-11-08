function mask = compute_conv2_valid_mask(imsz, fsz)

% compute a mask of size IMSZ with 1 entry at valid pixel positions and 0
% invalid oens
%
%   Author:  Stefan Roth, Department of Computer Science, TU Darmstadt
%   Contact: sroth@cs.tu-darmstadt.de
%   $Date:  $
%   $Revision: $

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

mask = ones(imsz);

h1   = floor(fsz(1)/2);
if fsz(1) == 1
    ;
elseif fsz(1) == 2
    mask(end,:) = 0;
elseif mod(fsz(1), 2) == 0
    mask(1:h1-1, :)         = 0;
    mask(end-h1+1:end, :)   = 0;
else
    mask(1:h1, :)    = 0;
    mask(end-h1+1:end, :)   = 0;    
end;

h2   = floor(fsz(2)/2);

if fsz(2) == 1
    ;
elseif fsz(2) == 2
    mask(:, end) = 0;
elseif mod(fsz(2), 2) == 0
    mask(:, 1:h2-1)    = 0;
    mask(:, end-h2+1:end)   = 0;    
else
    mask(:, 1:h2)    = 0;
    mask(:, end-h2+1:end)   = 0;        
end;