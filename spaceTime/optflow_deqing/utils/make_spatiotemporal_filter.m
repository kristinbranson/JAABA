function f = make_spatiotemporal_filter(type, direction, grid_spacing)
%MAKE_SPATIOTEMPORAL_FILTER   A number of spatiotemporal derivative filters
%   MAKE_SPATIOTEMPORAL_FILTER(TYPE, DIR, SPACING) creates a linear
%   filter of type TYPE.  DIR specifies the derivative axis (1 - X, 2 -
%   Y, 3 - Z), and SPACING is a 3x1 vector or grid spacings along the
%   various axes.  The following filter types are supported:
%    - 'central':      Central differences (3x1x1), no cross-smoothing.
%    - 'forward':      Forward differences (2x1x1), no cross-smoothing.
%    - 'horn':         Neighbor differences (2x2x2), cross-smoothing
%                      (from B.K.P. Horn's book).
%    - 'scharr_3x3x2': Optimized filters by H. Scharr, cross-smoothing
%    - 'scharr_4x4x2': Optimized filters by H. Scharr, cross-smoothing
%
%   See also FSPECIAL.
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


  switch (lower(type))
    case 'central'
      f = cat(direction, 1, 0, -1) / (2  *  grid_spacing(direction));
      
    case 'weickert'
      f = cat(direction, 1, -9, 45, 0, -45, 9, -1) / (60  *  grid_spacing(direction));
    
    case 'horn'
      switch(direction)
        case 1
          f = repmat([1; -1], [1 2 2]) / (4  *  grid_spacing(direction));
        case 2
          f = repmat([1, -1], [2 1 2]) / (4  *  grid_spacing(direction));
        case 3
          f = repmat(cat(3, 1, -1), [2 2 1]) / (4  *  grid_spacing(direction));
        otherwise
          error('Invalid direction!');
      end
    
    case {'forward', 'backward'}
      f = cat(direction, 1, -1) / grid_spacing(direction);
    
    case 'scharr_3x3x2'
      switch(direction)
        case 1
          f = [1; 0; -1] * [0.2202 0.5597 0.2202] / (2  *  grid_spacing(direction));
          f = cat(3, 0.5 * f, 0.5 * f);
        case 2
          f = [0.2202; 0.5597; 0.2202] * [1 0 -1] / (2  *  grid_spacing(direction));
          f = cat(3, 0.5 * f, 0.5 * f);
        case 3
          f = [0.2202; 0.5597; 0.2202] * [0.2202 0.5597 0.2202] / grid_spacing(direction);
          f = cat(3, f, -f);
        otherwise
          error('Invalid direction!');
      end
     
    case 'scharr_4x4x2'
      switch(direction)
        case 1
          f = [0.2619; 0.2143; -0.2143; -0.2619] * [0.1166 0.3834 0.3834 0.1166] / grid_spacing(direction);
          f = cat(3, 0.5 * f, 0.5 * f);
        case 2
          f = [0.1166; 0.3834; 0.3834; 0.1166] * [0.2619 0.2143 -0.2143 -0.2619] / grid_spacing(direction);
          f = cat(3, 0.5 * f, 0.5 * f);
        case 3
          f = [0.1166; 0.3834; 0.3834; 0.1166] * [0.1166 0.3834 0.3834 0.1166] / grid_spacing(direction);
          f = cat(3, f, -f);
        otherwise
          error('Invalid direction!');
      end
     
    otherwise
      error('Unknown filter type!')
      
  end
