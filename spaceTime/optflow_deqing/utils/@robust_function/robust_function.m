function this = robust_function(varargin)
%ROBUST_FUNCTION   Construct a robust function
%   ROBUST_FUNCTION(TYPE[, PARAM]) constructs a robust function of type
%   TYPE with optional parameters PARAM.
%   The following robust function types are available:
%    - 'geman_mcclure'
%    - 'huber'
%    - 'lorentzian'
%    - 'quadratic'
%    - 'tukey'
%    - 'spline'
%    - 'mixture'
%    - 'gaussian' (like quadratic, but normalized)
%    - 'tdist' (normalized)
%    - 'tdist_unnorm' (unnormalized)
%    - 'charbonnier'
%
%   ROBUST_FUNCTION(O) constructs a robust function by copying O.
%  
%   This is a member function of the class 'robust_function'. 
%
%   Author:  Stefan Roth, Department of Computer Science, TU Darmstadt
%   Contact: sroth@cs.tu-darmstadt.de
%   $Date $
%   $Revision$

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
   

  
  error(nargchk(0, 3, length(varargin)));
  
  switch (length(varargin))
   case 0
    this.type  = @quadratic;
    this.param = 1;
    this = class(this, 'robust_function');
    
   case 1
    if isa(varargin{1}, 'robust_function')
      this = other;
    else
      switch (varargin{1})
       case 'geman_mcclure'
        this.type = @geman_mcclure;
       case 'huber'
        this.type = @huber;
       case 'lorentzian'
        this.type = @lorentzian;
       case 'quadratic'
        this.type = @quadratic;
       case 'tukey'
        this.type = @tukey;
       case 'spline'
        this.type = @spline;
       case 'mixture'
        this.type = @mixture;
       case 'gaussian'
        this.type = @gaussian;
       case 'tdist'
        this.type = @tdist;
       case 'tdist_unnorm'
        this.type = @tdist_unnorm;
       case 'charbonnier'
        this.type = @charbonnier;
       case 'generalized_charbonnier'
        this.type = @generalized_charbonnier;
       otherwise
        error('Invalid robust function type.');
      end

      this.param = 1;
      
      if strcmp(varargin{1}, 'generalized_charbonnier')
          this.param = [1e-3, 1];
      end;
      
      
      this = class(this, 'robust_function');
    end

   case 2
    switch (varargin{1})
     case 'geman_mcclure'
      this.type = @geman_mcclure;
     case 'huber'
      this.type = @huber;
     case 'lorentzian'
      this.type = @lorentzian;
     case 'quadratic'
      this.type = @quadratic;
     case 'tukey'
      this.type = @tukey;
     case 'spline'
      this.type = @spline;
     case 'mixture'
      this.type = @mixture;
     case 'gaussian'
      this.type = @gaussian;
     case 'tdist'
      this.type = @tdist;
     case 'tdist_unnorm'
      this.type = @tdist_unnorm;
     case 'charbonnier'
      this.type = @charbonnier;
     otherwise
      error('Invalid robust function type.');
    end

    this.param = varargin{2};
    this = class(this, 'robust_function');
  
      case 3
          switch (varargin{1})
              case 'generalized_charbonnier'
                  this.type = @generalized_charbonnier;
                  
                  this.param = [varargin{2} varargin{3}];
                  this = class(this, 'robust_function');                  
                  
              otherwise
                  error('Invalid robust function type.');
          end;
  end