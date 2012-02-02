function y = weighted_prctile(x,p,w,dim,issorted)
% WEIGHTED_PRCTILE Percentiles of a sample.
%   Y = WEIGHTED_PRCTILE(X,P,W) returns percentiles of the values in X.  P is a scalar
%   or a vector of percent values.  When X is a vector, Y is the same size
%   as P, and Y(i) contains the P(i)-th percentile.  When X is a matrix,
%   the i-th row of Y contains the P(i)-th percentiles of each column of X.
%   For N-D arrays, WEIGHTED_PRCTILE operates along the first non-singleton
%   dimension.
%
%   Y = WEIGHTED_PRCTILE(X,P,W,DIM) calculates percentiles along dimension DIM.  The
%   DIM'th dimension of Y has length LENGTH(P).
%
%   Percentiles are specified using percentages, from 0 to 100.  For an N
%   element vector X, WEIGHTED_PRCTILE computes percentiles as follows:
%      1) The sorted values in X are taken as the 100*(0.5/N), 100*(1.5/N),
%         ..., 100*((N-0.5)/N) percentiles.
%      2) Linear interpolation is used to compute percentiles for percent
%         values between 100*(0.5/N) and 100*((N-0.5)/N)
%      3) The minimum or maximum values in X are assigned to percentiles
%         for percent values outside that range.
%
%   WEIGHTED_PRCTILE treats NaNs as missing values, and removes them.
%
%   Examples:
%      y = weighted_prctile(x,50); % the median of x
%      y = weighted_prctile(x,[2.5 25 50 75 97.5]); % a useful summary of x
%
%   See also IQR, MEDIAN, NANMEDIAN, QUANTILE.

%   Copyright 1993-2004 The MathWorks, Inc.
%   $Revision: 2.12.4.5 $  $Date: 2004/06/25 18:52:56 $

if ~isvector(p) || numel(p) == 0
    error('stats:prctile:BadPercents', ...
          'P must be a scalar or a non-empty vector.');
elseif any(p < 0 | p > 100) || ~isreal(p)
    error('stats:prctile:BadPercents', ...
          'P must take real values between 0 and 100');
end

% make sure sizes match
sz = size(x);
szw = size(w);
sz0 = sz;
if numel(x) ~= numel(w),
  error('Number of elements of x and w must match');
end
if length(sz) ~= length(szw) || ~all(sz == szw),
  if ~(isvector(x) && isvector(w)),
    warning('Number of elements of x and w match, but not sizes. Resizing w to size %s (was size %s)',mat2str(sz),mat2str(szw));
  end
end

% Figure out which dimension prctile will work along.
if nargin < 5,
  issorted = false;
end
if nargin < 4 
    dim = find(sz ~= 1,1);
    if isempty(dim)
        dim = 1; 
    end
    dimArgGiven = false;
else
    % Permute the array so that the requested dimension is the first dim.
    nDimsX = ndims(x);
    perm = [dim:max(nDimsX,dim) 1:dim-1];
    x = permute(x,perm);
    w = permute(w,perm);
    % Pad with ones if dim > ndims.
    if dim > nDimsX
        sz = [sz ones(1,dim-nDimsX)];
    end
    sz = sz(perm);
    dim0 = dim;
    dim = 1;
    dimArgGiven = true;
end

% If X is empty, return all NaNs.
if isempty(x)
    if isequal(x,[]) && ~dimArgGiven
        y = nan(size(p),class(x));
    else
        szout = sz; szout(dim) = numel(p);
        y = nan(szout,class(x));
    end

else
    % Drop X's leading singleton dims, and combine its trailing dims.  This
    % leaves a matrix, and we can work along columns.
    nrows = sz(dim);
    ncols = prod(sz) / nrows;
    x = reshape(x, nrows, ncols);
    w = reshape(w, nrows, ncols);
    y = zeros(numel(p), ncols, class(x));
    
    if ~issorted,
      % sort x
      [x,order] = sort(x,1);
      % apply ordering to w
      w = w(bsxfun(@plus,order,0:nrows:(ncols-1)*nrows));
      %w = w(sub2ind([nrows,ncols],order,repmat((1:ncols),[nrows,1])));
    end
    nans = isnan(x);
    % set weight to 0 for nans
    w(nans) = 0;
    
    % normalize w to sum to 1
    %w = w ./ repmat(sum(w,1),[nrows,1]);
    w = bsxfun(@rdivide,w,sum(w,1));
    
    % compute cumulative sum
    cw = cumsum(w,1);
    % for numerical errors
    cw = cw / cw(end);
    
    for i = 1:numel(p),
      % find the frame where cumsum of weight is at least p
      isbigger = [false(1,ncols);cw >= p(i)/100];
      ind = find(isbigger(1:end-1,:)==false & isbigger(2:end,:)==true);
      % check for border
      isexact = cw(ind) == p(i)/100 && ind < numel(cw);
      if isexact,
        y(i) = (x(ind) + x(ind+1))/2;
      else
        % store this index
        y(i) = x(ind);
      end
    end

end

% undo the DIM permutation
if dimArgGiven
     y = ipermute(y,perm);  
end

% If X is a vector, the shape of Y should follow that of P, unless an
% explicit DIM arg was given.
if ~dimArgGiven && isvector(x)
    y = reshape(y,size(p)); 
else
  sz2 = sz0;
  sz2(dim0) = 1;
  y = reshape(y,sz2);
end
