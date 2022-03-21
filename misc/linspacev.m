function y = linspacev(d1, d2, n)
%LINSPACEV Linearly spaced vector.
%   Same as linspace, but works on vectors.
%   LINSPACEV(X1, X2) generates a matrix of size size(X1) x 100 of
%   equally spaced points between X1 and X2.
%
%   LINSPACEV(X1, X2, N) generates a matrix of size size(X1) x N 
%   equally points between X1 and X2.
%   For N = 1, LINSPACE returns X2.
%
%   Class support for inputs X1,X2:
%      float: double, single
%
% D1 and D2 must be the same size, N must be a scalar. 

if nargin == 2
    n = 100;
else
    n = floor(double(n));
end
if ~isscalar(n)
    error(message('MATLAB:linspacev:scalarInputs'));
end
if ~isequal(size(d1),size(d2)),
  error('Inputs d1 and d2 must be the same size');
end
n1 = n-1;
c = (d2 - d1).*(n1-1); %check intermediate value for appropriate treatment

sz = size(d1);
if numel(sz) == 2 && sz(2) == 1,
  sz = sz(1);
end

d1 = d1(:);
d2 = d2(:);
nel = numel(d1);

y = nan(nel,n);

cinf = isinf(c);
dinf = isinf(d2-d1);

if any(cinf&dinf),
  y(cinf&dinf,:) = d1(cinf&dinf) + (d2(cinf&dinf)./n1).*(0:n1) - (d1(cinf&dinf)./n1).*(0:n1);
end
if any(cinf&~dinf),
  y(cinf&~dinf,:) = d1(cinf&~dinf) + (0:n1).*((d2(cinf&~dinf) - d1(cinf&~dinf))./n1);
end
if any(~cinf),
  y(~cinf,:) = d1(~cinf) + (0:n1).*(d2(~cinf) - d1(~cinf))./n1;
end

deq = d1 == d2;
if any(deq),
  y(deq,:) = d1(deq);
end
if any(~deq),
  y(~deq,1) = d1;
  y(~deq,end) = d2;
end

y = reshape(y,[sz,n]);

