function res = MinWindowCore(x,r)

if isempty(x)
  res = zeros(size(x));
  return;
end
w = 2*r+1;
ntrans = size(x,1);
% pad with infs for boundary conditions
x(isnan(x)) = inf;
x_pad = [-inf(ntrans,r),x,-inf(ntrans,r)];
% maximum: use imdilate
se = strel(ones(1,w));
res = imerode(x_pad,se);
