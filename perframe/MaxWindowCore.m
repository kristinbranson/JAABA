function res = MaxWindowCore(x,r)

w = 2*r+1;
ntrans = size(x,1);
% pad with infs for boundary conditions
x_pad = [-inf(ntrans,r),x,-inf(ntrans,r)];
% maximum: use imdilate
se = strel(ones(1,w));
res = imdilate(x_pad,se);
