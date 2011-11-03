% absolute value of smoothed change in orientation
function [data,units] = compute_abssmoothdtheta(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = abs(trx(fly).smoothdtheta);
end
units = parseunits('rad/s');

