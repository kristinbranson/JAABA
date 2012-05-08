% absolute change in polar angle to closest point on wall
function [data,units] = compute_absdarena_angle(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  data{i} = abs(trx(fly).darena_angle);
end
units = parseunits('rad/s');

