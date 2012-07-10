% change in polar angle to closest point on wall
function [data,units] = compute_darena_angle(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  data{i} = modrange(diff(trx(fly).arena_angle,1,2),-pi,pi)./trx(fly).dt;
end
units = parseunits('rad/s');

