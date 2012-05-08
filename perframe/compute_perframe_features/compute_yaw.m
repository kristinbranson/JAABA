% difference between orientation and velocity direction
function [data,units] = compute_yaw(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  data{i} = modrange(trx(fly).phi - trx(fly).theta_mm,-pi,pi);
end
units = parseunits('rad');

