% absolute difference between orientation and velocity direction
function [data,units] = compute_absyaw(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  data{i} = abs(trx(fly).yaw);
end
units = parseunits('rad');

