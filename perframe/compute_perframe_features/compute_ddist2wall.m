% change in distance to wall
function [data,units] = compute_ddist2wall(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  data{i} = diff(trx(fly).dist2wall)./trx(fly).dt;
end
units = parseunits('mm/s');

