% absolute change in angle to closest point on wall in fly's coordinate system
function [data,units] = compute_absdangle2wall(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  data{i} = abs(trx(fly).dangle2wall);
end
units = parseunits('rad/s');

