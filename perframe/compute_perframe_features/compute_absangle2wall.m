% absolute angle to closest point on the wall in the fly's coordinate system
function [data,units] = compute_absangle2wall(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  
  data{i} = abs(trx(fly).angle2wall);
    
end
units = parseunits('rad');

