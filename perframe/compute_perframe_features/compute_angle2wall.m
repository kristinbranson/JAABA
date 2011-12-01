% angle to closest point on the wall in the fly's coordinate system
function [data,units] = compute_angle2wall(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  
  data{i} = modrange(trx(fly).arena_angle-trx(fly).theta_mm,-pi,pi);
    
end
units = parseunits('rad');

