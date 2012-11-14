% distance to wall
function [data,units] = compute_dist2wall(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  
  data{i} = trx.landmark_params{n}.arena_radius_mm - trx(fly).arena_r;

end
units = parseunits('mm');

