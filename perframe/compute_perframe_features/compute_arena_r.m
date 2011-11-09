% distance from arena center
function [data,units] = compute_arena_r(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  
  data{i} = sqrt((trx(fly).x_mm - trx.landmark_params.arena_center_mm_x).^2 + ...
    (trx(fly).y_mm - trx.landmark_params.arena_center_mm_y).^2);

end
units = parseunits('mm');

