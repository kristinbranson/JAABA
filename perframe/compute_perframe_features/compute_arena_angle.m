% polar angle of closest point on the wall
function [data,units] = compute_arena_angle(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  
  data{i} = atan2(trx(fly).y_mm-trx.landmark_params{n}.arena_center_mm_y(n),...
    trx(fly).x_mm-trx.landmark_params{n}.arena_center_mm_x(n));
    
end
units = parseunits('rad');

