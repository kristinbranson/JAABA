% distance to wall
function [data,units] = compute_dist2corner_rect(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  x = trx(fly).x_mm;
  y = trx(fly).y_mm;
  
  dtl = sqrt( (trx.landmark_params{n}.tl_x(fly)-x).^2 + (trx.landmark_params{n}.tl_y(fly)-y).^2);
  dtr = sqrt( (trx.landmark_params{n}.tr_x(fly)-x).^2 + (trx.landmark_params{n}.tr_y(fly)-y).^2);
  dbl = sqrt( (trx.landmark_params{n}.bl_x(fly)-x).^2 + (trx.landmark_params{n}.bl_y(fly)-y).^2);
  dbr = sqrt( (trx.landmark_params{n}.br_x(fly)-x).^2 + (trx.landmark_params{n}.br_y(fly)-y).^2);
  data{i} = min([dtl;dtr; dbl; dbr],[],1);

end
units = parseunits('mm');


