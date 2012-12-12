% distance to wall
function [data,units] = compute_angle2corner_rect(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  x = trx(fly).x_mm;
  y = trx(fly).y_mm;
  
  cornera = zeros(1,4);
  cornera(1) = atan2(trx.landmark_params{n}.tl_y(fly)-trx.landmark_params{n}.br_y(fly),trx.landmark_params{n}.tl_x(fly)-trx.landmark_params{n}.br_x(fly));
  cornera(2) = atan2(trx.landmark_params{n}.tr_y(fly)-trx.landmark_params{n}.bl_y(fly),trx.landmark_params{n}.tr_x(fly)-trx.landmark_params{n}.bl_x(fly));
  cornera(3) = atan2(trx.landmark_params{n}.br_y(fly)-trx.landmark_params{n}.tl_y(fly),trx.landmark_params{n}.br_x(fly)-trx.landmark_params{n}.tl_x(fly));
  cornera(4) = atan2(trx.landmark_params{n}.bl_y(fly)-trx.landmark_params{n}.tr_y(fly),trx.landmark_params{n}.bl_x(fly)-trx.landmark_params{n}.tr_x(fly));
  
  dtl = sqrt( (trx.landmark_params{n}.tl_x(fly)-x).^2 + (trx.landmark_params{n}.tl_y(fly)-y).^2);
  dtr = sqrt( (trx.landmark_params{n}.tr_x(fly)-x).^2 + (trx.landmark_params{n}.tr_y(fly)-y).^2);
  dbr = sqrt( (trx.landmark_params{n}.br_x(fly)-x).^2 + (trx.landmark_params{n}.br_y(fly)-y).^2);
  dbl = sqrt( (trx.landmark_params{n}.bl_x(fly)-x).^2 + (trx.landmark_params{n}.bl_y(fly)-y).^2);

  
  [~,closest] = min([dtl;dtr; dbr; dbl],[],1);
  
  closest_cornera = cornera(closest);
  targeta = trx(fly).theta;
  
  a = mod(targeta-closest_cornera+pi,2*pi)-pi;
  
  data{i} = a;
end
units = parseunits('rad');


