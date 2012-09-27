% distance to wall
function [data,units] = compute_angle2corner_rect(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  x = trx(fly).x;
  y = trx(fly).y;
  
  cornera = zeros(1,4);
  cornera(1) = atan2(trx.tl_y(fly)-trx.br_y(fly),trx.tl_x(fly)-trx.br_x(fly));
  cornera(2) = atan2(trx.tr_y(fly)-trx.bl_y(fly),trx.tr_x(fly)-trx.bl_x(fly));
  cornera(3) = atan2(trx.br_y(fly)-trx.tl_y(fly),trx.br_x(fly)-trx.tl_x(fly));
  cornera(4) = atan2(trx.bl_y(fly)-trx.tr_y(fly),trx.bl_x(fly)-trx.tr_x(fly));
  
  dtl = sqrt( (trx.tl_x(fly)-x).^2 + (trx.tl_y(fly)-y).^2);
  dtr = sqrt( (trx.tr_x(fly)-x).^2 + (trx.tr_y(fly)-y).^2);
  dbr = sqrt( (trx.br_x(fly)-x).^2 + (trx.br_y(fly)-y).^2);
  dbl = sqrt( (trx.bl_x(fly)-x).^2 + (trx.bl_y(fly)-y).^2);

  
  [~,closest] = min([dtl;dtr; dbr; dbl],[],1);
  
  closest_cornera = cornera(closest);
  targeta = trx(fly).theta;
  
  a = mod(targeta-closest_cornera+pi,2*pi)-pi;
  
  data{i} = a;
end
units = parseunits('rad');


