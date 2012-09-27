% distance to wall
function [data,units] = compute_dist2corner_rect(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  x = trx(fly).x;
  y = trx(fly).y;
  
  dtl = sqrt( (trx.tl_x(fly)-x).^2 + (trx.tl_y(fly)-y).^2);
  dtr = sqrt( (trx.tr_x(fly)-x).^2 + (trx.tr_y(fly)-y).^2);
  dbl = sqrt( (trx.bl_x(fly)-x).^2 + (trx.bl_y(fly)-y).^2);
  dbr = sqrt( (trx.br_x(fly)-x).^2 + (trx.br_y(fly)-y).^2);
  data{i} = min([dtl;dtr; dbl; dbr],[],1);

end
units = parseunits('mm');


