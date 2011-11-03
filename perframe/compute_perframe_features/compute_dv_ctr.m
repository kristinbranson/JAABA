% sideways motion of body center
function [data,units] = compute_dv_ctr(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  % change in center position
  dx = diff(trx(fly).x_mm);
  dy = diff(trx(fly).y_mm);
  
  if trx(fly).nframes < 2,
    data{i} = [];
  else
    data{i} = (dx.*cos(trx(fly).theta_mm(1:end-1)+pi/2) + dy.*sin(trx(fly).theta_mm(1:end-1)+pi/2))./trx(fly).dt;
  end
end
units = parseunits('mm/s');

