% magnitude of velocity of center
function [data,units] = compute_velmag_ctr(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  % change in center position
  dx = diff(trx(fly).x_mm,1,2);
  dy = diff(trx(fly).y_mm,1,2);
  
  if trx(fly).nframes < 2,
    data{i} = [];
  else
    % magnitude of velocity vector
    data{i} = sqrt(dx.^2 + dy.^2)./trx(fly).dt;
  end
end
units = parseunits('mm/s');

