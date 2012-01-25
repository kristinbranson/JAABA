% magnitude of velocity of tail
function [data,units] = compute_velmag_tail(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  % location of tail
  tailx = trx(fly).x_mm + 2*cos(-trx(fly).theta).*trx(fly).a_mm;
  taily = trx(fly).y_mm + 2*sin(-trx(fly).theta).*trx(fly).a_mm;
  dx = diff(tailx);
  dy = diff(taily);
  
  if trx(fly).nframes < 2,
    data{i} = [];
  else
    % magnitude of velocity vector
    data{i} = sqrt(dx.^2 + dy.^2)./trx(fly).dt;
  end
end
units = parseunits('mm/s');

