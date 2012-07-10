% change in body orientation
function [data,units] = compute_dtheta(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  if trx(fly).nframes <= 1,
    data{i} = [];
  else
    data{i} = modrange(diff(trx(fly).theta_mm,1,2),-pi,pi)./trx(fly).dt;
  end
end
units = parseunits('rad/s');

