% change in body orientation
function [data,units] = compute_dtheta(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = modrange(diff(trx(fly).theta_mm),-pi,pi)./trx(fly).dt;
end
units = parseunits('rad/s');

