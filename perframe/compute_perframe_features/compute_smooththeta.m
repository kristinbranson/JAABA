% smoothed orientation
function [data,units] = compute_smooththeta(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  % smoothed orientation
  smooththeta = myconv(unwrap(trx(fly).theta_mm),trx.perframe_params.thetafil,'replicate','same');
  data{i} = modrange(smooththeta,-pi,pi);
end
units = parseunits('rad');

