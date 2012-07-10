% smoothed change in orientation
function [data,units] = compute_smoothdtheta(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  if trx(fly).nframes <= 1,
    data{i} = [];
  else
    smooththeta = myconv(unwrap(trx(fly).theta_mm),trx.perframe_params.thetafil,'replicate','same');
    data{i} = diff(smooththeta,1,2)./trx(fly).dt;
  end
end
units = parseunits('rad/s');

