function [data,units] = compute_dmean_wing_angle(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = diff(trx(fly).mean_wing_angle) ./ trx(fly).dt;
  
end
units = parseunits('rad/s');

