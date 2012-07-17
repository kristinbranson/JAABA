function [data,units] = compute_dwing_angle_diff(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = diff(trx(fly).wing_angle_diff) ./ trx(fly).dt;
  
end
units = parseunits('rad/s');

