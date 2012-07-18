function [data,units] = compute_dmean_wing_area(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = diff(trx(fly).mean_wing_area) ./ trx(fly).dt;
  
end
units = parseunits('mm^2/s');

