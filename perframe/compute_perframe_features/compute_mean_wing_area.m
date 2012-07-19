function [data,units] = compute_mean_wing_area(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = (trx(fly).wing_areal_mm+trx(fly).wing_arear_mm)/2;
  
end
units = parseunits('mm^2');

