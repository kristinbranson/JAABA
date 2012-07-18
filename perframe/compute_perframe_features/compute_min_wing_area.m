function [data,units] = compute_min_wing_area(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = min(trx(fly).wing_areal_mm,trx(fly).wing_arear_mm);
  
end
units = parseunits('mm^2');

