function [data,units] = compute_dmax_wing_area(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  dareal = diff(trx(fly).wing_areal_mm);
  darear = diff(trx(fly).wing_arear_mm);
  data{i} = dareal;
  idx = trx(fly).wing_arear_mm(1:end-1) > trx(fly).wing_areal_mm(1:end-1);
  data{i}(idx) = darear(idx);
  
  data{i} = data{i} ./ trx(fly).dt;
    
end
units = parseunits('mm^2/s');

