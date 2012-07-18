function [data,units] = compute_darea_outmost_wing(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  dareal = diff(trx(fly).wing_areal_mm);
  darear = diff(trx(fly).wing_arear_mm);
  data{i} = dareal;
  idx = trx(fly).wing_angler(1:end-1) > -trx(fly).wing_anglel(1:end-1);
  data{i}(idx) = darear(idx);
  
  data{i} = data{i} ./ trx(fly).dt;  
  
end
units = parseunits('mm^2/s');

