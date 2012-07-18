function [data,units] = compute_angle_biggest_wing(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  % left is negative
  data{i} = -trx(fly).wing_anglel;
  idx = trx(fly).wing_arear_mm > trx(fly).wing_areal_mm;
  data{i}(idx) = trx(fly).wing_angler(idx);
  
end
units = parseunits('rad');

