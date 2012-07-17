function [data,units] = compute_wing_angle_diff(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  % difference: left wing will be negative, right wing positive
  data{i} = trx(fly).wing_angler - trx(fly).wing_anglel;
  
end
units = parseunits('rad');

