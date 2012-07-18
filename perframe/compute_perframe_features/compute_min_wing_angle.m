function [data,units] = compute_min_wing_angle(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  % signed minimum: left wing will be negative, right wing positive
  data{i} = min(-trx(fly).wing_anglel,trx(fly).wing_angler);
  
end
units = parseunits('rad');

