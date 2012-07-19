function [data,units] = compute_mean_wing_angle(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  % signed mean: left wing will be negative, right wing positive
  data{i} = (-trx(fly).wing_anglel+trx(fly).wing_angler)/2;
  
end
units = parseunits('rad');

