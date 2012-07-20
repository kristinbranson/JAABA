function [data,units] = compute_wing_angle_imbalance(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  % difference: left wing will be negative, right wing positive
  data{i} = abs(trx(fly).wing_angler+trx(fly).wing_anglel);
  % same as:
%   max(trx(fly).wing_angler+trx(fly).wing_anglel,...
%     -trx(fly).wing_anglel-trx(fly).wing_angler)
  
end
units = parseunits('rad');

