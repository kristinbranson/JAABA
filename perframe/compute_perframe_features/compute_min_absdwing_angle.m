function [data,units] = compute_min_absdwing_angle(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = min(abs(diff(trx(fly).wing_anglel)),abs(diff(trx(fly).wing_angler))) ./ trx(fly).dt;
  
end
units = parseunits('rad/s');

