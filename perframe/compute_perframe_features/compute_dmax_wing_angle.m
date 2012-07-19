function [data,units] = compute_dmax_wing_angle(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  danglel = diff(-trx(fly).wing_anglel);
  dangler = diff(trx(fly).wing_angler);
  data{i} = danglel;
  idx = trx(fly).wing_angler(1:end-1) > trx(fly).wing_anglel(1:end-1);
  data{i}(idx) = dangler(idx);
  
  data{i} = data{i} ./ trx(fly).dt;
  
end
units = parseunits('rad/s');

