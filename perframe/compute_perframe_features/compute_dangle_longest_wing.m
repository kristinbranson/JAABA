function [data,units] = compute_dangle_longest_wing(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  danglel = diff(-trx(fly).wing_anglel,1,2);
  dangler = diff(trx(fly).wing_angler,1,2);
  data{i} = danglel;
  idx = trx(fly).wing_lengthr_mm(1:end-1) > trx(fly).wing_lengthl_mm(1:end-1);
  data{i}(idx) = dangler(idx);
  data{i} = data{i} ./ trx(fly).dt;
  
end
units = parseunits('rad/s');

