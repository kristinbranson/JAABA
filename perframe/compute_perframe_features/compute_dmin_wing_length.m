function [data,units] = compute_dmin_wing_length(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  dlengthl = diff(trx(fly).wing_lengthl_mm,1,2);
  dlengthr = diff(trx(fly).wing_lengthr_mm,1,2);
  data{i} = dlengthl;
  idx = trx(fly).wing_lengthr_mm(1:end-1) <= trx(fly).wing_lengthl_mm(1:end-1);
  data{i}(idx) = dlengthr(idx);
  
  data{i} = data{i} ./ trx(fly).dt;
  
end
units = parseunits('mm/s');

