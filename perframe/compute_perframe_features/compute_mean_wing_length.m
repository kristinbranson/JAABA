function [data,units] = compute_mean_wing_length(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = (trx(fly).wing_lengthl_mm+trx(fly).wing_lengthr_mm)/2;
  
end
units = parseunits('mm');

