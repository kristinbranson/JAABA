function [data,units] = compute_length_outmost_wing(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = trx(fly).wing_lengthl_mm;
  idx = trx(fly).wing_angler > -trx(fly).wing_anglel;
  data{i}(idx) = trx(fly).wing_lengthr_mm(idx);
  
end
units = parseunits('mm');

