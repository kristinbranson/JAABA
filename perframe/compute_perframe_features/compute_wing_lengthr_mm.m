function [data,units] = compute_wing_lengthr_mm(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = trx(fly).wing_lengthr ./ trx.pxpermm;
  
end
units = parseunits('mm');

