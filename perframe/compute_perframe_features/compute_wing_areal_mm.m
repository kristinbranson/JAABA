function [data,units] = compute_wing_areal_mm(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = trx(fly).wing_areal ./ trx.pxpermm^2;
  
end
units = parseunits('mm^2');

