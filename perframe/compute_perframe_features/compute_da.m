% change in quarter-major axis length
function [data,units] = compute_da(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = diff(trx(fly).a_mm)./trx(fly).dt;
end
units = parseunits('mm/s');

