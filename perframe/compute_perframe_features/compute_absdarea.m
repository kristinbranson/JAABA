% absolute change in area
function [data,units] = compute_absdarea(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = abs(trx(fly).darea);
end
units = parseunits('mm^2/s');