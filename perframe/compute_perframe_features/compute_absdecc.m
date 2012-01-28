% absolute change in ecc
function [data,units] = compute_absdecc(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = abs(trx(fly).decc);
end
units = parseunits('unit/s');