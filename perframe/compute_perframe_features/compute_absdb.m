% absolute change in quarter-minor axis length
function [data,units] = compute_absdb(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = abs(trx(fly).db);
end
units = parseunits('mm/s');