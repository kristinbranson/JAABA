% change in quarter-minor axis length
function [data,units] = compute_db(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = diff(trx(fly).b_mm)./trx(fly).dt;
end
units = parseunits('mm/s');

