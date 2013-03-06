% change in quarter-major axis length
function [data,units] = compute_da(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  if numel(trx(fly).dt)>0
    data{i} = diff(trx(fly).a_mm)./trx(fly).dt;
  else
    data{i} = [];
  end
end
units = parseunits('mm/s');

