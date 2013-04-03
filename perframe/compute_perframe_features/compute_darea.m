% change in area
function [data,units] = compute_darea(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  if numel(trx(fly).dt)>0
    data{i} = diff(trx(fly).area)./trx(fly).dt;
  else
    data{i} = [];
  end
end
units = parseunits('mm*mm/s');

