function [data,units] = compute_dnwingsdetected(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = diff(trx(fly).nwingsdetected) ./ trx(fly).dt;
  
end
units = parseunits('unit/s');

