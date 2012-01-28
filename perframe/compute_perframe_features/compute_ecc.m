% eccentricity
function [data,units] = compute_ecc(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = trx(fly).b_mm ./ trx(fly).a_mm;
end
units = parseunits('unit');
