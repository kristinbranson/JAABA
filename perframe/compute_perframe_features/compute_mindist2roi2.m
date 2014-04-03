% distance to closest roi
function [data,units] = compute_mindist2roi2(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i = 1:nflies,
  fly = flies(i);
  % access closestroi2 to ensure that mindist2roi2 is computed
  trx(fly).closestroi2;
  data{i} = trx(fly).mindist2roi2;
end
units = parseunits('mm');
