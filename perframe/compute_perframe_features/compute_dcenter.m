% center distance to closest fly
function [data,units] = compute_dcenter(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  % access closestfly to ensure that dcenter is computed
  trx(fly1).closestfly_center;
  data{i1} = trx(fly1).dcenter;
end
units = parseunits('mm');
