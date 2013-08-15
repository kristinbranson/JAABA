% center distance to closest larva
function [data,units] = compute_dcenter(trx,n)

larvae = trx.exp2flies{n};
nlarvae = numel(larvae);
data = cell(1,nlarvae);

for i1 = 1:nlarvae,
  larva1 = larvae(i1);
  % access closestfly to ensure that dcenter is computed
  trx(larva1).closestlarva_center;
  data{i1} = trx(larva1).dcenter;
end
units = parseunits('mm');
