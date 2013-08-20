% center distance to closest larva
function [data,units] = compute_dcentral(trx,n)

larvae = trx.exp2flies{n};
nlarvae = numel(larvae);
data = cell(1,nlarvae);

for i1 = 1:nlarvae,
  larva1 = larvae(i1);
  % access closestlarva to ensure that dcentral is computed
  trx(larva1).closestlarva_central;
  data{i1} = trx(larva1).dcentral;
end
units = parseunits('mm');
