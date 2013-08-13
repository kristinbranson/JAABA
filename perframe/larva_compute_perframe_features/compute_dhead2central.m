% minimum distance from head of fly to central of any other fly
function [data,units] = compute_dhead2central(trx,n)

larvae = trx.exp2flies{n};
nlarvae = numel(larvae);
data = cell(1,nlarvae);

for i1 = 1:nlarvae,
  larva1 = larvae(i1);
  % access closestfly to ensure that dhead2central is computed
  trx(larva1).closestlarva_head2central;
  data{i1} = trx(larva1).dhead2central;
end
units = parseunits('mm');