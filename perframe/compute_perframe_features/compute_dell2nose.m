% minimum distance from point on ellipse of fly to nose of any other fly
function [data,units] = compute_dell2nose(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  % access closestfly to ensure that dell2nose is computed
  trx(fly1).closestfly_ell2nose;
  data{i1} = trx(fly1).dell2nose;
end
units = parseunits('mm');