% minimum distance from nose of fly to point on ellipse of any other fly
function [data,units] = compute_dnose2ell(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  % access closestfly to ensure that dcenter is computed
  trx(fly1).closestfly_nose2ell;
  data{i1} = trx(fly1).dnose2ell;
end
units = parseunits('mm');
