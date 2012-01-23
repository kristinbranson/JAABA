% minimum distance from nose of fly to tail of any other fly
function [data,units] = compute_dnose2tail(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  % access closestfly to ensure that dnose2tail is computed
  trx(fly1).closestfly_nose2tail;
  data{i1} = trx(fly1).dnose2tail;
end
units = parseunits('mm');
