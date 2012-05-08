% sign of change in body orientation
function [data,units] = compute_signdtheta(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = sign(trx(fly).dtheta);
end
units = parseunits('unit');

