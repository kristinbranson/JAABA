% absolute rotation of nose around mean tail location
function [data,units] = compute_absdtheta_tail(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  data{i} = abs(trx(fly).dtheta_tail);
end
units = parseunits('rad/s');

