% change in distance to closest fly center
function [data,units] = compute_ddcenter(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = diff(trx(fly).dcenter) ./ trx(fly).dt;
end
units = parseunits('mm/s');