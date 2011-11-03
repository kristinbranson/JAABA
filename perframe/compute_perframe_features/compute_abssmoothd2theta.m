% smoothed acceleration in orientation
function [data,units] = compute_abssmoothd2theta(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = abs(trx(fly).smoothd2theta);
end
units = parseunits('rad/s/s');

