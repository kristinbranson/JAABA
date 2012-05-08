% change in max anglesub
function [data,units] = compute_danglesub(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = diff(trx(fly).anglesub) ./ trx(fly).dt;
end
units = parseunits('rad/s');