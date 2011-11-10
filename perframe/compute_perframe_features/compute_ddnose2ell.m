% change in min distance from nose of one fly to any point on the ellipse
% of any other fly
function [data,units] = compute_ddnose2ell(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = diff(trx(fly).dnose2ell) ./ trx(fly).dt;
end
units = parseunits('mm/s');