% area of ellipse
function [data,units] = compute_area(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = (2*trx(fly).a_mm).*(2*trx(fly).b_mm)*pi;
end
units = parseunits('mm^2');

