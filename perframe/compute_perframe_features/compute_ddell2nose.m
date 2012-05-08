% change in min distance from any point on the ellipse of one fly to the
% nose of another fly
function [data,units] = compute_ddell2nose(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  data{i} = diff(trx(fly).dell2nose) ./ trx(fly).dt;
end
units = parseunits('mm/s');