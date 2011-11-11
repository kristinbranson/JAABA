% absolute value of the projection of the center of rotation on the minor
% axis
function [data,units] = compute_abscorfrac_min(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  
  data{i} = abs(trx(fly).corfrac_min);

end
units = parseunits('unit');

