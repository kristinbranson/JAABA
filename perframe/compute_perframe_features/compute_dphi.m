% change in velocity direction
function [data,units] = compute_dphi(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  if numel(trx(fly).dt)>0
  data{i} = modrange(diff(trx(fly).phi),-pi,pi)./trx(fly).dt;
  else
    data{i} = [];
  end
end
units = parseunits('rad/s');

