function [data,units] = compute_max_dwing_angle_out(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);

  data{i} = max(-diff(trx(fly).wing_anglel),diff(trx(fly).wing_angler)) ./ trx(fly).dt;
  
end
units = parseunits('rad/s');

