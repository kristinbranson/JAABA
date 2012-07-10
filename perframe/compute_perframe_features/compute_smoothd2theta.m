% smoothed acceleration in orientation
function [data,units] = compute_smoothd2theta(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);
  if trx(fly).nframes < 2,
    data{i} = [];
  elseif trx(fly).nframes == 2,
    data{i} = 0;  
  else
    data{i} = [0,diff(trx(fly).smoothdtheta,1,2)]./trx(fly).dt;
  end
end
units = parseunits('rad/s/s');

