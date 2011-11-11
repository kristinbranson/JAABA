% absolute difference between orientation and velocity direction, where
% sideways has a large value and forward/backward is 0
function [data,units] = compute_phisideways(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  % how sideways is the velocity direction?
  if trx(fly).nframes < 2,
    data{i} = [];
  else
    data{fly} = abs(modrange(trx(fly).phi-trx(fly).theta_mm,-pi/2,pi/2));
  end
end
units = parseunits('rad');

