% absolute difference in orientation between a fly and the closest fly
% according to type
function [data,units] = compute_absthetadiff(trx,n,type)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  
  % fly closest to fly1 according to type
  closestfly = trx(fly1).(['closestfly_',type]);
  
  % orientation of fly1
  theta_mm1 = trx(fly1).theta_mm;

  % loop over all flies
  for i2 = 1:nflies,
    
    fly2 = flies(i2);
    if i1 == i2, continue; end
    
    % frames where this fly is closest
    idx = find(closestfly == fly2);
    if isempty(idx), continue; end
    
    % orientation of fly2
    off = trx(fly1).firstframe - trx(fly2).firstframe;
    theta_mm2 = trx(fly2).theta_mm(off+idx);
    
    % absolute difference in orientation
    data{i1}(idx) = abs(modrange(theta_mm2 - theta_mm1(idx),-pi,pi));

  end
end

units = parseunits('rad');