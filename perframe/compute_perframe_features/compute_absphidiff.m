% absolute difference in velocity direction between a fly and the closest
% fly according to type
function [data,units] = compute_absphidiff(trx,n,type)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  
  % fly closest to fly1 according to type
  closestfly = trx(fly1).(['closestfly_',type]);
  
  % velocity direction of fly1
  phi1 = trx(fly1).phi;

  % loop over all flies
  for i2 = 1:nflies,
    
    fly2 = flies(i2);
    if i1 == i2, continue; end
    
    % frames where this fly is closest
    idx = find(closestfly(1:end-1) == fly2);
    
    % don't use the last frame of fly2
    off = trx(fly1).firstframe - trx(fly2).firstframe;
    idx(idx+off == trx(fly2).nframes) = [];
    
    if isempty(idx), continue; end
    
    % orientation of fly2
    off = trx(fly1).firstframe - trx(fly2).firstframe;
    phi2 = trx(fly2).phi(off+idx);
    
    % absolute difference in orientation
    data{i1}(idx) = abs(modrange(phi2 - phi1(idx),-pi,pi));

  end
end

units = parseunits('rad');