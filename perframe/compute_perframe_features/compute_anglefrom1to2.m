% azimuthal angle to from fly to closest fly according to type
function [data,units] = compute_anglefrom1to2(trx,n,type)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  
  % fly closest to fly1 according to type
  closestfly = trx(fly1).(['closestfly_',type]);
  
  % position of fly1
  xnose_mm1 = trx(fly1).xnose_mm;
  ynose_mm1 = trx(fly1).ynose_mm;
  theta_mm1 = trx(fly1).theta_mm;

  % loop over all flies
  for i2 = 1:nflies,
    
    fly2 = flies(i2);
    if i1 == i2, continue; end
    
    % frames where this fly is closest
    idx = find(closestfly == fly2);
    if isempty(idx), continue; end
    
    % angle to fly2
    off = trx(fly1).firstframe - trx(fly2).firstframe;
    dx2 = trx(fly2).x_mm(off+idx)-xnose_mm1(idx);
    dy2 = trx(fly2).y_mm(off+idx)-ynose_mm1(idx);
    theta2 = atan2(dy2,dx2);
    
    % angle relative to fly1's orientation
    data{i1}(idx) = modrange(theta2 - theta_mm1(idx),-pi,pi);

  end
end

units = parseunits('rad');