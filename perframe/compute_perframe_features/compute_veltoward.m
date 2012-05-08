% signed velocity in the direction of the closest fly, according to type
function [data,units] = compute_veltoward(trx,n,type)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,
  fly1 = flies(i1);
  
  % fly closest to fly1 according to type
  closestfly = trx(fly1).(['closestfly_',type]);
  
  % velocity of fly1
  dx1 = diff(trx(fly1).x_mm);
  dy1 = diff(trx(fly1).y_mm);
  x_mm1 = trx(fly1).x_mm;
  y_mm1 = trx(fly1).y_mm;

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
    
    % unit vector in direction of fly2 from fly1
    off = trx(fly1).firstframe - trx(fly2).firstframe;
    dx2 = trx(fly2).x_mm(off+idx)-x_mm1(idx);
    dy2 = trx(fly2).y_mm(off+idx)-y_mm1(idx);
    dz2 = sqrt(dx2.^2 + dy2.^2);
    dx2 = dx2 ./ dz2;
    dy2 = dy2 ./ dz2;
    dx2(dz2==0) = 0;
    dy2(dz2==0) = 0;
    
    % project velocity of fly1 onto this vector
    data{i1}(idx) = dx1(idx).*dx2 + dy1(idx).*dy2;

  end
end

units = parseunits('mm/s');