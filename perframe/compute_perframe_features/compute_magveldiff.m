% magnitude of difference in velocity between fly and closest fly based on
% type. 
function [data,units] = compute_magveldiff(trx,n,type)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);

for i1 = 1:nflies,

  fly1 = flies(i1);
  
  % fly closest to fly1 according to type
  closestfly = trx(fly1).(['closestfly_',type]);
  
  % velocity of fly1
  dx1 = diff(trx(fly1).x_mm,1,2);
  dy1 = diff(trx(fly1).y_mm,1,2);

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

    % velocity from these frames to the next frame
    dx2 = trx(fly2).x_mm(off+idx+1)-trx(fly2).x_mm(off+idx);
    dy2 = trx(fly2).y_mm(off+idx+1)-trx(fly2).y_mm(off+idx);

    % magnitude of difference between velocity vectors
    data{i1}(idx) = sqrt((dx1(idx)-dx2).^2 + (dy1(idx)-dy2).^2);
  end
end

units = parseunits('mm/s');