function dcentral = dcentral_pair(trx,larva1,larva2)

% initialize
dcentral = nan(1,trx(larva1).nframes);

% get start and end frames of overlap
t0 = max(trx(larva1).firstframe,trx(larva2).firstframe);
t1 = min(trx(larva1).endframe,trx(larva2).endframe);
  
% no overlap
if t1 < t0, 
  return;
end
  
% indices for these frames
i0 = t0 + trx(larva1).off;
i1 = t1 + trx(larva1).off;
j0 = t0 + trx(larva2).off;
j1 = t1 + trx(larva2).off;

% centroid distance
dx = trx(larva2).xcentral_mm(j0:j1)-trx(larva1).xcentral_mm(i0:i1);
dy = trx(larva2).ycentral_mm(j0:j1)-trx(larva1).ycentral_mm(i0:i1);
z = sqrt(dx.^2 + dy.^2);
dcentral(i0:i1) = z;
