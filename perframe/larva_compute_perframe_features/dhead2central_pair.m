function d = dhead2central_pair(trx,larva1,larva2)

% initialize
d = nan(1,trx(larva1).nframes);

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

% position of head1
xhead = trx(larva1).xspine_mm(1,i0:i1); 
yhead = trx(larva1).yspine_mm(1,i0:i1);

% position of central2
xcentral = trx(larva2).xspine_mm(6,j0:j1); 
ycentral = trx(larva2).yspine_mm(6,j0:j1); 

dx = xcentral-xhead;
dy = ycentral-yhead;
z = sqrt(dx.^2 + dy.^2);
d(i0:i1) = z;