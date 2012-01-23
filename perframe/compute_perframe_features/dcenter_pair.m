function dcenter = dcenter_pair(trx,fly1,fly2)

% initialize
dcenter = nan(1,trx(fly1).nframes);

% get start and end frames of overlap
t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
t1 = min(trx(fly1).endframe,trx(fly2).endframe);
  
% no overlap
if t1 < t0, 
  return;
end
  
% indices for these frames
i0 = t0 + trx(fly1).off;
i1 = t1 + trx(fly1).off;
j0 = t0 + trx(fly2).off;
j1 = t1 + trx(fly2).off;

% centroid distance
dx = trx(fly2).x_mm(j0:j1)-trx(fly1).x_mm(i0:i1);
dy = trx(fly2).y_mm(j0:j1)-trx(fly1).y_mm(i0:i1);
z = sqrt(dx.^2 + dy.^2);
dcenter(i0:i1) = z;
