function d = dnose2tail_pair(trx,fly1,fly2)

% initialize
d = nan(1,trx(fly1).nframes);

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

% position of nose1
xnose = trx(fly1).x_mm(i0:i1) + 2*trx(fly1).a_mm(i0:i1).*cos(trx(fly1).theta_mm(i0:i1));
ynose = trx(fly1).y_mm(i0:i1) + 2*trx(fly1).a_mm(i0:i1).*sin(trx(fly1).theta_mm(i0:i1));

% position of tail2
xtail = trx(fly2).x_mm(j0:j1) - 2*trx(fly2).a_mm(j0:j1).*cos(trx(fly2).theta_mm(j0:j1));
ytail = trx(fly2).y_mm(j0:j1) - 2*trx(fly2).a_mm(j0:j1).*sin(trx(fly2).theta_mm(j0:j1));

dx = xtail-xnose;
dy = ytail-ynose;
z = sqrt(dx.^2 + dy.^2);
d(i0:i1) = z;