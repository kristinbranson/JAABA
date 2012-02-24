function [d,i0,i1,j0,j1,t0,t1,anglefrom1to2] = dnose2center_pair(trx,fly1,fly2)

% initialize
d = nan(1,trx(fly1).nframes);
anglefrom1to2 = [];

% get start and end frames of overlap
t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
t1 = min(trx(fly1).endframe,trx(fly2).endframe);
  
% indices for these frames
i0 = t0 + trx(fly1).off;
i1 = t1 + trx(fly1).off;
j0 = t0 + trx(fly2).off;
j1 = t1 + trx(fly2).off;

% no overlap
if t1 < t0, 
  return;
end

% position of nose1
xnose1 = trx(fly1).x_mm(i0:i1) + 2*trx(fly1).a_mm(i0:i1).*cos(trx(fly1).theta_mm(i0:i1));
ynose1 = trx(fly1).y_mm(i0:i1) + 2*trx(fly1).a_mm(i0:i1).*sin(trx(fly1).theta_mm(i0:i1));

% position of center2
x_mm2 = trx(fly2).x_mm(j0:j1);
y_mm2 = trx(fly2).y_mm(j0:j1);

dx = x_mm2-xnose1;
dy = y_mm2-ynose1;

z = sqrt(dx.^2 + dy.^2);
d(i0:i1) = z;

% anglefrom1to2
if nargout >= 8,
  theta1 = trx(fly1).theta_mm(i0:i1);
  theta2 = atan2(dy,dx);
  anglefrom1to2 = modrange(theta2-theta1,-pi,pi);
end
