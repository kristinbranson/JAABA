function d = dnose2ell_anglerange_pair(trx,fly1,fly2,anglerange)

nsamples = 20;

anglerange(2) = modrange(anglerange(2),anglerange(1),anglerange(1)+2*pi);

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
xnose1 = trx(fly1).x_mm(i0:i1) + 2*trx(fly1).a_mm(i0:i1).*cos(trx(fly1).theta_mm(i0:i1));
ynose1 = trx(fly1).y_mm(i0:i1) + 2*trx(fly1).a_mm(i0:i1).*sin(trx(fly1).theta_mm(i0:i1));
theta1 = trx(fly1).theta_mm(i0:i1);

% position of center2
x_mm2 = trx(fly2).x_mm(j0:j1);
y_mm2 = trx(fly2).y_mm(j0:j1);

% anglefrom1to2
dx2 = x_mm2-xnose1;
dy2 = y_mm2-ynose1;
theta2 = atan2(dy2,dx2);
anglefrom1to2 = modrange(theta2-theta1,anglerange(1),anglerange(1)+2*pi);

% indices within range
idx = anglefrom1to2 >= anglerange(1) & anglefrom1to2 <= anglerange(2);

if ~any(idx),
  return;
end

% ellipse of fly 2
a_mm2 = trx(fly2).a_mm(j0:j1);
b_mm2 = trx(fly2).b_mm(j0:j1);
theta_mm2 = trx(fly2).theta_mm(j0:j1);

for i = find(idx(:)'),
  j = i + i0 - 1;
  d(j) = ellipsedist_hack(x_mm2(i),y_mm2(i),...
    2*a_mm2(i),2*b_mm2(i),theta_mm2(i),...
    xnose1(i),ynose1(i),nsamples);
end