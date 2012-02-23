function d = dell2nose_pair(trx,fly1,fly2,istry)

nsamples = 20;

% initialize
d = nan(1,trx(fly1).nframes);

% get start and end frames of overlap
t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
t1 = min(trx(fly1).endframe,trx(fly2).endframe);
  
% no overlap
if t1 < t0, 
  return;
end

% position of nose2
xnose = trx(fly2).x_mm + 2*trx(fly2).a_mm.*cos(trx(fly2).theta_mm);
ynose = trx(fly2).y_mm + 2*trx(fly2).a_mm.*sin(trx(fly2).theta_mm);

% ellipse 1
x_mm1 = trx(fly1).x_mm;
y_mm1 = trx(fly1).y_mm;
a_mm1 = trx(fly1).a_mm;
b_mm1 = trx(fly1).b_mm;
theta_mm1 = trx(fly1).theta_mm;

off1 = trx(fly1).off;
off2 = trx(fly2).off;

if nargin < 4,
  tstry = t0:t1;
else
  tstry = istry(:)' - off1;
end

for t = tstry,
  i = t + off1;
  j = t + off2;
  d(i) = ellipsedist_hack(x_mm1(i),y_mm1(i),...
    2*a_mm1(i),2*b_mm1(i),theta_mm1(i),...
    xnose(j),ynose(j),nsamples);
end