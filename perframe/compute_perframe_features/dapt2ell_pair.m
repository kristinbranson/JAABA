function [d,angle] = dapt2ell_pair(trx,aptdata,aptLD, fly1,fly2,istry)

nsamples = 20;

% initialize
d = nan(1,trx(fly1).nframes);
angle = nan(1,trx(fly1).nframes);

% get start and end frames of overlap
t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
t1 = min(trx(fly1).endframe,trx(fly2).endframe);
  
% no overlap
if t1 < t0, 
  return;
end

% position of apt landmark aptLD
pTrk = aptdata.pTrk{fly1};
xapt = squeeze(pTrk(aptLD,1,:));
yapt = squeeze(pTrk(aptLD,2,:));

% % ellipse 2
% x_mm2 = trx(fly2).x_mm;
% y_mm2 = trx(fly2).y_mm;
% a_mm2 = trx(fly2).a_mm;
% b_mm2 = trx(fly2).b_mm;
% theta_mm2 = trx(fly2).theta_mm;

% ellipse 2
x_mm2 = trx(fly2).x;
y_mm2 = trx(fly2).y;
a_mm2 = trx(fly2).a;
b_mm2 = trx(fly2).b;
theta_mm2 = trx(fly2).theta;

off1 = trx(fly1).off;
off2 = trx(fly2).off;

if nargin < 6,
  tstry = t0:t1;
else
  tstry = istry(:)' - off1;
end

for t = tstry,
  i = t + off1;
  j = t + off2;
  [d(i),~,~,angle(i)] = ellipsedist_hack(x_mm2(j),y_mm2(j),...
    2*a_mm2(j),2*b_mm2(j),theta_mm2(j),...
    xapt(i),yapt(i),nsamples);
end
angle = angle - pi;