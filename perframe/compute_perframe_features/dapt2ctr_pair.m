function [d] = dapt2ctr_pair(trx,aptdata,aptLD,fly1,fly2,istry)

%initalize
d = nan(1,trx(fly1).nframes);

% get start and end frames of overlap
t0 = max(trx(fly1).firstframe,trx(fly2).firstframe);
t1 = min(trx(fly1).endframe,trx(fly2).endframe);
  
% no overlap
if t1 < t0, 
  return;
end

% position of apt landmark aptLD (fly1)
pTrk = aptdata.pTrk{fly1};
xLD = squeeze(pTrk(aptLD,1,:));
yLD = squeeze(pTrk(aptLD,2,:));

% fly 2 
a = (trx(fly2).a);
b = (trx(fly2).b);
theta = trx(fly2).theta;
x = trx(fly2).x;
y = trx(fly2).y;

off1 = trx(fly1).off;
off2 = trx(fly2).off;

% compute malahanobis distance only where the two trajecories overlap
if nargin < 6,
  tstry = t0:t1;
else
  tstry = istry(:)' - off1;
end

% should rewrite to use malah_dist vector computation! 
for t = tstry
    j = t + off2; % fly2
    i = t + off1; % fly1
    [d(i)] = malah_dist(x(j),y(j),a(j),b(j),theta(j),xLD(i),yLD(i));
end