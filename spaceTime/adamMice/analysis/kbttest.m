function pval_ttest = kbttest(X,y)

y = y(:)';

n0 = sum(~isnan(X(y==0,:)),1);
n1 = sum(~isnan(X(y==1,:)),1);
goodidx = find(n0 > 1 & n1 > 1); % AL: n0>=1 & n1>=1?

m0 = nanmean(X(y==0,:),1);
v0 = nanvar(X(y==0,:),1,1); % AL: nanvar(~,0,1)?
m1 = nanmean(X(y==1,:),1);
v1 = nanvar(X(y==1,:),1,1); % AL: nanvar(~,0,1)?
df = (v0./n0 + v1./n1).^2 ./ ...
  max(eps,( (v0./n0).^2./(n0-1) + (v1./n1).^2./(n1-1) ));

tstat = abs((m1-m0)./ max(eps,sqrt(v1./n1 + v0./n0)));
pval_ttest = 2 * tcdf(-abs(tstat),df);
pval_ttest(tstat==0) = 1;