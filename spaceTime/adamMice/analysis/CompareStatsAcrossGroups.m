function hfigs = CompareStatsAcrossGroups(X,y,datafns,groupnames,varargin)

[ti,ax] = myparse(varargin,...
  'title','',...
  'axes',[]);

[ntrials,nstats] = size(X);

assert(isvector(y) && numel(y)==ntrials);
assert(all(y==0 | y==1),'y must be a binary grouping vector');
assert(iscellstr(datafns) && numel(datafns)==nstats);
assert(iscellstr(groupnames)); % currently groupnames must be two-el vec corresponding to y==[0;1] in that order

y = y(:)'; % needs to be a row vec

% need at least one example per class
n0 = sum(~isnan(X(y==0,:)),1);
n1 = sum(~isnan(X(y==1,:)),1);
goodidx = find(n0 > 1 & n1 > 1);

if isempty(goodidx)
  warningNoTrace('CompareStats:noData','Two data points are required for each group/class.');
  hfigs = [];
  return;
end

% t-test
m0 = nanmean(X(y==0,:),1);
v0 = nanvar(X(y==0,:),1,1);
m1 = nanmean(X(y==1,:),1);
v1 = nanvar(X(y==1,:),1,1);
df = (v0./n0 + v1./n1).^2 ./ ...
  max(eps,( (v0./n0).^2./(n0-1) + (v1./n1).^2./(n1-1) ));
m = nanmean(X,1);
s = max(eps,nanstd(X,1,1));
z0 = (m0 - m) ./ s;
z1 = (m1 - m) ./ s;
s0 = sqrt(v0)./s;
s1 = sqrt(v1)./s;
Z = bsxfun(@rdivide,bsxfun(@minus,X,m),s);

tstat = abs((m1-m0)./ max(eps,sqrt(v1./n1 + v0./n0)));
pval_ttest = 2 * tcdf(-abs(tstat),df);
pval_ttest(tstat==0) = 1;
[~,order_ttest] = sort(pval_ttest(goodidx));
assert(~any(isnan(m0(goodidx))));
assert(~any(isnan(m1(goodidx))));

% Wilcoxon rank-sum test
wilstat = nan(1,nstats);
pval_wil = nan(1,nstats);
pval_kw = nan(1,nstats);
for i = goodidx,
  tfGrp0 = y'==0&~isnan(X(:,i));
  tfGrp1 = y'==1&~isnan(X(:,i));
  [pval_wil(i),~,pstats_wil] = ranksum(X(tfGrp0,i),X(tfGrp1,i));
  wilstat(i) = pstats_wil.ranksum;
  pval_kw(i) = kruskalwallis(X(tfGrp0|tfGrp1,i),y(tfGrp0|tfGrp1)','off');
end
pval_wil(isnan(pval_wil)) = 1;
pval_kw(isnan(pval_kw)) = 1;
% [~,order_wil] = sort(pval_wil(goodidx));
% [~,order_kw] = sort(pval_kw(goodidx));

if isempty(ax)
  hfigs(1) = gcf;
  subplot(2,1,1);
else
  hfigs = [];
  axes(ax(1));
end
hold off;
plot([0,numel(goodidx)+1],[.05,.05],'c--');
hold on;
httest = plot(1:numel(goodidx),pval_ttest(goodidx(order_ttest)),'k.-');
hwil = plot(1:numel(goodidx),pval_wil(goodidx(order_ttest)),'x-','Color',[.5,0,0]);
hkw = plot(1:numel(goodidx),pval_kw(goodidx(order_ttest)),'o-','Color',[0,0,1]);
set(gca,'XTick',1:numel(goodidx),'XTickLabel',{},'XLim',[0,numel(goodidx)+1]);
ylim = get(gca,'YLim');
ylim(1) = min(pval_ttest(goodidx));
set(gca,'YLim',ylim);
ylabel('p-value');
box off;
legend([httest,hwil,hkw],{'t-test','Wilcoxon rank-sum','KruskalWallis'},'location','northwest');
title(ti,'interpreter','none','fontsize',12);

if isempty(ax)
  hax = subplot(2,1,2);
  pos = get(hax,'Position');
  pos(2) = pos(2) + .1;
  set(hax,'Position',pos);
else
  hax = ax(2);
  axes(hax);
end
hold off;
% plot(repmat(1:numel(goodidx),[nnz(y==0),1])+rand([nnz(y==0),numel(goodidx)])*.1-.175,Z(y==0,goodidx(order_ttest)),'.','Color',[.5,0,.5]+.5);
% hold on;
% plot(repmat(1:numel(goodidx),[nnz(y==1),1])+rand([nnz(y==1),numel(goodidx)])*.1+.075,Z(y==1,goodidx(order_ttest)),'.','Color',[0,.5,.5]+.5);
errorbar((1:numel(goodidx))-.125,z0(goodidx(order_ttest)),s0(goodidx(order_ttest))./sqrt(n0(goodidx(order_ttest))),'o','Color',[.5,0,.5],'MarkerFaceColor',[.5,0,.5]);
hold on;
errorbar((1:numel(goodidx))+.125,z1(goodidx(order_ttest)),s1(goodidx(order_ttest))./sqrt(n1(goodidx(order_ttest))),'o','Color',[0,.5,.5],'MarkerFaceColor',[0,.5,.5]);
set(gca,'XTick',1:numel(goodidx),'XTickLabel',datafns(goodidx(order_ttest)),...
  'XLim',[0,numel(goodidx)+1]);
rotateticklabel(gca);
ylabel('z-score');
box off;
legend(groupnames,'interpreter','none');

