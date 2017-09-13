function hfigs = MeasureStatDependencies(X,y,datafns,ti)

if nargin < 4,
  ti = '';
end

[ntrials,nstats] = size(X);

% Pearson's correlation coefficient
rho_spearman = nan(1,nstats);
pval_spearman = nan(1,nstats);
for i = 1:nstats,
  if ~any(~isnan(X(:,i))),
    continue;
  end
  [rho_spearman(i),pval_spearman(i)] = corr(X(~isnan(X(:,i)),i),y(~isnan(X(:,i)))','type','Spearman');
end
[~,order_spearman] = sort(pval_spearman);

% mean per session
idx = ~isnan(X);
n0 = hist(y,1:max(y));
m0 = nan(max(y),nstats);
s0 = nan(max(y),nstats);
for i = 1:nstats,
  m0(:,i) = accumarray(y(idx(:,i))',X(idx(:,i),i),[max(y),1],@mean,nan);
  s0(:,i) = accumarray(y(idx(:,i))',X(idx(:,i),i),[max(y),1],@std,nan);
end

statidxsig = find(pval_spearman <= .05);
nr = ceil((numel(statidxsig)+2)/2);

colors = jet(numel(statidxsig))*.7;

hfigs(1) = gcf;
clf;
hax = subplot(nr,2,1);
% pos = get(hax,'Position');
% pos(2) = pos(2) + .1;
% pos(4) = pos(4)-.1;
% set(hax,'Position',pos);
hold off;
plot([0,nstats+1],[.05,.05],'c--');
hold on;
hspearman = plot(1:nstats,pval_spearman(order_spearman),'k.-');
set(gca,'XTick',1:nstats,'XTickLabel',datafns(order_spearman),...
  'XLim',[0,nstats+1]);
htick = rotateticklabel(gca);
ylim = get(gca,'YLim');
ylim(1) = min(pval_spearman);
set(gca,'YLim',ylim);
ylabel('p-value, Spearman correlation coefficient');
box off;
title(ti);

[~,order] = sort(pval_spearman(statidxsig));

for i = 1:numel(statidxsig),
  set(htick(i),'Color',colors(i,:));
end

axi = reshape(1:nr*2,[2,nr]);
axi(1,1:2) = nan;
putxtick = false(2,nr);
putxtick(:,end) = true;
putxtick = putxtick(~isnan(axi));
axi = axi(~isnan(axi));


for i = 1:numel(statidxsig),
  subplot(nr,2,axi(i));
  hold off;
  errorbar(1:max(y),m0(:,statidxsig(order(i))),s0(:,statidxsig(order(i)))./sqrt(n0'),'o','Color',colors(i,:),'MarkerFaceColor',colors(i,:));
  axisalmosttight;
  set(gca,'XLim',[0,max(y)+1]);
  if ~putxtick(i),
    set(gca,'XTickLabel',{});
  end
  box off;
  ylabel(datafns{statidxsig(order(i))},'Interpreter','none');
end


