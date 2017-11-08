function hfigs = MeasureStatDependenciesPerGroup(X,y,groupidx,datafns,groupnames,ti)

if nargin < 4,
  ti = '';
end

[ntrials,nstats] = size(X);
ngroups = max(groupidx);

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
m0 = nan(max(y),nstats);
s0 = nan(max(y),nstats);
n0 = hist(y,1:max(y));
for i = 1:nstats,
  idxcurr = idx(:,i);
  m0(:,i) = accumarray(y(idxcurr)',X(idxcurr,i),[max(y),1],@mean,nan);
  s0(:,i) = accumarray(y(idxcurr)',X(idxcurr,i),[max(y),1],@std,nan);
end

% Pearson's correlation coefficient per group
rho_spearman_pergroup = nan(nstats,ngroups);
pval_spearman_pergroup = nan(nstats,ngroups);
for groupi = 1:ngroups,
  for i = 1:nstats,
    idxcurr = ~isnan(X(:,i)) & (groupidx(:)==groupi);
    if ~any(idxcurr),
      continue;
    end
    [rho_spearman_pergroup(i,groupi),pval_spearman_pergroup(i,groupi)] = corr(X(idxcurr,i),y(idxcurr)','type','Spearman');
  end
end

% mean per session, per group
idx = ~isnan(X);
n0_pergroup = nan(max(y),ngroups);
m0_pergroup = nan(max(y),nstats,ngroups);
s0_pergroup = nan(max(y),nstats,ngroups);
for groupi = 1:ngroups,
  groupidxcurr = groupidx==groupi;
  n0_pergroup(:,groupi) = hist(y(groupidxcurr),1:max(y));
  for i = 1:nstats,
    idxcurr = idx(:,i)&groupidxcurr(:);
    m0_pergroup(:,i,groupi) = accumarray(y(idxcurr)',X(idxcurr,i),[max(y),1],@mean,nan);
    s0_pergroup(:,i,groupi) = accumarray(y(idxcurr)',X(idxcurr,i),[max(y),1],@std,nan);
  end
end

doshowstd = max(s0_pergroup(:)) > 1;

statidxsig = find(pval_spearman <= .05);
nr = ceil((numel(statidxsig)+2)/2);

colors = jet(ngroups)*.7;
%markers0 = {'x','+','*','o','s','d','^','v','<','>',};
%markers = [repmat(markers0,[1,floor(ngroups/numel(markers0))]),markers0(1:mod(ngroups,numel(markers0)))];
markers = repmat({'s'},[1,ngroups]);
offx = linspace(-.125,.125,ngroups);

hfigs(1) = gcf;
clf;
hax = subplot(nr,2,1);
% pos = get(hax,'Position');
% pos(2) = pos(2) + .1;
% pos(4) = pos(4)-.1;
% set(hax,'Position',pos);
hold off;
for i = 0:nstats,
  plot(i+.5+[0,0],[0,1],':','Color',[.5,.5,.5]);
  hold on;
end
plot([0,nstats+1],[.05,.05],'c--');
hspearman = nan(1,ngroups+1);
for groupi = 1:ngroups,
  hspearman(groupi+1) = plot(1:nstats,pval_spearman_pergroup(order_spearman,groupi),markers{groupi},'Color',colors(groupi,:));
end
hspearman(1) = plot(1:nstats,pval_spearman(order_spearman),'k.-');
ylim = [0,1];
set(hax,'XTick',1:nstats,'XTickLabel',datafns(order_spearman),...
  'XLim',[0,nstats+1],'YLim',ylim);
htick = rotateticklabel(gca);
ylabel('p-value, Spearman correlation coefficient');
box off;
title(ti);

[~,order] = sort(pval_spearman(statidxsig));

for i = 1:numel(statidxsig),
  set(htick(i),'Color',[.7,0,0]);
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
  h = nan(1,ngroups);
  errorbar((1:max(y)),m0(:,statidxsig(order(i))),s0(:,statidxsig(order(i)))./sqrt(n0'),'ko-','MarkerFaceColor','k');
  hold on;
  for groupi = 1:ngroups,
    if doshowstd,
      h(groupi) = errorbar((1:max(y))+offx(groupi),m0_pergroup(:,statidxsig(order(i)),groupi),s0_pergroup(:,statidxsig(order(i)),groupi)./sqrt(n0_pergroup(:,groupi)),[markers{groupi}],'Color',colors(groupi,:));
    else
      h(groupi) = plot((1:max(y))+offx(groupi),m0_pergroup(:,statidxsig(order(i)),groupi),[markers{groupi}],'Color',colors(groupi,:));
    end
  end
  axisalmosttight;
  ylim = get(gca,'YLim');
  for j = 0:max(y),
    plot(j+.5+[0,0],ylim,':','Color',[.5,.5,.5]);
  end

  set(gca,'XLim',[0,max(y)+1]);
  if ~putxtick(i),
    set(gca,'XTickLabel',{});
  end
  box off;
  ylabel(datafns{statidxsig(order(i))},'Interpreter','none');
  if i == 1,
    legend(h,groupnames);
  end
  title(sprintf('rho = %f',rho_spearman(i)));
end


