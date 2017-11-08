function hfig = CompareStatsAcrossMultipleGroups(X0,ys,ynames,idxs,idxnames,datafns,pairsignore)

[ntrials,nstats] = size(X0);
ny = numel(ynames);
nidx = numel(idxnames);

if ~exist('pairsignore','var'),
  pairsignore = cell(0,2);
end

pval_wil = nan(nstats,ny,nidx);
wilstat = nan(nstats,ny,nidx);
doignore = false(ny,nidx);

for yi = 1:ny,
  y0 = ys(yi,:);
  for idxi = 1:nidx,
    idx = idxs(idxi,:);
    
    doignore(yi,idxi) = any(strcmpi(ynames{yi},pairsignore(:,1))&strcmpi(idxnames{idxi},pairsignore(:,2)));
    if doignore(yi,idxi),
      continue;
    end
    
    X = X0(idx,:);
    y = y0(idx);

    % need at least one example per class
    n0 = sum(~isnan(X(y==0,:)),1);
    n1 = sum(~isnan(X(y==1,:)),1);
    goodidx = find(n0 > 1 & n1 > 1);

    % Wilcoxon rank-sum test
    for i = goodidx,
      [pval_wil(i,yi,idxi),~,pstats_wil] = ranksum(X(y'==0&~isnan(X(:,i)),i),X(y'==1&~isnan(X(:,i)),i));
      wilstat(i,yi,idxi) = pstats_wil.ranksum;
    end
    pval_wil(isnan(pval_wil)) = 1;

  end
end

pval_wil = reshape(pval_wil,[nstats,ny*nidx]);
pval_wil(:,doignore) = [];

pval_color = nan(nstats*nnz(~doignore),3);
binedges = [0,.0001,.001,.01,.05,1+eps];
nbins = numel(binedges)-1;
colors = [linspace(.8,0,nbins)',zeros(nbins,2)];
for i = 1:nbins,
  for c = 1:3,
    pval_color(pval_wil>=binedges(i)&pval_wil<binedges(i+1),c) = colors(i,c);
  end
end
pval_color = reshape(pval_color,[nstats,nnz(~doignore),3]);

hfig = gcf;
clf;
image(pval_color);
set(gca,'CLim',[0,3]);

xticklabels = cell(ny,nidx);
for yi = 1:ny,
  for idxi = 1:nidx,
    xticklabels{yi,idxi} = [ynames{yi},', ',idxnames{idxi}];
  end
end
xticklabels(doignore) = [];

set(gca,'XTick',1:ny*nidx,'XTickLabel',xticklabels,'YTick',1:nstats,'YTickLabel',datafns);
box off;