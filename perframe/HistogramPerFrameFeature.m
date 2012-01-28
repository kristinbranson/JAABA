function [hhist,hax,hfig,hleg,hxlabel,hylabel,frac,frac_outside,edges,centers_plot] = ...
  HistogramPerFrameFeature(data,perframefn,varargin)

[hax,hfig,...
  minprctile,maxprctile,minv,maxv,...
  edges,nbins,...
  binmode,...
  labelcolors,unknowncolor,...
  frac,...
  frac_outside,...
  centers_plot] = ...
  myparse(varargin,'axes',[],'figure',[],...
  'minprctile',[],'maxprctile',[],...
  'minv',[],'maxv',[],...
  'edges',[],...
  'nbins',100,...
  'binmode','linear',...
  'labelcolors',[],...
  'unknowncolor',[1,1,1],...
  'frac',[],...
  'frac_outside',[],...
  'centers_plot',[]);

nbehaviors = data.nbehaviors;

if isempty(hax),
  if isempty(hfig),
    hfig = figure;
  end
  hax = gca;
end

holdstate = ishold(hax);
hold(hax,'off');

docompute = isempty(frac) || isempty(frac_outside) || ...
  isempty(edges) || isempty(centers_plot);

if docompute,
  
  perframeidx = find(strcmpi(perframefn,data.allperframefns),1);
  if isempty(perframeidx),
    error('Unknown perframe field %s',perframefn);
  end
  
  % use current fly to choose hist edges
  if isempty(edges),
    x = data.perframedata{perframeidx};
    if isempty(minv),
      if isempty(minprctile),
        minv = min(x);
      else
        minv = prctile(x,minprctile);
      end
    end
    if isempty(maxv),
      if isempty(minprctile),
        maxv = max(x);
      else
        maxv = prctile(x,maxprctile);
      end
    end
    [edges,centers] = SelectHistEdges(nbins,[minv,maxv],binmode);
  else
    centers = (edges(1:end-1)+edges(2:end))/2;
    nbins = numel(centers);
  end
  
  counts = zeros(nbehaviors+1,nbins);
  counts_outside = zeros(nbehaviors+1,2);
  n = zeros(nbehaviors+1,1);
  
  for expi = 1:data.nexps,
    
    % load per-frame data for this experiment
    perframedir = data.GetFile('perframedir',expi);
    file = fullfile(perframedir,[perframefn,'.mat']);
    if ~exist(file,'file'),
      warning('Per-frame data file %s does not exist',file);
      continue
    end
    perframedata = load(file);
    
    % TODO: extend to multiple flies
    for fly = 1:data.nflies_per_exp(expi),
      
      x = perframedata.data{fly};
      
      % load labels for this fly and experiment
      if expi == data.expi && all(fly == data.flies),
        % use current labelidx
        labelidx = data.labelidx.vals;
      else
        labelidx = data.GetLabelIdx(expi,fly);
        labelidx = labelidx.vals;
      end
      noff = numel(labelidx) - numel(x);
      noff_start = floor(noff/2);
      noff_end = noff - noff_start;
      labelidx = labelidx(1+noff_start:end-noff_end);
      
      % histogram
      for i = 0:nbehaviors,
        idx = ~isnan(x) & labelidx == i;
        n(i+1) = n(i+1) + nnz(idx);
        counts_curr = histc(x(idx),edges);
        counts(i+1,:) = counts(i+1,:) + counts_curr(1:end-1);
        counts_outside(i+1,1) = counts_outside(i+1,1) + nnz(x(idx)<edges(1));
        counts_outside(i+1,2) = counts_outside(i+1,1) + nnz(x(idx)>edges(end));
      end
    end
  end
  
  frac = bsxfun(@rdivide,counts,n);
  frac_outside = bsxfun(@rdivide,counts_outside,n);
  
  centers_plot = [2*centers(1)-centers(2),centers,2*centers(end)-centers(end-1)];
  
end

hhist = zeros(1,nbehaviors+1);

if size(labelcolors,1) < nbehaviors,
  labelcolors(size(labelcolors,1)+1:nbehaviors,:) = jet(nbehaviors-size(labelcolors,1));
end

for i = 1:nbehaviors+1,
  hhist(i) = plot(hax,centers_plot,...
    [frac_outside(i,1),frac(i,:),frac_outside(i,2)],...
    '-','linewidth',2);
  if i == 1,
    set(hhist(i),'Color',unknowncolor);
    hold(hax,'on');
  else
    set(hhist(i),'Color',labelcolors(i-1,:));
  end
end

axisalmosttight(1/20,hax);
xtick = get(hax,'XTick');
xticklabel = get(hax,'XTickLabel');
xticklabel = cellstr(xticklabel);
i = find(xtick >= centers_plot(5),1);
if ~isempty(i),
  xtick = xtick(i:end);
  xticklabel = xticklabel(i:end);
end
i = find(xtick <= centers_plot(end-4),1,'last');
if ~isempty(i),
  xtick = xtick(1:i);
  xticklabel = xticklabel(1:i);
end
xtick = [centers_plot(1),xtick,centers_plot(end)];
xticklabel = [{sprintf('<%.1f',edges(1))};xticklabel;{sprintf('>=%.1f',edges(end))}];
set(hax,'XTick',xtick,'XTickLabel',xticklabel);

s = [{'Unknown'},data.labelnames];
hleg = legend(hax,hhist(n>0),s(n>0));
hxlabel = xlabel(hax,perframefn,'Interpreter','none');
hylabel = ylabel(hax,'Fraction of frames');

if ~holdstate,
  hold(hax,'off');
end