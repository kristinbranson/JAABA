function [hweight,hscore,hax,hfig,hylabel,hticks,hcolorbar,...
  sorted_weights,feature_order,bins,scores] = ...
  ShowWindowFeatureWeights(data,varargin)

[hax,hfig,figpos,...
  nfeatures_show] = ...
  myparse(varargin,'axes',[],'figure',[],...
  'figpos',[],...
  'nfeatures_show',numel(data.classifier));

didcreateaxes = false;
if numel(hax) < 2,
  if isempty(hfig),
    hfig = figure;
  else
    if ~ishandle(hfig),
      figure(hfig);
    else
      clf(hfig);
    end
  end
  hax = [0,0];
  axes_height = .375;
  axes_border = .025;
  bottom_border = .2;
  hax(1) = axes('Position',[.05,bottom_border+axes_height+axes_border,.925,axes_height]);
  hax(2) = axes('Position',[.05,bottom_border,.925,axes_height]);
  didcreateaxes = true;
else
  hfig = get(hax(1),'Parent');
end
hylabel = zeros(1,2);
if ~isempty(figpos),
  set(hfig,'Position',figpos);
end

holdstate = false(1,2);
for i = 1:2,
  holdstate(i) = ishold(hax(i));
  hold(hax(i),'off');
end

[sorted_weights,feature_order,bins,scores] = ...
  sortWindowFeaturesByWeight(data.classifier,data.windowdata.X);

nfeatures_show = min(nfeatures_show,find(sorted_weights <= 0,1)-1);

hweight = plot(hax(1),1:nfeatures_show,sorted_weights(1:nfeatures_show),'k-','linewidth',3);
hylabel(1) = ylabel(hax(1),'Total absolute weight');
set(hax(1),'XLim',[.5,nfeatures_show+.5],'XTick',1:nfeatures_show,'XTickLabel',{});

maxneg = max(max(-scores(:,1:nfeatures_show)));
maxpos = max(max(scores(:,1:nfeatures_show)));
ncolors = 128;
npos = round(maxpos / (maxpos + maxneg) * (ncolors-1));
nneg = ncolors - 1 - npos;
cm = [[zeros(nneg,2),linspace(.7,0,nneg)']
  [0,0,0]
  [linspace(0,.7,npos)',zeros(npos,2)]];

hscore = imagesc([1,nfeatures_show],[0,100],scores(:,1:nfeatures_show),'parent',hax(2));
axis(hax(2),'xy');
colormap(hax(2),cm);

print_names = cell(1,nfeatures_show);
for ii = 1:nfeatures_show,
  i = feature_order(ii);
  print_names{ii} = data.windowdata.featurenames{i}{1};
  for j = 3:2:numel(data.windowdata.featurenames{i}),
    if ischar(data.windowdata.featurenames{i}{j}),
      print_names{ii} = [print_names{ii},'_',data.windowdata.featurenames{i}{j}];
    else
      print_names{ii} = [print_names{ii},'_',data.windowdata.featurenames{i}{j-1},num2str(data.windowdata.featurenames{i}{j})];
    end
  end
end

set(hax(2),'XTick',1:nfeatures_show,'XTickLabel',print_names);
hticks = rotateticklabel(hax(2),90);
hylabel(2) = ylabel('Prctile of training data');

hcolorbar = colorbar('peer',hax(2),'Location','East');
set(hcolorbar,'XColor','w','YColor','w')

linkaxes(hax,'x');
set(hax,'XLim',[.5,nfeatures_show+.5],'Units','normalized');

% make labels fit
if didcreateaxes,
  axes_pos = get(hax,'Position');
  set(hticks,'Units','normalized');
  miny1 = min(cellfun(@(x) x(2), get(hticks,'Extent')));
  % put in figure's units
  miny = axes_pos{2}(2) + miny1*axes_pos{2}(4);
  
  if miny < 0,
    new_height = axes_height + miny/2;
    new_bottom_border = bottom_border - miny;
    axes_pos{1} = [axes_pos{1}(1),new_bottom_border+new_height+axes_border,axes_pos{1}(3),new_height];
    axes_pos{2} = [axes_pos{2}(1),new_bottom_border,axes_pos{2}(3),new_height];
    for i = 1:2,
      set(hax(i),'Position',axes_pos{i});
    end
  end
end
for i = 1:2,
  if holdstate(i),
    hold(hax(i),'on');
  end
end