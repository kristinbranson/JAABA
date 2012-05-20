function hfig = PlotSampleFrames(trx,readframe,predictions,mainfly,otherflies,ts,varargin)

colorpos = [.7,0,0];
colorneg = [0,0,.7];
border = 20;

[colorpos,colorneg,border,hfig,figpos] = ...
  myparse(varargin,'colorpos',colorpos,...
  'colorneg',colorneg,...
  'border',border,...
  'hfig',[],...
  'figpos',[]);

if isempty(hfig),
  hfig = figure;
end

figure(hfig);
clf;
hax = createsubplots(1,numel(ts),.01);

xlim = [inf,-inf];
ylim = [inf,-inf];
for fly = [mainfly,otherflies],
  x = trx(fly).x(ts+trx(fly).off);
  y = trx(fly).y(ts+trx(fly).off);
  xlim(1) = min([xlim(1),x]);
  xlim(2) = max([xlim(2),x]);
  ylim(1) = min([ylim(1),y]);
  ylim(2) = max([ylim(2),y]);
end

ax = [xlim(1)-border,xlim(2)+border,ylim(1)-border,ylim(2)+border];

fly = mainfly;
x = trx(fly).x(ts(1)+trx(fly).off:ts(end)+trx(fly).off);
y = trx(fly).y(ts(1)+trx(fly).off:ts(end)+trx(fly).off);
idxpos = predictions{fly}(ts(1):ts(end));
for i = 1:numel(ts),
  t = ts(i);
  im = readframe(t);
  imagesc(im,'Parent',hax(i),[0,255]);
  axis(hax(i),'image','off');
  axis(hax(i),ax);
  hold(hax(i),'on');
  plot(hax(i),x,y,'k.-');
  plot(hax(i),x(idxpos),y(idxpos),'.','Color',colorpos);
  plot(hax(i),x(~idxpos),y(~idxpos),'.','Color',colorneg);
  if predictions{fly}(t),
    colorcurr = colorpos;
  else
    colorcurr = colorneg;
  end
  plot(hax(i),trx(fly).x(t+trx(fly).off),trx(fly).y(t+trx(fly).off),'o','color',colorcurr,'markerfacecolor',colorcurr);
end
colormap gray;

if isempty(figpos),
  truesize;
else
  set(hfig,'Position',figpos);
end