function [hfig] = OverlayFrames_Ellipses(trx,readframe,predictions,mainfly,otherflies,ts,varargin)

colorpos = [.7,0,0];
colorneg = [0,0,.7];
border = 10;

[colorpos,colorneg,border,hfig,figpos,cm] = ...
  myparse(varargin,'colorpos',colorpos,...
  'colorneg',colorneg,...
  'border',border,...
  'hfig',[],...
  'figpos',[],...
  'colormap',{@jet,@jet});

if isempty(otherflies),
  flygroups = {mainfly};
else
  flygroups = {mainfly,otherflies};
end
nflygroups = numel(flygroups);
nframes = numel(ts);
allts = ts(1):ts(end);
nframes_total = numel(allts);
allflies = [flygroups{:}];

im0 = readframe(ts(1));
[nr0,nc0,~] = size(im0);

%% rotate everything so that things are as horizontal as possible

x = nan(trx.nflies,nframes_total);
y = nan(trx.nflies,nframes_total);
a = nan(trx.nflies,nframes_total);
b = nan(trx.nflies,nframes_total);
theta = nan(trx.nflies,nframes_total);
imcenter = [(nc0+1)/2,(nr0+1)/2];
for fly = 1:trx.nflies,
  idx = allts+trx(fly).off;
  i0 = find(idx>1,1);
  i1 = find(idx<trx(fly).nframes,1,'last');
  x(fly,i0:i1) = trx(fly).x(idx(i0:i1));
  y(fly,i0:i1) = trx(fly).y(idx(i0:i1));
  a(fly,i0:i1) = trx(fly).a(idx(i0:i1));
  b(fly,i0:i1) = trx(fly).b(idx(i0:i1));
  theta(fly,i0:i1) = trx(fly).theta(idx(i0:i1));
end

% % rotate based on covariance matrix of positions
% xflies = x(allflies,:);
% yflies = y(allflies,:);
% S = cov([xflies(:),yflies(:)]);
% [~,~,rotateangle] = cov2ell(S);
% rotateangle = modrange(-rotateangle,-pi/2,pi/2);
% % 
% % %find max distance between any pair of positions
% % d = squareform(pdist([xflies(:),yflies(:)]));
% % [maxd,k] = max(d(:));
% % [i,j] = ind2sub([numel(x),numel(x)],k);
% % rotateangle = modrange(-atan2(y(j)-y(i),x(j)-x(i)),-pi/2,pi/2);
% 
% % rotate around the center of the image
% R = [cos(rotateangle),sin(rotateangle)
%   -sin(rotateangle),cos(rotateangle)];
% T = [1,0,0;0,1,0;-imcenter,1]*[R,[0;0];[0,0,1]]*[1,0,0;0,1,0;imcenter,1];
% p0 = [x(:),y(:),ones(numel(x),1)];
% p = p0*T;
% x = reshape(p(:,1),[trx.nflies,nframes_total]);
% y = reshape(p(:,2),[trx.nflies,nframes_total]);
% theta = theta+rotateangle;
% 
% 
% tform = maketform('affine', T);

% p0 = [x(:),y(:),ones(numel(x),1)];
% p1 = p0*[1,0,0;0,1,0;-imcenter,1];
% q1 = bsxfun(@minus,p0,[imcenter,0]);
% max(abs(p1-q1))
% p2 = p1*[R,[0;0];[0,0,1]];
% q2 = [q1(:,1:2)*R,ones(numel(x),1)];
% max(abs(p2-q2))
% p3 = p2*[1,0,0;0,1,0;imcenter,1];
% q3 = bsxfun(@plus,q2,[imcenter,0]);
% max(abs(p3-q3))
% r3 = p0*T;
% max(abs(p3-r3))

% bkgdim = imtransform(bkgdim0,tform,'XData',[1,nc0],'YData',[1,nr0],'XYScale',1);

%% figure out the axes limits
xlim = [inf,-inf];
ylim = [inf,-inf];
for fly = allflies,
  for j = 1:nframes_total,
    [x1,x2,y1,y2] = ellipse_to_bounding_box(x(fly,j),y(fly,j),a(fly,j),b(fly,j),theta(fly,j));
    xlim(1) = min([xlim(1),x1,x2]);
    xlim(2) = max([xlim(2),x1,x2]);
    ylim(1) = min([ylim(1),y1,y2]);
    ylim(2) = max([ylim(2),y1,y2]);
  end
end
xlim = [floor(xlim(1))-border,ceil(xlim(2))+border];
ylim = [floor(ylim(1))-border,ceil(ylim(2))+border];
%bkgdim = bkgdim(ylim(1):ylim(2),xlim(1):xlim(2));
x = x - xlim(1)+1;
y = y - ylim(1)+1;

%% colors for each frame
colors_all = nan(ts(end)-ts(1)+1,3,nflygroups);
for i = 1:nflygroups,
  cmcurr = cm{min(numel(cm),i)};
  if isempty(cmcurr),
    colors_all(:,:,i) = 1;
  elseif isnumeric(cmcurr),
    colors_all(:,:,i) = cmcurr;
  else
    colors_all(:,:,i) = cmcurr(ts(end)-ts(1)+1);
  end
end
colors = colors_all(ts-ts(1)+1,:,:);


%% read in the frames, convert to grayscale
imcurr = im0;
if size(imcurr,3) > 1,
  imcurr = rgb2gray(imcurr);
end
%imcurr = imtransform(imcurr,tform,'XData',[1,nc0],'YData',[1,nr0],'XYScale',1);
imcurr = imcurr(ylim(1):ylim(2),xlim(1):xlim(2));
im = im2double(imcurr);

%% plot

if isempty(hfig),
  hfig = figure;
end

figure(hfig);
clf;
if ~isempty(figpos),
  set(hfig,'Units','pixels','Position',figpos);
end
hax = gca;
imagesc(im,[0,1]);
axis image off;
hold on;
colormap gray;

hell = nan(1,nframes);
for i = 1:nframes,
  j = ts(i)-ts(1)+1;
  hell(i) = drawellipse(x(mainfly,j),y(mainfly,j),theta(mainfly,j),2*a(mainfly,j),2*b(mainfly,j));
  set(hell(i),'LineWidth',2,'Color',colors(i,:,1));
end

linecolors = colors_all(:,:,1)*.75;
hmain = patch([x(mainfly,:),nan]',[y(mainfly,:),nan]',permute([linecolors;0,0,0],[1,3,2]),...
  'EdgeColor','interp','FaceColor','none',...
  'Marker','.','EdgeAlpha',.25,'MarkerFaceColor','k','MarkerEdgeColor',[.25,.25,.25],....
  'LineWidth',1.5);
idxpos = predictions{mainfly}(ts(1):ts(end));
plot(x(mainfly,idxpos),y(mainfly,idxpos),'.','color',colorpos);
plot(x(mainfly,~idxpos),y(mainfly,~idxpos),'.','color',colorneg);

hothers = nan(1,numel(otherflies));
for flyi = 1:numel(otherflies),
  fly = otherflies(flyi);
  %plot(x,y,'.-','color',[.5,.5,.5]);
  linecolors = colors_all(:,:,flyi+1)*.5+.25;
    
  hothers(flyi) = patch([x(fly,:),nan]',[y(fly,:),nan]',permute([linecolors;0,0,0],[1,3,2]),...
  'EdgeColor','interp','FaceColor','none',...
  'Marker','.','EdgeAlpha',.25,'MarkerFaceColor','k','MarkerEdgeColor',[.5,.5,.5],...
  'LineStyle','-','LineWidth',1.5);
%   idxpos = predictions{fly}(ts(1):ts(end));
%   plot(x(fly,idxpos),y(fly,idxpos),'.','color',colorpos);
%   plot(x(fly,~idxpos),y(fly,~idxpos),'.','color',colorneg);
end

% colormap(colors_all(:,:,1)*.75);
% tlim = trx(fly).timestamps(ts(end)+trx(fly).off)-...
%   trx(fly).timestamps(ts(1)+trx(fly).off);
% set(hax,'CLim',[0,tlim]);
% hcb = colorbar;
% tticks = 0:.25:tlim;
% realylim = get(hcb,'YLim');
% yticks = realylim(1)+(tticks/tlim)*(realylim(2)-realylim(1));
% set(hcb,'YTick',yticks,'YTickLabel',tticks');


% 
% 
% axes(hax(2));
% image([0,ts(end)-ts(1)],[1,3],repmat(permute(colors*(1-max_weight_color)+.1,[3,1,2]),[2,1,1]));
% set(hax(2),'Box','off');
% set(hax(2),'XTick',ts-ts(1),'YTick',[]);
