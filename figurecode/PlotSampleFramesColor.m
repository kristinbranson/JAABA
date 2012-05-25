function [hfig,finalim,final_pforefly] = PlotSampleFramesColor(trx,readframe,bkgdim0,predictions,mainfly,otherflies,ts,varargin)

colorpos = [.7,0,0];
colorneg = [0,0,.7];
border = 20;

[colorpos,colorneg,border,hfig,figpos,...
  cm,fg_thresh,bg_thresh,sigma_bkgd,wmah,frac_a_back,dist_epsilon] = ...
  myparse(varargin,'colorpos',colorpos,...
  'colorneg',colorneg,...
  'border',border,...
  'hfig',[],...
  'figpos',[],...
  'colormap',{@jet,@jet},...
  'fg_thresh',115/255,...
  'bg_thresh',10/255,...
  'sigma_bkgd',115/255-10/255,...
  'wmah',.5,...
  'frac_a_back',1,...
  'dist_epsilon',.1);

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

[nr0,nc0,~] = size(bkgdim0);

%% grab out relevant time points
% not rotating currently, maybe rotate in the future?

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

bkgdim = bkgdim0;

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
bkgdim = bkgdim(ylim(1):ylim(2),xlim(1):xlim(2));
x = x - xlim(1)+1;
y = y - ylim(1)+1;

%% colors for each frame
colors_all = nan(ts(end)-ts(1)+1,3,nflygroups);
for i = 1:nflygroups,
  cmcurr = cm{min(numel(cm),i)};
  if isnumeric(cmcurr),
    colors_all(:,:,i) = cmcurr;
  else
    colors_all(:,:,i) = cmcurr(ts(end)-ts(1)+1);
  end
end
colors = colors_all(ts-ts(1)+1,:,:);

%% read in the frames, convert to grayscale
for i = 1:nframes,
  t = ts(i);
  imcurr = readframe(t);
  if size(imcurr,3) > 1,
    imcurr = rgb2gray(imcurr);
  end
  %imcurr = imtransform(imcurr,tform,'XData',[1,nc0],'YData',[1,nr0],'XYScale',1);
  imcurr = imcurr(ylim(1):ylim(2),xlim(1):xlim(2));
  if i == 1,
    gray_ims = repmat(im2double(imcurr),[1,1,nframes]);
  else
    gray_ims(:,:,i) = im2double(imcurr);
  end
end

%% probability that each pixel is foreground
dbkgd = abs(bsxfun(@minus,gray_ims,bkgdim));
pback = exp(-dbkgd/sigma_bkgd^2/2);
pback = pback / max(pback(:));
minv = exp(-fg_thresh/sigma_bkgd^2/2);
maxv = exp(-bg_thresh/sigma_bkgd^2/2);
pback = min(1,max(0,(pback-minv) / (maxv-minv)));
pfore = 1-pback;

%% make a connectivity graph for computing distances within foreground
% pixels
nr = ylim(2)-ylim(1)+1;
nc = xlim(2)-xlim(1)+1;
np = nr*nc;

[cgrid,rgrid] = meshgrid(1:nc,1:nr);
rconn1 = repmat(rgrid,[1,1,4]);
cconn1 = repmat(cgrid,[1,1,4]);
rconn2 = nan(size(rconn1));
cconn2 = nan(size(cconn1));

rconn2(:,:,1) = rgrid+1; % above
cconn2(:,:,1) = cgrid; 
rconn2(:,:,2) = rgrid-1; % below
cconn2(:,:,2) = cgrid;
rconn2(:,:,3) = rgrid; % right
cconn2(:,:,3) = cgrid+1;
rconn2(:,:,4) = rgrid; % left
cconn2(:,:,4) = cgrid-1; 
badidx = rconn2 < 1 | rconn2 > nr | cconn2 < 1 | cconn2 > nc;
rconn1(badidx) = [];
cconn1(badidx) = [];
rconn2(badidx) = [];
cconn2(badidx) = [];
conn1 = sub2ind([nr,nc],rconn1,cconn1);
conn2 = sub2ind([nr,nc],rconn2,cconn2);

%isdy = conn2 == conn1+1 | conn2 == conn1-1;
%w = min(maxdist,-log(pfore)+dist_epsilon)/maxdist;
w = reshape((pback+dist_epsilon)/(1+dist_epsilon),[np,nframes]);

%% assign pixels to flies

% probability that each pixel belongs to each fly
pfly = zeros(nr*nc,trx.nflies,nframes);
%allflies = [flygroups{:}];

for ii = 1:nframes,
  t = ts(ii);
  i = t-ts(1)+1;
  for fly = 1:trx.nflies,
    if t > trx(fly).endframe || t < trx(fly).firstframe,
      continue;
    end
    mu = [x(fly,i)-frac_a_back*a(fly,i)*cos(theta(fly,i)),...
      y(fly,i)-frac_a_back*a(fly,i)*sin(theta(fly,i))];
%     mu = [trx(fly).x(j)-frac_a_back*trx(fly).a(j)*cos(trx(fly).theta(j)),...
%       trx(fly).y(j)-frac_a_back*trx(fly).a(j)*sin(trx(fly).theta(j))];
    
    s = sub2ind([nr,nc],min(nr,max(1,round(mu(2)))),...
      min(nc,max(1,round(mu(1)))));
    wcurr = w(conn2,ii);
    S = axes2cov(a(fly,i),b(fly,i),theta(fly,i));
    diffs = bsxfun(@minus,[cgrid(:),rgrid(:)],mu);
    c = chol(S);
    dmah = sum((diffs/c).^2',1); %#ok<UDIM>
%     ratio = S(1,1)/S(2,2);
%     rationorm = 2 / (ratio+1);
%     wcurr(isdy) = wcurr(isdy)*ratio*rationorm;
%     wcurr(~isdy) = wcurr(~isdy)*rationorm;
    G = sparse(conn1,conn2,wcurr,np,np);
    d = graphshortestpath(G,s,'Directed',true);
    %pfly(:,fly,i) = exp(-d.^2/(trx(fly).a(j)*dist_epsilon)^2);
    pfly(:,fly,ii) = exp(-(1-wmah)*d.^2/(a(fly,i)*dist_epsilon)^2 - wmah*dmah);
  end
end
Z = sum(pfly,2);
Z(Z==0) = 1;
pfly = bsxfun(@rdivide,pfly,Z);

gray_ims = reshape(gray_ims,[nr*nc,nframes]);
pfore = reshape(pfore,[nr*nc,nframes]);
psomefly = reshape(1-prod(1-pfly(:,[flygroups{:}],:),2),[nr*nc,nframes]);
pforefly = pfore.*psomefly;
%pfore_any = 1 - prod(1-pfore.*psomefly,2);

%% color the images

im = zeros(nr*nc,3,nframes);
for i = 1:nframes,
  imcurr = zeros(nr*nc,3);
  for flygroupi = 1:nflygroups,
    flygroup = flygroups{flygroupi};
    pcurr = 1 - prod(1-pfly(:,flygroup,i),2);
    colorcurr = colors(i,:,flygroupi);
    imcurr = imcurr + bsxfun(@times,colorcurr,pcurr);
  end
  imcurr = bsxfun(@plus,bsxfun(@times,imcurr,gray_ims(:,i).*pforefly(:,i)),...
    (1-pforefly(:,i)).*reshape(bkgdim,[nr*nc,1]));
  im(:,:,i) = imcurr;
end
finalim = reshape(im,[nr,nc,3,nframes]);

final_pforefly = reshape(pforefly,[nr,nc,nframes]);


%%

if isempty(hfig),
  hfig = figure;
end

figure(hfig);
clf;
hax = createsubplots(1,numel(ts),.01);

fly = mainfly;
idxpos = predictions{mainfly}(ts(1):ts(end));
timestamps = trx(mainfly).timestamps(ts+trx(mainfly).off)-...
  trx(mainfly).timestamps(ts(1)+trx(mainfly).off);

for i = 1:numel(ts),
  t = ts(i);
  image(finalim(:,:,:,i),'Parent',hax(i));
  axis(hax(i),'image','off');
  hold(hax(i),'on');
  plot(hax(i),x(mainfly,:),y(mainfly,:),'k.-');
  plot(hax(i),x(mainfly,idxpos),y(mainfly,idxpos),'.','Color',colorpos);
  plot(hax(i),x(mainfly,~idxpos),y(mainfly,~idxpos),'.','Color',colorneg);
  if predictions{mainfly}(t),
    colorcurr = colorpos;
  else
    colorcurr = colorneg;
  end
  plot(hax(i),x(mainfly,t-ts(1)+1),y(mainfly,t-ts(1)+1),'o','color',colorcurr,'markerfacecolor',colorcurr);
  text(1,1,sprintf('t = %.2fs',timestamps(i)),'HorizontalAlignment','left','VerticalAlignment','top','Parent',hax(i));
end

axes(hax(end));
colormap(colors_all(:,:,1)*.75);
tlim = timestamps(end);
set(hax,'CLim',[0,tlim]);
hcb = colorbar('East');
tticks = 0:.25:tlim;
realylim = get(hcb,'YLim');
yticks = realylim(1)+(tticks/tlim)*(realylim(2)-realylim(1));
set(hcb,'YTick',yticks,'YTickLabel',tticks');


if isempty(figpos),
  truesize;
else
  set(hfig,'Units','pixels','Position',figpos);
end