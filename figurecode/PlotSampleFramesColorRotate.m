function [hfig,finalim,final_pforefly] = PlotSampleFramesColorRotate(trx,readframe,bkgdim0,predictions,mainfly,otherflies,ts,varargin)

colorpos = [.7,0,0];
colorneg = [0,0,.7];
border = 20;

[colorpos,colorneg,border,hfig,figpos,...
  cm,fg_thresh,bg_thresh,sigma_bkgd,wmah,frac_a_back,dist_epsilon,...
  ncolors_reference,solidbkgdcolor,...
  flipim,intensityoff] = ...
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
  'dist_epsilon',.1,...
  'ncolors_reference',[],...
  'solidbkgdcolor',[.85,.85,.85],...
  'flipim',false,...
  'intensityoff',0);

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

[nr0,nc0,ncolors] = size(bkgdim0);

%% grab out relevant time points
% not rotating currently, maybe rotate in the future?

x = nan(trx.nflies,nframes_total);
y = nan(trx.nflies,nframes_total);
a = nan(trx.nflies,nframes_total);
b = nan(trx.nflies,nframes_total);
theta = nan(trx.nflies,nframes_total);
imcenter = [(nc0+1)/2,(nr0+1)/2]; %#ok<NASGU>
for fly = 1:trx.nflies,
  idx = allts+trx(fly).off;
  i0 = find(idx>=1,1);
  i1 = find(idx<=trx(fly).nframes,1,'last');
  if isempty(i0) || isempty(i1),
    continue;
  end
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
bkgdim = bkgdim(ylim(1):ylim(2),xlim(1):xlim(2),:);
x = x - xlim(1)+1;
y = y - ylim(1)+1;

%% figure out how to rotate

imid = ceil(numel(ts)/2);
% rotate by orientation of fly in middle frame
rotateby = -theta(mainfly,ts(imid)-ts(1)+1);
rotatearoundx = x(mainfly,ts(imid)-ts(1)+1);
rotatearoundy = y(mainfly,ts(imid)-ts(1)+1);

T1 = [1,0,0;0,1,0;-rotatearoundx,-rotatearoundy,1];
T2 = [1,0,0;0,1,0;rotatearoundx,rotatearoundy,1];
R = [cos(rotateby),sin(rotateby),0;-sin(rotateby),cos(rotateby),0;0,0,1];
A = T1*R*T2;
tform = maketform('affine',A);

tmp = [x(:),y(:),ones(numel(x),1)]*A;
xr = reshape(tmp(:,1),size(x));
yr = reshape(tmp(:,2),size(x));
thetar = theta - rotateby;

xlimr = [inf,-inf];
ylimr = [inf,-inf];

for fly = allflies,
  for j = 1:nframes_total,
    [x1,x2,y1,y2] = ellipse_to_bounding_box(xr(fly,j),yr(fly,j),a(fly,j),b(fly,j),thetar(fly,j));
    xlimr(1) = min([xlimr(1),x1,x2]);
    xlimr(2) = max([xlimr(2),x1,x2]);
    ylimr(1) = min([ylimr(1),y1,y2]);
    ylimr(2) = max([ylimr(2),y1,y2]);
  end
end
xlimr = [floor(xlimr(1))-border,ceil(xlimr(2))+border];
ylimr = [floor(ylimr(1))-border,ceil(ylimr(2))+border];

xr = xr - xlimr(1)+1;
yr = yr - ylimr(1)+1;

%% colors for each frame
colors_all = nan(ts(end)-ts(1)+1,3,nflygroups);
for i = 1:nflygroups,
  cmcurr = cm{min(numel(cm),i)};
  if isempty(cmcurr),
    colors_all(:,:,i) = 1;
  elseif isnumeric(cmcurr),
    colors_all(:,:,i) = cmcurr;
  elseif ~isempty(ncolors_reference),
    colors_reference = cmcurr(ncolors_reference);
    idxcurr = round(linspace(1,ncolors_reference,ts(end)-ts(1)+1));
    colors_all(:,:,i) = colors_reference(idxcurr,:);
  else
    colors_all(:,:,i) = cmcurr(ts(end)-ts(1)+1);
  end
end
colors = colors_all(ts-ts(1)+1,:,:);

%% read in the frames, convert to grayscale
for i = 1:nframes,
  t = ts(i);
  imcurr = readframe(t);
  imcurr = im2double(imcurr(ylim(1):ylim(2),xlim(1):xlim(2),:));
  if ncolors > 1,
    grayim_curr = rgb2gray(imcurr);
  else
    grayim_curr = imcurr;
    imcurr = repmat(imcurr,[1,1,3]);
  end
  if i == 1,
    ims = repmat(imcurr,[1,1,1,nframes]);
    gray_ims = repmat(grayim_curr,[1,1,nframes]);
  else
    ims(:,:,:,i) = imcurr;
    gray_ims(:,:,i) = grayim_curr;
  end
end

%% probability that each pixel is foreground
dbkgd = permute(sum(abs(bsxfun(@minus,ims,bkgdim)),3),[1,2,4,3])/3;
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
if flipim,
  if intensityoff > 0,
    gray_ims = intensityoff+(1 - gray_ims)*(1-intensityoff);
  end
end
for i = 1:nframes,
  imcurr = zeros(nr*nc,3);
  for flygroupi = 1:nflygroups,
    flygroup = flygroups{flygroupi};
    pcurr = 1 - prod(1-pfly(:,flygroup,i),2);
    colorcurr = colors(i,:,flygroupi);
    imcurr = imcurr + bsxfun(@times,colorcurr,pcurr);
  end
  imcurr = bsxfun(@times,imcurr,gray_ims(:,i).*pforefly(:,i)) + ...
    bsxfun(@times,1-pforefly(:,i),repmat(reshape(solidbkgdcolor,[1,3]),[nr*nc,1]));
  im(:,:,i) = imcurr;
end
finalim = reshape(im,[nr,nc,3,nframes]);

final_pforefly = reshape(pforefly,[nr,nc,nframes]);

%% rotate images

for i = 1:size(finalim,4),
  [imrotcurr,xdata,ydata] = imtransform(finalim(:,:,:,i),tform,'bicubic','FillValues',solidbkgdcolor','xdata',xlimr,'ydata',ylimr);
  if i == 1,
    finalimr = repmat(imrotcurr,[1,1,1,size(finalim,4)]);
  else
    finalimr(:,:,:,i) = imrotcurr;
  end
  if ~all(xdata == xlimr) || ~all(ydata == ylimr),
    error('data and limr do not match');
  end
end

finalimr = min(1,max(0,finalimr));

%%

if isempty(hfig),
  hfig = figure;
end

figure(hfig);
clf;
hax = createsubplots(1,numel(ts),.01);

idxpos = predictions{mainfly}(ts(1):ts(end))>0;
idxneg = predictions{mainfly}(ts(1):ts(end))==0;
idxunlabeled = predictions{mainfly}(ts(1):ts(end))<0;

try
  timestamps = trx(mainfly).timestamps(ts+trx(mainfly).off)-...
    trx(mainfly).timestamps(ts(imid)+trx(mainfly).off);
catch
  timestamps = [0,cumsum(trx(mainfly).dt)];
  timestamps = timestamps(ts+trx(mainfly).off)-...
    timestamps(ts(imid)+trx(mainfly).off);
end

for i = 1:numel(ts),
  t = ts(i);
  image(finalimr(:,:,:,i),'Parent',hax(i));
  axis(hax(i),'image','off');
  hold(hax(i),'on');
  plot(hax(i),xr(mainfly,:),yr(mainfly,:),'k.-','MarkerSize',20);
  plot(hax(i),xr(mainfly,idxunlabeled),yr(mainfly,idxunlabeled),'.','Color','w','MarkerSize',20);
  plot(hax(i),xr(mainfly,idxneg),yr(mainfly,idxneg),'.','Color',colorneg,'MarkerSize',20);
  plot(hax(i),xr(mainfly,idxpos),yr(mainfly,idxpos),'.','Color',colorpos,'MarkerSize',20);
  if predictions{mainfly}(t)>0
    colorcurr = colorpos;
  elseif predictions{mainfly}(t)==0,
    colorcurr = colorneg;
  else
    colorcurr = [1,1,1];
  end
  plot(hax(i),xr(mainfly,t-ts(1)+1),yr(mainfly,t-ts(1)+1),'o','color',colorcurr,'markerfacecolor',colorcurr,'MarkerSize',12);
  text(1,1,sprintf('t = %.2fs',timestamps(i)),'HorizontalAlignment','left','VerticalAlignment','top','Parent',hax(i),'FontSize',12);
end

axes(hax(end));
colormap(colors_all(:,:,1)*.75);
tlim = timestamps(end)-timestamps(1);
set(hax,'CLim',[0,tlim]);
hcb = colorbar('East');
tticks = [0,tlim];
realylim = get(hcb,'YLim');
yticks = realylim(1)+(tticks/tlim)*(realylim(2)-realylim(1));
set(hcb,'YTick',yticks,'YTickLabel',tticks');


if isempty(figpos),
  truesize;
else
  set(hfig,'Units','pixels','Position',figpos);
end