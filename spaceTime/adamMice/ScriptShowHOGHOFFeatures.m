setuppaths;
nbins = 8;
psize = 40;
expdir = '/groups/hantman/hantmanlab/Jay/videos/M119_CNO_G6/2nd_T/20140409/M119_20140409_v01';
%expdir = '/groups/hantman/hantmanlab/Jay/videos/M119_CNO_G6/20140326/M119_20140326_v01';
%frames = [293,314,339,347,357];
frames = [321,348,362,371,375];
scoresfilestr = 'scores_Grabm119.mat';
scoresfile = fullfile(expdir,scoresfilestr);
moviefilestr = 'movie_comb.avi';
moviefile = fullfile(expdir,moviefilestr);
trxfilestr = 'trx.mat';
trxfile = fullfile(expdir,trxfilestr);



[readframe,nframes] = get_readframe_fcn(moviefile);
sd = load(scoresfile);
td = load(trxfile);

for i = 1:numel(frames),
  
  im1(:,:,i) = rgb2gray(readframe(frames(i)-1));
  im2(:,:,i) = rgb2gray(readframe(frames(i)));
    
end

[nr,nc,~] = size(im1);

% figure out the histogram bins
tmptheta = (0:179)*pi/180;
res = nan(nbins,numel(tmptheta));
m = single(ones(10,10));
o = single(zeros(10,10));
for i = 1:numel(tmptheta),
  fprintf('i = %d, theta = %f\n',i,tmptheta(i));
  o(:) = single(tmptheta(i));
  rescurr = gradientHist(m,o,1,nbins,1);
  res(:,i) = rescurr(1,1,:);
end

bincenters = nan(1,nbins);
for i = 1:nbins,
  bincenters(i) = tmptheta(argmax(res(i,:)));
end

% this seems to be what the centers correspond to
bincenters = linspace(0,pi,nbins+1);
bincenters = bincenters(1:nbins);

for i = 1:numel(frames),
  [Vx,Vy,~] = optFlowLk(im1(:,:,i),im2(:,:,i),[],3);

  M = sqrt(Vx.^2 + Vy.^2);
  [x,y] = meshgrid(1:nc,1:nr);
  
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  [~,thetaidx2bin] = max(res,[],1);
  [~,thetaidx] = min(abs(bsxfun(@minus,O(:),tmptheta)),[],2);
  thetaidx = reshape(thetaidx,size(O));
  binidx = thetaidx2bin(thetaidx);
  
  minm = prctile(M(:),90);
  idxplot = M>=minm;
  
  colors = hsv(nbins);
  figure(i);
  clf;
  imagesc(im1(:,:,i),[0,255]); axis image;
  colormap gray;
  hold on;
  for j = 1:nbins,
    
    idxcurr = idxplot & binidx == j;
    h = quiver(x(idxcurr),y(idxcurr),Vx(idxcurr),Vy(idxcurr),0,'Color',colors(j,:));
    
  end
  
end

% mayank divides by 2
bincenters2 = bincenters*2;
dt = mean(diff(bincenters2));
binedges2 = [bincenters2(1)-dt/2,(bincenters2(1:end-1)+bincenters2(2:end))/2,bincenters2(end)+dt/2];

%%

maxv = 0;
for imi = 1:numel(frames),

  [Vx,Vy,~] = optFlowLk(im1(:,:,imi),im2(:,:,imi),[],3);
  M = sqrt(Vx.^2 + Vy.^2);
  [x,y] = meshgrid(1:nc,1:nr);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  H = gradientHist(single(M),single(O),psize,nbins,1);
  maxv = max(maxv,max(H(:)));
  
end


%%
hfig0 = 10;
for imi = 1:numel(frames),

  hfig = hfig0 + imi;
  figure(hfig);
  clf;
  set(hfig,'Units','pixels','Position',get(0,'ScreenSize'));

  imagesc(im1(:,:,imi));
  axis image;
  colormap gray;
  hold on;

  [Vx,Vy,~] = optFlowLk(im1(:,:,imi),im2(:,:,imi),[],3);
  M = sqrt(Vx.^2 + Vy.^2);
  [x,y] = meshgrid(1:nc,1:nr);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  H = gradientHist(single(M),single(O),psize,nbins,1);
    
  
  for xi = 1:ceil(nc/psize),
    cx = psize/2 + 1 + (xi-1)*psize;
    if cx+psize/2 > nc,
      break;
    end
    for yi = 1:ceil(nr/psize),
      cy = psize/2 + 1 + (yi-1)*psize;
      if cy+psize/2 > nr,
        break;
      end
      
      for bini = 1:nbins,
        tmp = linspace(binedges2(bini),binedges2(bini+1),20);
        xcurr = cx + [0,psize/2*cos(tmp),0];
        ycurr = cy + [0,psize/2*sin(tmp),0];
        h(yi,xi,bini) = patch(xcurr,ycurr,colors(bini,:),'LineStyle','none','FaceAlpha',H(yi,xi,bini)/maxv);
      end
      
    end
  end

  fprintf('t = %.2f s, prediction = %d\n',(frames(imi)-frames(1))/200,sd.allScores.postprocessed{1}(frames(imi)));
  axis off;
  
  drawnow;
  
  %saveas(hfig,sprintf('GrabOverlayHOF_%05d.png',frames(imi)),'png');
  
end
  

%% make a video
hfig = 100;
t0 = 201;
t1 = 500;
colorpos = [1,0,0];
colorneg = [0,.3,1];
scores = max(-1,min(1,sd.allScores.scores{1}/sd.allScores.scoreNorm));

maxv2 = maxv;
for t = t0:20:t1,

  im1curr = rgb2gray(readframe(t-1));
  im2curr = rgb2gray(readframe(t));

  [Vx,Vy,~] = optFlowLk(im1curr,im2curr,[],3);
  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  H = gradientHist(single(M),single(O),psize,nbins,1);
  maxv2 = max(maxv2,max(H(:)));
end

vid = VideoWriter(sprintf('GrabHOF_%s.avi',datestr(now,'yyyymmdd')));
open(vid);

for t = t0:t1,

  im1curr = rgb2gray(readframe(t-1));
  im2curr = rgb2gray(readframe(t));
  pred = sd.allScores.postprocessed{1}(t);
  score = scores(t);
  if t == t0,
    figure(hfig);
    clf;
    hax = axes('Position',[0,0,1,1]);
    set(hfig,'Units','pixels','Position',get(0,'ScreenSize'));

    him = imagesc(im1curr);
    axis image;
    truesize;
    colormap gray;
    hold on;
    axis off;
    hpred = plot(15,15,'o','MarkerFaceColor',colorneg,'Color',colorneg,'MarkerSize',20);
    hpredtext = text(33,18,'Grab','Color','w','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',18);
    htext = text(1,nr-1,sprintf('%.2f s',(t-1)/200),'Color','w','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',18);

  else
    delete(h(ishandle(h)));
    set(him,'CData',im1curr);
    set(htext,'String',sprintf('%.2f s',(t-1)/200));
  end

  if pred > 0,
    set(hpred,'Color',colorpos,'MarkerFaceColor',max(0,score)*colorpos);
    set(hpredtext,'String','Grab');
  else
    set(hpred,'Color',colorneg,'MarkerFaceColor',max(0,-score)*colorneg);
    set(hpredtext,'String','');
  end
  [Vx,Vy,~] = optFlowLk(im1curr,im2curr,[],3);
  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  H = gradientHist(single(M),single(O),psize,nbins,1);

  h = [];
  for xi = 1:ceil(nc/psize),
    cx = psize/2 + 1 + (xi-1)*psize;
    if cx+psize/2 > nc,
      break;
    end
    for yi = 1:ceil(nr/psize),
      cy = psize/2 + 1 + (yi-1)*psize;
      if cy+psize/2 > nr,
        break;
      end
      
      for bini = 1:nbins,
        tmp = linspace(binedges2(bini),binedges2(bini+1),20);
        xcurr = cx + [0,psize/2*cos(tmp),0];
        ycurr = cy + [0,psize/2*sin(tmp),0];
        h(yi,xi,bini) = patch(xcurr,ycurr,colors(bini,:),'LineStyle','none','FaceAlpha',min(1,H(yi,xi,bini)/maxv2));
      end
      
    end
  end
  
  drawnow;
  fr = getframe(hax);
  writeVideo(vid,fr);
  
end
close(vid);