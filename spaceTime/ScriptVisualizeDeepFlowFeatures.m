% nbins = 8;
% psize = 10;
% npatches = 8;
% optflowwinsig = 3;
% optflowsig = 2;
% optreliability = 1e-4;
% patchsz = psize*npatches;

params = getParams;
npatches = params.npatches;
psize = params.psize;
nbins = params.nbins; 
patchsz = params.patchsz;

optflowwinsig = params.optflowwinsig ;
optflowsig = params.optflowsig ;
optreliability = params.optreliability ;
scale = params.scale;

maxflow = 5;
fly_thres = 90;
flow_thres = 1; 
% error is double when we compute flow in half the size

makeVideo = false;

% expdir = '../mated10_20140714T131113';
% frames = 10900 +(1:200);
% fly = 1;
% fly = 1;
% frames = 11043; % large middle leg movement
% frames = 10925; % large middle leg movement

expdir = '/home/mayank/Work/FlySpaceTime/walkMovies/SS03500';
fly = 1;
% frames = 10200 + (1:100);
% frames = 1:200;
frames = 7320:7330;
% frames = 735:780;
% fly = 3;
% frames = 1886:1888;
% frames = 1477:1477;

% fly = 2; % This one shows the failure of HS-brightness on rear leg.
% frames = 5503;

[~,expname] = fileparts(expdir);
moviefilestr = 'movie.ufmf';
moviefile = fullfile(expdir,moviefilestr);
trxfilestr = 'trx.mat';
trxfile = fullfile(expdir,trxfilestr);


[readframe,nframes] = get_readframe_fcn(moviefile);
td = load(trxfile);
tracks = td.trx;
dx = []; dy = [];
theta = [];
rdx = [];rdy = [];
im1 = []; im2 = []; im3 = [];
for i = 1:numel(frames),
  curf = readframe(frames(i));
  nextf = readframe(frames(i)+1);
  
  trackndx = frames(i) - tracks(fly).firstframe + 1;
  locy = round(tracks(fly).y(trackndx));
  locx = round(tracks(fly).x(trackndx));
  curpatch = extractPatch(curf,...
    locy,locx,tracks(fly).theta(trackndx),patchsz);
  locy = round(tracks(fly).y(trackndx+1));
  locx = round(tracks(fly).x(trackndx+1));
  curpatch2 = extractPatch(nextf,...
    locy,locx,tracks(fly).theta(trackndx+1),patchsz);

  curpatch3 = extractPatch(nextf,...
    locy,locx,tracks(fly).theta(trackndx),patchsz);
  
  im1(:,:,i) = curpatch;
  im2(:,:,i) = curpatch2;
  im3(:,:,i) = curpatch3;
  dx(i) = round(tracks(fly).x(trackndx+1))-round(tracks(fly).x(trackndx));
  dy(i) = round(tracks(fly).y(trackndx+1))-round(tracks(fly).y(trackndx));
  rdx(i) = tracks(fly).x(trackndx+1)-tracks(fly).x(trackndx)-dx(i);
  rdy(i) = tracks(fly).y(trackndx+1)-tracks(fly).y(trackndx)-dy(i);
  dtheta(i) = tracks(fly).theta(trackndx+1)-tracks(fly).theta(trackndx);
  theta(i) = tracks(fly).theta(trackndx);
end



[nr,nc,~] = size(im1);
%
% figure out the histogram bins
tmptheta = (0:179)*pi/180;
res = nan(nbins,numel(tmptheta));
m = single(ones(10,10));
o = single(zeros(10,10));
for i = 1:numel(tmptheta),
%   fprintf('i = %d, theta = %f\n',i,tmptheta(i));
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

% mayank divides by 2
bincenters2 = bincenters*2;
dt = mean(diff(bincenters2));
binedges2 = [bincenters2(1)-dt/2,(bincenters2(1:end-1)+bincenters2(2:end))/2,bincenters2(end)+dt/2];


%% make a video

tr = 191;
hfig = 100;
t0 = frames(1);
t1 = frames(end);
colorpos = [1,0,0];
colorneg = [0,.3,1];

colors = hsv(nbins);
sz = round(size(im1)/2);
bwimg = zeros(sz(1),sz(2));
ctr = [ceil( (sz(1)+1)/2),ceil( (sz(2)+1)/2)];
bwimg(ctr(1),ctr(2))=1;
dimg = bwdist(bwimg,'euclidean')*2;
[xx,yy]= meshgrid(1:sz(2),1:sz(1));
aimg = atan2(-(yy-ctr(1)),-(xx-ctr(2)));

maxv2 = 1.5;
%%
if makeVideo
  vid = VideoWriter(sprintf('%s_Fly%d_From%d_To%d_DeepHOF_%s.avi',expname,fly,t0,t1,datestr(now,'yyyymmdd')));
  open(vid);
end

dd_err = dimg/patchsz*2;
% There is some flow towards the center in certain cases.
% d_err is for that.

for t = t0:t1,

  im1curr = im1(:,:,t-t0+1);
  im2curr = im2(:,:,t-t0+1);

  imsz = size(im1curr);
  pairimg = zeros(3*imsz(1),3*imsz(2),3);
  ttimg = [];
  ttimg(:,:,1) = im2curr;
  ttimg(:,:,2) = im1curr;
  ttimg(:,:,3) = im2curr;
  pairimg(1:imsz(1),:,:) = repmat(ttimg,[1 3 1]);
  
  [uv,Vs] = computeDeepFlow(im1curr,im2curr);
  [xx,yy] = meshgrid((-5:5),(-5:5));
  
  uv = uv(2:2:end,2:2:end,:);
  Vss = Vs(2:2:end,2:2:end);
  selpx = Vss<3.5;
  im1sm = imresize(im1curr,0.5);
  im2sm = imresize(im2curr,0.5);
  
  uvdeep = uv;
  uvold = estimate_flow_interface(im1sm,im2sm,'hs-brightness',...
    {'max_warping_iters',2});
  uvold = uvold*2;
  for ndx = 1:2
    tt = uv(:,:,ndx);
    ttold = uvold(:,:,ndx);
    tt(selpx) = ttold(selpx);
    uv(:,:,ndx) = tt;
  end
  uvorig = uv;
  uvorigo = uvold;
  
  % BKG flow
  cdx = dx(t-t0+1);
  cdy = dy(t-t0+1);
  curt = theta(t-t0+1);
  rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[cdx;cdy];
  cdx = -rotd(1); cdy = -rotd(2);
  ctheta = dtheta(t-t0+1);
  rotflowu = dimg.*(cos(aimg+ctheta)-cos(aimg));
  rotflowv = dimg.*(sin(aimg+ctheta)-sin(aimg));
  uv = uvorig;
  uvo = uvorigo;
  uvd = uvdeep;
  dd1 = sqrt( (uv(:,:,1)-cdx-rotflowu).^2 + (uv(:,:,2)-cdy-rotflowv).^2);
  dd1o = sqrt( (uvo(:,:,1)-cdx-rotflowu).^2 + (uvo(:,:,2)-cdy-rotflowv).^2);
  dd1d = sqrt( (uvd(:,:,1)-cdx-rotflowu).^2 + (uvd(:,:,2)-cdy-rotflowv).^2);

  %FLY flow
  crdx = rdx(t-t0+1);  crdy = rdy(t-t0+1);
  rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[crdx;crdy];
  fly_flow = rotd;
  
  uv = uvorig; uvo = uvorigo; uvd = uvdeep;
  dd2 = sqrt( (uv(:,:,1)-fly_flow(1)).^2 + (uv(:,:,2)-fly_flow(2)).^2);
  dd2o = sqrt( (uvo(:,:,1)-fly_flow(1)).^2 + (uvo(:,:,2)-fly_flow(2)).^2);
  dd2d = sqrt( (uvd(:,:,1)-fly_flow(1)).^2 + (uvd(:,:,2)-fly_flow(2)).^2);
  dd = min(dd1,dd2);ddo = min(dd1o,dd2o); ddd = min(dd1d,dd2d);
  for ndx = 1:2
    tt = uv(:,:,ndx);
    tt( dd<(dd_err+flow_thres)) = 0;
    uv(:,:,ndx) = tt;
    tt = uvo(:,:,ndx);
    tt( ddo<(dd_err+flow_thres)) = 0;
    uvo(:,:,ndx) = tt;
    tt = uvd(:,:,ndx);
    tt( ddd<(dd_err+flow_thres)) = 0;
    uvd(:,:,ndx) = tt;
  end
  uvall = imresize(uv,2); 
  uvallo = imresize(uvo,2);
  uvalld = imresize(uvd,2);
  
  pairimg(imsz(1)+(1:imsz(1)),(1:imsz(2)),:) = flowToColor(uvall,maxflow);
  pairimg(imsz(1)+(1:imsz(1)),imsz(2)+(1:imsz(2)),:) = flowToColor(uvalld,maxflow);
  pairimg(imsz(1)+(1:imsz(1)),2*imsz(2)+(1:imsz(2)),:) = flowToColor(uvallo,maxflow);

  pairimg(2*imsz(1)+(1:imsz(1)),(1:imsz(2)),:) = flowToColor(imresize(uvorig,2),maxflow);
  pairimg(2*imsz(1)+(1:imsz(1)),imsz(2)+(1:imsz(2)),:) = flowToColor(imresize(uvdeep,2),maxflow);
  pairimg(2*imsz(1)+(1:imsz(1)),2*imsz(2)+(1:imsz(2)),:) = flowToColor(imresize(uvorigo,2),maxflow);

  pairimg = uint8(pairimg);
  pairimg = imresize(pairimg,scale);

  Vx = uvall(:,:,1); Vy = uvall(:,:,2);
  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  Hall = gradientHist(single(M),single(O),psize,nbins,1);
  Vx = uvallo(:,:,1); Vy = uvallo(:,:,2);
  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  Hallo = gradientHist(single(M),single(O),psize,nbins,1);
  Vx = uvalld(:,:,1); Vy = uvalld(:,:,2);
  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  Halld = gradientHist(single(M),single(O),psize,nbins,1);
  
  if t == t0,
    figure(hfig);
    clf;
    hax = axes('Position',[0,0,1,1]);
    set(hfig,'Units','pixels','Position',get(0,'ScreenSize'));

    him = imshow(pairimg);
    axis image;
    truesize;
    colormap gray;
    hold on;
    axis off;
    htext = text(1,nr-1,sprintf('%.2f s',(t-1)/200),'Color','k','HorizontalAlignment','left','VerticalAlignment','bottom','FontSize',18);

  else
    delete(hall(ishandle(hall)));
    delete(hallo(ishandle(hallo)));
    delete(halld(ishandle(halld)));
    set(him,'CData',pairimg);
    set(htext,'String',sprintf('%.2f s',(t-1)/200));
  end

  
  hall = plotHofFeatures(Hall,hax,maxv2,colors,binedges2,psize*scale,[0 0]*scale);
  halld = plotHofFeatures(Halld,hax,maxv2,colors,binedges2,psize*scale,[imsz(2) 0]*scale);
  hallo = plotHofFeatures(Hallo,hax,maxv2,colors,binedges2,psize*scale,[2*imsz(2) 0]*scale);
  
  drawnow;
  if makeVideo
    fr = getframe(hax);
    writeVideo(vid,fr);
  elseif t0<t1
    pause;
  end
end

if makeVideo,  close(vid); end
