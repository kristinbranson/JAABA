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
scale = params.scale;

optflowwinsig = params.optflowwinsig ;
optflowsig = params.optflowsig ;
optreliability = params.optreliability ;

maxflow = 5;
fly_thres = 90;
flow_thres = params.flow_thres;

makeVideo = false;

% Antennal grooming for andy
expdir = '/home/mayank/Dropbox/ForMayankFromAlice/grooming_GMR_30B01_AE_01_CsChr_RigB_20150903T161828/';
fly = 1;
frames = 5961:5980; 

% expdir = '../mated10_20140714T131113';
% % frames = 10900 +(1:200);
% % fly = 1;
% fly = 1;
% % frames = 11043; % large middle leg movement
% frames = 10925; % large middle leg movement

% expdir = '/home/mayank/Work/FlySpaceTime/walkMovies/SS03500';
% fly = 1;
% frames = 10200 + (1:100);
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

maxv2 = 0;
colors = hsv(nbins);
sz = size(im1);
bwimg = zeros(sz(1),sz(2));
ctr = [ceil( (sz(1)+1)/2),ceil( (sz(2)+1)/2)];
bwimg(ctr(1),ctr(2))=1;
dimg = bwdist(bwimg,'euclidean');
[xx,yy]= meshgrid(1:sz(2),1:sz(1));
aimg = atan2(-(yy-ctr(1)),-(xx-ctr(2)));

for t = t0:20:t1,

  im1curr = im1(:,:,t-t0+1);
  im2curr = im2(:,:,t-t0+1);
%   [Vx,Vy,~] = optFlowLk(im1(:,:,i),im2(:,:,i),[],optflowwinsig,optflowsig,optreliability);
%  [Vx,Vy,] = optFlowHorn(im1curr,im2curr,optflowsig);
  uv = estimate_flow_interface(im1curr,im2curr,'hs-brightness',{'max_warping_iters',10});
%   uv = estimate_flow_interface(im1curr,im2curr,'ba-brightness');
  Vx = uv(:,:,1); Vy = uv(:,:,2);

%   Vx = Vx-dx(t-t0+1);
%   Vy = Vy-dy(t-t0+1);
  
  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  H = gradientHist(single(M),single(O),psize,nbins,1);
  maxv2 = max(maxv2,max(H(:)));
end

%%
if makeVideo
  vid = VideoWriter(sprintf('%s_Fly%d_From%d_To%d_HOF_%s.avi',expname,fly,t0,t1,datestr(now,'yyyymmdd')));
  open(vid);
end

dd_err = dimg/patchsz;
% There is some flow towards the center in certain cases.
% d_err is for that.

for t = t0:t1,

  im1curr = im1(:,:,t-t0+1);
  im2curr = im2(:,:,t-t0+1);

  imsz = size(im1curr);
  pairimg = zeros(2*imsz(1),3*imsz(2),3);
  ttimg = [];
  tt(:,:,1) = im2curr;
  tt(:,:,2) = im1curr;
  tt(:,:,3) = im2curr;
  pairimg(1:imsz(1),:,:) = repmat(tt,[1 3 1]);
  
  
%   [Vx,Vy,~] = optFlowLk(im1curr,im2curr,[],3);
%   [Vx,Vy,~] = optFlowLk(im1curr,im2curr,[],optflowwinsig,optflowsig,optreliability);
%   [Vx,Vy,] = optFlowHorn(im1curr,im2curr,optflowsig);
  uv = estimate_flow_interface(im1curr,im2curr,'hs-brightness',...
    {'max_warping_iters',2 });
%   uv = estimate_flow_interface(im1curr,im2curr,'ba-brightness',{'max_warping_iters',2 });
%   uv = estimate_flow_interface(im1curr,im2curr,'classic++');
  uvorig = uv;
  pairimg(imsz(1)+(1:imsz(1)),1:imsz(2),:) = flowToColor(uv,maxflow);
  
  cdx = dx(t-t0+1);
  cdy = dy(t-t0+1);
  curt = theta(t-t0+1);
  rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[cdx;cdy];
  cdx = -rotd(1); cdy = -rotd(2);

  ctheta = dtheta(t-t0+1);
  rotflowu = dimg.*(cos(aimg+ctheta)-cos(aimg));
  rotflowv = dimg.*(sin(aimg+ctheta)-sin(aimg));

  
  uv = uvorig;
  dd1 = sqrt( (uv(:,:,1)-cdx-rotflowu).^2 + (uv(:,:,2)-cdy-rotflowv).^2);
  for ndx = 1:2
    tt = uv(:,:,ndx);
    tt(dd1< (dd_err+flow_thres)) = 0;
    uv(:,:,ndx) = tt;
  end
  uvmotion = uv;
  pairimg(imsz(1)+(1:imsz(1)),imsz(2)+(1:imsz(2)),:) = flowToColor(uvmotion,maxflow);

  % Using bkg/fkg.
%   fly_bod = im1curr<fly_thres;
%   fly_flow = zeros(1,2);
%   for ndx = 1:2
%     tt = uv(:,:,ndx);
%     fly_flow(ndx) = median(tt(fly_bod));
%   end

  %using computation.
  crdx = rdx(t-t0+1);  crdy = rdy(t-t0+1);
  rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[crdx;crdy];
  fly_flow = rotd;
  
  uv = uvorig;
  dd2 = sqrt( (uv(:,:,1)-fly_flow(1)).^2 + (uv(:,:,2)-fly_flow(2)).^2);
  dd = min(dd1,dd2);
  for ndx = 1:2
    tt = uv(:,:,ndx);
    tt( dd<(dd_err+flow_thres)) = 0;
    uv(:,:,ndx) = tt;
  end
  uvall = uv;
  
  
  pairimg(imsz(1)+(1:imsz(1)),2*imsz(2)+(1:imsz(2)),:) = flowToColor(uvall,maxflow);

  pairimg = uint8(pairimg);
  pairimg = imresize(pairimg,scale);
  Vx = uvorig(:,:,1); Vy = uvorig(:,:,2);
  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  Horig = gradientHist(single(M),single(O),psize,nbins,1);

  Vx = uvmotion(:,:,1); Vy = uvmotion(:,:,2);
  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  Hmotion = gradientHist(single(M),single(O),psize,nbins,1);

  Vx = uvall(:,:,1); Vy = uvall(:,:,2);
  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  Hall = gradientHist(single(M),single(O),psize,nbins,1);
  
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
    delete(horig(ishandle(horig)));
    delete(hmotion(ishandle(hmotion)));
    delete(hall(ishandle(hall)));
    set(him,'CData',pairimg);
    set(htext,'String',sprintf('%.2f s',(t-1)/200));
  end

  
  horig = plotHofFeatures(Horig,hax,maxv2,colors,binedges2,psize*scale,[0 0]);
  hmotion = plotHofFeatures(Hmotion,hax,maxv2,colors,binedges2,psize*scale,[imsz(2) 0]*scale);
  hall = plotHofFeatures(Hall,hax,maxv2,colors,binedges2,psize*scale,[2*imsz(2) 0]*scale);
  
  drawnow;
  if makeVideo
    fr = getframe(hax);
    writeVideo(vid,fr);
  elseif t0<t1
    pause;
  end
end

if makeVideo,  close(vid); end
