% nbins = 8;
% psize = 10;
% npatches = 8;
% optflowwinsig = 3;
% optflowsig = 2;
% optreliability = 1e-4;
% patchsz = psize*npatches;

expdir = '/groups/branson/home/kabram/bransonlab/data/katie/20170308_133311/20170308_133311_08';
moviefilestr = 'movie.avi';
frames = 2831:2900; 
vidoutfile = sprintf('%s_From%d_To%d_HOF_%s.avi',expname,fly,t0,t1,datestr(now,'yyyymmdd'));
makeVideo = true;


%%
params = getParamsKatie;
npatchesx = params.npatchesx;
npatchesy = params.npatchesy;
psize = params.psize;
nbins = params.nbins; 
scale = params.scale;

optflowwinsig = params.optflowwinsig ;
optflowsig = params.optflowsig ;
optreliability = params.optreliability ;

maxflow = 5;
fly_thres = 90;
flow_thres = params.flow_thres;


% Antennal grooming for andy


[~,expname] = fileparts(expdir);
moviefile = fullfile(expdir,moviefilestr);


[readframe,nframes] = get_readframe_fcn(moviefile);
for i = 1:numel(frames)
  curf = readframe(frames(i));
  nextf = readframe(frames(i)+1);
  
  curpatch = curf;
  curpatch2 = nextf;

  im1(:,:,i) = curpatch;
  im2(:,:,i) = curpatch2;
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

  im1curr = gaussSmooth(im1(:,:,t-t0+1),1,'same');
  im2curr = gaussSmooth(im2(:,:,t-t0+1),1,'same');
  uv = estimate_flow_interface(im1curr,im2curr,'hs-brightness',{'max_warping_iters',2});
  Vx = uv(:,:,1); Vy = uv(:,:,2);

  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  H = gradientHist(single(M),single(O),psize,nbins,1);
  maxv2 = max(maxv2,max(H(:)));
end

%%

maxv2 = 4;

if makeVideo
  vid = VideoWriter(vidoutfile);
  open(vid);
end


for t = t0:t1

  im1curr = gaussSmooth(im1(:,:,t-t0+1),1,'same');
  im2curr = gaussSmooth(im2(:,:,t-t0+1),1,'same');

  imsz = size(im1curr);
  pairimg = zeros(2*imsz(1),imsz(2),3);
  ttimg = [];
  tt(:,:,1) = im2curr;
  tt(:,:,2) = im1curr;
  tt(:,:,3) = im2curr;
  pairimg(1:imsz(1),:,:) = tt;
  
  uv = estimate_flow_interface(im1curr,im2curr,'hs-brightness',...
    {'max_warping_iters',2 });
  uvorig = uv;
  pairimg(imsz(1)+(1:imsz(1)),1:imsz(2),:) = flowToColor(uv,maxflow);
  
  pairimg = uint8(pairimg);
  pairimg = imresize(pairimg,scale);
  Vx = uvorig(:,:,1); Vy = uvorig(:,:,2);
  M = sqrt(Vx.^2 + Vy.^2);
  O = mod(atan2(Vy,Vx)/2,pi);
  O = min(O,pi-1e-6);
  Horig = gradientHist(single(M),single(O),psize,nbins,1);
  
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
    set(him,'CData',pairimg);
    set(htext,'String',sprintf('%.2f s',(t-1)/200));
  end

  
  horig = plotHofFeatures(Horig,hax,maxv2,colors,binedges2,psize*scale,[0 0]);
  
  drawnow;
  if makeVideo
    fr = getframe(hax);
    writeVideo(vid,fr);
  elseif t0<t1
    pause;
  end
end

if makeVideo,  close(vid); end
