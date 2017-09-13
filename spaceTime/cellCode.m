
%% histogram of flow values to find the clipping range.
rr = randsample(headerinfo.nframes-1,50);

vals = [];
for ndx = 1:numel(rr)
  I1 = readframefcn(rr(ndx));
  I2 = readframefcn(rr(ndx)+1);
  [Vx,Vy,~] = optFlowLk(I1,I2,[],3);
  xx = [Vx(:) Vy(:)];
  xx(abs(xx)<1) = [];
  vals = [vals xx];
  
end

vals(vals<3) = [];
figure; hist(vals);

%% create an ufmf movie

hh = ufmf_read_header('../mated10_20140714T131113/movie.ufmf');
outfile = '../temp/test.ufmf';
fig = figure('Visible','off');
fid = fopen(outfile,'w');

% Setup the ufmf.
if strcmp(hh.coding,'MONO8'),
  ncolors = 1;
end

gui = get(fig,'UserData');
[im,hh1,timestamp,bb,mu] = ufmf_read_frame(hh,1);

gui.bg.params.UpdatePeriod = inf;
gui.bg.params.NFrames = inf;
gui.bg.params.boxBasedCompression = true;
gui.bg.params.KeyframePeriod = inf;

nroi = 1;
roi = [1 1 hh.max_height hh.max_width];
gui.bg.lastkeyframetime = -inf;
tmp = struct('loc',cast([],'int64'),'timestamp',[]);
index = struct;
index.frame = tmp;
index.keyframe.mean = tmp;
gui.bg.index = repmat(index,[1,nroi]);
% (re)initialize the bg models
if ~isfield(gui.bg,'model') || length(gui.bg.model) ~= nroi,  
  gui.bg.model = struct('roi',cell(1,nroi),'mu',cell(1,nroi),'nframes',cell(1,nroi));
  for i = 1:nroi,
    gui.bg.model(i) = struct('roi',roi(i,:),...
      'mu',mu,...
      'nframes',0);
%      'mu',zeros([roi(i,4),roi(i,3),ncolors],'uint8'),...
  end
end
gui.bg.lastupdatetime = -inf;
gui.isBGModel = true;
set(fig,'UserData',gui);

writeUFMFHeader(fid,fig,1);

% Write few frames
writeUFMFKeyFrame(fig,timestamp,fid,1);

for ndx = 1:10
  [im,hh1,timestamp,bb] = ufmf_read_frame(hh,ndx);
  %bb = bb(:,[2 1 4 3]);
  mywriteUFMFFrame(fig,im,timestamp,fid,1,bb);
  
end

% wrap up the ufmf


wrapupUFMF(fid,fig,1);
fclose(fid);
fclose(hh.fid);
close(fig);
% read the created ufmf

hh1 = ufmf_read_header(outfile);
hh = ufmf_read_header('../mated10_20140714T131113/movie.ufmf');

%
figure;
for ndx = 1:10
  [im1,hh2,timestamp,bb] = ufmf_read_frame(hh1,ndx);
  subplot(1,3,1); imshow(im1);
  [im2,hh2,timestamp,bb] = ufmf_read_frame(hh,ndx);
  subplot(1,3,2); imshow(im2);
  dd = abs(double(im1)-double(im2));
  subplot(1,3,3); imshow(uint8(dd));
  title(sprintf('%d %d',max(dd(:)),nnz(dd(:))));
  pause;
end
fclose all
%%

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn('../mated10_20140714T131113/movie.ufmf');
im = uint8([]);
fstarts= [540,610,96705,97577,98161,105925,106136,106507,106600,115402];
flies = [2 4 9 11 9 11 6 10 6 8];
names = {'front_leg_grooming'
  'wing_grooming'
  'back_leg_grooming'
  'eye_grooming'
  'back_leg_grooming'
  'abdomen_grooming'
  'mid_abdomen_grooming'
  'under_wing_grooming'
  'front_leg_grooming'
  'abdomen_grooming'
  };

for clipno = 1:numel(fstarts)
  fstart = fstarts(clipno);
  flynum = flies(clipno);
  name = names{clipno};
  
  fend = fstart+10;
  count = 1;
  for ndx =fstart:(fend+2)
    im(:,:,count) = readfcn(ndx);
    count = count+1;
  end
  
  % visualize the hog/hof features
  
  I = [];
  
  ftrs = genFeatures(readfcn,headerinfo,fstart,fend,tracks);
  hogftrs = ftrs.hogftrs;
  flowftrs = ftrs.flowftrs;
  
  npatches = 8;
  psize = 10; % patch size for hog/hof
  nbins = 8; % number of bins in hog/hof
  patchsz = psize*npatches;
  
  limshog = prctile(hogftrs{flynum}(:),[1 99]);
  limsflow = prctile(flowftrs{flynum}(:),[1 99]);
  
  for fno = fstart:fend;
    curpatch = [];
    trackndx = fno - tracks(flynum).firstframe + 1;
    locy = round(tracks(flynum).y(trackndx));
    locx = round(tracks(flynum).x(trackndx));
    theta = tracks(flynum).theta(trackndx);
    curpatch(:,:,1) = imsharpen(extractPatch(im(:,:, fno-fstart+1),...
      locy,locx,theta,patchsz));
    % locy = round(tracks(flynum).y(trackndx+1));
    % locx = round(tracks(flynum).x(trackndx+1));
    % theta = tracks(flynum).theta(trackndx+1);
    curpatch(:,:,2) = imsharpen(extractPatch(im(:,:, fno-fstart+2),...
      locy,locx,theta,patchsz));
    
    ftrndx = fno-fstart+1;
    curpatch = imresize(curpatch,2,'bilinear');
    hogI = myhogDraw(hogftrs{flynum}(:,:,:,ftrndx),psize*2,false);
    flowI = myhogDraw(flowftrs{flynum}(:,:,:,ftrndx),psize*2,true);
    
    hogI = (hogI-limshog(1))/(limshog(2)-limshog(1));
    flowI = (flowI-limsflow(1))/(limsflow(2)-limsflow(1));
    
    hogI = uint8(hogI*256);
    flowI =uint8(flowI*256);
    hogI = repmat(hogI,[1 1 3]);
    flowI = repmat(flowI,[1 1 3]);
    curI = repmat(curpatch(:,:,1),[1 1 3]);
    curI = addgrid(curI,[psize psize]*2);
    curdI = imfuse(curpatch(:,:,1),curpatch(:,:,2));
    curdI = addgrid(curdI,[psize psize]*2);
    I(:,:,:,ftrndx,1) = curI;
    I(:,:,:,ftrndx,2) = hogI;
    I(:,:,:,ftrndx,3) = curdI;
    I(:,:,:,ftrndx,4) = flowI;
  end
  
  sz = size(I);
  I = reshape(I,[sz(1) sz(2) sz(3) sz(4)*sz(5)]);
  figure; h = montage2(uint8(I),struct('mm',4,'hasChn',1));
  pI = get(h,'CData');
  tI = cell(1,4);
  count = 1;
  for ndx = 1:size(I,4)
    idx = mod(ndx-1,4)+1;
    tI{idx} = cat(2,tI{idx},I(:,:,:,count));
    count = count + 1;
  end
  imwrite(pI,sprintf('../Figures/mated_frame%d_fly%d_%s.png',fstart,flynum,name))
  % outI = [];
  % for ndx = 1:sz(5)
  %   outI = cat(1,outI,tI{ndx});
  % end
  % figure; imshow(uint8(outI))
end
fclose(fid);


%% compare two opt flow methods side by side

[read1,nframes1] = get_readframe_fcn('HS-brightness_GrabHOF_20150717.avi');
[read2,nframes2] = get_readframe_fcn('HS-brightness_compensated_GrabHOF_20150717.avi');
[read3,nframes3] = get_readframe_fcn('HS-brightness_suppressed_GrabHOF_20150717.avi');
[read4,nframes4] = get_readframe_fcn('LK_GrabHOF_20150717.avi');
assert(nframes1==nframes2);

vid = VideoWriter(sprintf('Allcomparison_%s.avi',datestr(now,'yyyymmdd')));
open(vid);
for ndx = 1:nframes1
  i1 = read1(ndx);
  i2 = read2(ndx);
  i3 = read3(ndx);
  i4 = read4(ndx);
  i2 = imresize(i2,[size(i1,1) size(i1,2)]);
  i3 = imresize(i3,[size(i1,1) size(i1,2)]);
  i4 = imresize(i4,[size(i1,1) size(i1,2)]);
  i1 = insertText(i1,[10 10],'HS');
  i2 = insertText(i2,[10 10],'HS-compensated');
  i3 = insertText(i3,[10 10],'HS-suppressed');
  i4 = insertText(i4,[10 10],'LK');
  allii = [i1 i2; i3 i4];
  writeVideo(vid,allii);
  
end

close(vid);

%% Background Flow compensation

fly_thres = 90;

sz = size(im1);
bwimg = zeros(sz(1),sz(2));
ctr = [round( (sz(1)+1)/2),round( (sz(2)+1)/2)];
bwimg(ctr(1),ctr(2))=1;
dimg = bwdist(bwimg,'euclidean');
[xx,yy]= meshgrid(1:sz(2),1:sz(1));
aimg = atan2(-(yy-ctr(1)),-(xx-ctr(2)));

for t = t0:10:t1;
%%
im1curr = im1(:,:,t-t0+1);
im2curr = im2(:,:,t-t0+1);
fig = figure; subplot(1,4,1);

him = imshowpair(imresize(im1curr,1),imresize(im2curr,1));
ctheta = dtheta(t-t0+1);
title(sprintf('%d:%.2f',t-t0+1,180/pi*ctheta));

uv = estimate_flow_interface(im1curr,im2curr,'hs-brightness');
uvdd = uv;
cdx = dx(t-t0+1);
cdy = dy(t-t0+1);
curt = theta(t-t0+1);
rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[cdx;cdy];
cdx = -rotd(1); cdy = -rotd(2);
ctheta = dtheta(t-t0+1);
uvdd(1:10,1:10,1) = cdx;
uvdd(1:10,1:10,2) = cdy;

rotflowu = (dimg).*(cos(aimg+ctheta)-cos(aimg));
rotflowv = (dimg).*(sin(aimg+ctheta)-sin(aimg));

fly_bod = (im1curr<fly_thres)|(im2curr<fly_thres);
fly_flow = zeros(1,2);
for ndx = 1:2
  tt = uv(:,:,ndx);
  fly_flow(ndx) = median(tt(fly_bod));
end
uvdd((end-9:end)-10,1:10,1) = fly_flow(1);
uvdd((end-9:end)-10,1:10,2) = fly_flow(2);

crdx = rdx(t-t0+1);
crdy = rdy(t-t0+1);
rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[crdx;crdy];
crdx = rotd(1); crdy = rotd(2);
uvdd(end-9:end,1:10,1) = crdx(1);
uvdd(end-9:end,1:10,2) = crdy(1);

figure(fig); subplot(1,4,2);
imshow(flowToColor(uvdd)); 

uvflow = cat(3,rotflowu+cdx,rotflowv+cdy);
% dd = abs(uv(:,:,1)-cdx)+abs(uv(:,:,2)-cdy);
% subplot(1,4,3);
% imagesc(dd); colorbar;axis equal;

dd1 = sqrt( (uv(:,:,1)-cdx-rotflowu).^2 + (uv(:,:,2)-cdy-rotflowv).^2);
dd2 = sqrt( (uv(:,:,1)-fly_flow(1)).^2 + (uv(:,:,2)-fly_flow(2)).^2);
dd = min(dd1,dd2);
subplot(1,4,3);
imshow(flowToColor(uvflow));

uvtt = uv;
for ndx = 1:2
  tt = uvtt(:,:,ndx);
  tt(dd< (sqrt(2)/2)) = 0;
  uvtt(:,:,ndx) = tt;
end
subplot(1,4,4);
imshow(flowToColor(uvtt)); 
%%
pause;

end

%% Debug flow estiamte

nc = 4; nr = 3;

method = 'ba-brightness';

t = t0;
% t = t1-1;
% t = t0+93;
im1curr = im1(:,:,t-t0+1);
im2curr = im2(:,:,t-t0+1);
im3curr = im3(:,:,t-t0+1);
fig = figure; 
subplot(nc,nr,1);
him = imshowpair(imresize(im1curr,1),imresize(im3curr,1));
subplot(nc,nr,2);
him = imshowpair(imresize(im3curr,1),imresize(im2curr,1));
subplot(nc,nr,3);
him = imshowpair(imresize(im1curr,1),imresize(im2curr,1));

cdx = dx(t-t0+1);
cdy = dy(t-t0+1);
curt = theta(t-t0+1);
rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[cdx;cdy];
cdx = -rotd(1); cdy = -rotd(2);

ctheta = dtheta(t-t0+1);
title(sprintf('%d:dx:%.1f,dy:%.1f,theta:%.2f',t-t0+1,cdx,cdy,180/pi*ctheta));

% uv = estimate_flow_interface(im1curr,im2curr,'hs-brightness',{'max_warping_iters',1});
% uv13 = estimate_flow_interface(im1curr,im3curr,'hs-brightness',{'max_warping_iters',1});
% uv23 = estimate_flow_interface(im3curr,im2curr,'hs-brightness',{'max_warping_iters',1});
uv = estimate_flow_interface(im1curr,im2curr,method);
uv13 = estimate_flow_interface(im1curr,im3curr,method);
uv23 = estimate_flow_interface(im3curr,im2curr,method);

figure(fig); subplot(nc,nr,4);
imshow(flowToColor(uv13,maxflow)); 
figure(fig); subplot(nc,nr,5);
imshow(flowToColor(uv23,maxflow)); 
figure(fig); subplot(nc,nr,6);
imshow(flowToColor(uv,maxflow)); 

subplot(nc,nr,7);
flowtrans = cat(3,repmat(cdx,size(im1curr)),repmat(cdy,size(im1curr)));
imshow(flowToColor(flowtrans,maxflow));

subplot(nc,nr,8)
rotflowu = (dimg).*(cos(aimg+ctheta)-cos(aimg));
rotflowv = (dimg).*(sin(aimg+ctheta)-sin(aimg));
flowrot = cat(3,rotflowu,rotflowv);
imshow(flowToColor(flowrot,maxflow));

uvflow = cat(3,rotflowu+cdx,rotflowv+cdy);
subplot(nc,nr,9);
imshow(flowToColor(uvflow,maxflow));

subplot(nc,nr,10)
imshow(flowToColor(uv13-flowtrans));
subplot(nc,nr,11);
imshow(flowToColor(uv23-flowrot));
subplot(nc,nr,12);
imshow(flowToColor(uv-uvflow));

%%

ss = 1./[0.9:0.05:1.1];
figure;count = 1;
for cc = ss(:)'
rotflowu = (dimg).*(cos(aimg+ctheta)-cos(aimg))*cc;
rotflowv = (dimg).*(sin(aimg+ctheta)-sin(aimg))*cc;
flowrot = cat(3,rotflowu,rotflowv);

gg = uv23-flowrot;
subplot(3,numel(ss),count);quiver(gg(:,:,1),gg(:,:,2));
subplot(3,numel(ss),count+2*numel(ss));imshow(flowToColor(gg,4));
count = count + 1;
end

%% Test computeFlow

t = t0 + 4;
% t = t1-1;
% t = t0+93;
im1curr = im1(:,:,t-t0+1);
im2curr = im2(:,:,t-t0+1);

uv = estimate_flow_interface(im1curr,im2curr,'hs-brightness',{'max_warping_iters',1});

cdx = dx(t-t0+1);
cdy = dy(t-t0+1);
curt = theta(t-t0+1);
rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[cdx;cdy];
cdx = -rotd(1); cdy = -rotd(2);

ctheta = dtheta(t-t0+1);
rotflowu = dimg.*(cos(aimg+ctheta)-cos(aimg));
rotflowv = dimg.*(sin(aimg+ctheta)-sin(aimg));

crdx = rdx(t-t0+1);  crdy = rdy(t-t0+1);
rotd = [cos(curt) sin(curt); -sin(curt) cos(curt)]*[crdx;crdy];
fly_flow = rotd;

dd1 = sqrt( (uv(:,:,1)-cdx-rotflowu).^2 + (uv(:,:,2)-cdy-rotflowv).^2);
dd2 = sqrt( (uv(:,:,1)-fly_flow(1)).^2 + (uv(:,:,2)-fly_flow(2)).^2);
dd = min(dd1,dd2);
for ndx = 1:2
  tt = uv(:,:,ndx);
  tt( dd<(dd_err+flow_thres)) = 0;
  uv(:,:,ndx) = tt;
end
uvLocal = uv;


params.dx = cdx; params.dy = cdy;
params.rdx = crdx; params.rdy = crdy;
params.theta = curt; params.dtheta = ctheta;
params.flow_thres = flow_thres; 
params.dd_err = dd_err;
params.dimg = dimg; params.aimg = aimg;

[Vx,Vy] = computeFlowBkgSup(im1curr,im2curr,params);
uvC = cat(3,Vx,Vy);
disp(max(abs(uvLocal(:)-uvC(:))))

%% testing patch match

cores = 24; algo = 'cputiled';
A = im1curr;
B = im2curr;
psz = size(A);
patch_w = 5;
hsz = (patch_w-1)/2;
ann0 = nnmex(repmat(uint8(A),1,1,3), repmat(uint8(B),1,1,3), ...
  algo, patch_w, [], [], [], [], [], cores,[],[15 15]);   % Warm up
dd = zeros(psz(1),psz(2));
dd(hsz+1:end-hsz,hsz+1:end-hsz) = double(ann0(1:end-patch_w+1,1:end-patch_w+1,3));
dd(dd>100000) = 0;
% dd(dd>5000)=5000;
[xx yy] = meshgrid(1:size(A,1),1:size(A,2));
kk = zeros(size(A,1),size(A,2),2);
kk(hsz+1:end-hsz,hsz+1:end-hsz,:) = double(ann0(1:end-patch_w+1,1:end-patch_w+1,1:2));
kk(:,:,1) = kk(:,:,1)-xx;
kk(:,:,2) = kk(:,:,2)-yy;
figure;
ax = [];
ax(1) = subplot(1,3,1); imshowpair(A,B);
ax(2) = subplot(1,3,2); imshow(flowToColor(kk,5));
ax(3) = subplot(1,3,3); imagesc(dd); axis equal; axis image
linkaxes(ax);


%% debug deepmex and offline deepmatching
im1curr = im1(:,:,1);
im2curr = im2(:,:,1);
scale = 4;
tic;
for ndx = 1:4,
[uvo,uvso] = computeDeepFlow(im1curr,im2curr,scale);
end
toc
tic;
for ndx =1:4
k = deepmex(single(im1curr),single(im2curr),80*scale,80*scale,round(80*scale*0.23));
end
toc
A = k';
A(:,1:4) = A(:,1:4)/scale;
Vx = zeros(sz); 
Vy = Vx; Vs = Vx;
idx = sub2ind(sz,round(A(:,1)-0.0006)+1,round(A(:,2)-0.0006)+1);
Vx(idx) = A(:,4)-A(:,2);
Vy(idx) = A(:,3)-A(:,1);
Vs(idx) = A(:,5);
uv = cat(3,Vx,Vy);

nc = 2; nr =1;
figure; 
subplot(nr,nc,1); imshow(flowToColor(uvo));
subplot(nr,nc,2); imshow(flowToColor(uv));


%% Make a short movie to test flow feature computation.

Q = load('../walkMovies/SS03500_test/trx.mat');
selflds = {'x','y','theta','a','b','timestamps',...
  'x_mm','y_mm','theta_mm','a_mm','b_mm'};
maxframe = 10000;
for ndx = 1:numel(Q.trx)
  for fnum = 1:numel(selflds)
    Q.trx(ndx).(selflds{fnum}) = Q.trx(ndx).(selflds{fnum})(1:maxframe);
  end
  Q.trx(ndx).dt = Q.trx(ndx).dt(1:maxframe-1);
  Q.trx(ndx).endframe = maxframe;
  Q.trx(ndx).nframes = maxframe;
end
save('../walkMovies/SS03500_test/trx.mat','-struct','Q');

%%

figure;
ax = [];
nr = 3; nc = 3;
ax(1) = subplot(nr,nc,1); 
imshow(flowToColor(uvdeep,4));
ax(2) = subplot(nr,nc,2); 
imshow(flowToColor(uvorigo,4));
ax(3) = subplot(nr,nc,3); 
imshow(flowToColor(uvorig,4));
jj = cat(3,rotflowu,rotflowv);
jj(:,:,1) = jj(:,:,1)+cdx;
jj(:,:,2) = jj(:,:,2)+cdy;
ax(4) = subplot(nr,nc,4);
imagesc(dd1d); axis image
ax(5) = subplot(nr,nc,5);
imagesc(dd1o); axis image
ax(6) = subplot(nr,nc,6);
imagesc(dd1); axis image
ax(7) = subplot(nr,nc,7);
imshow(flowToColor(uvd,4));
ax(8) = subplot(nr,nc,8);
imshow(flowToColor(uvo,4));
ax(9) = subplot(nr,nc,9);
imshow(flowToColor(uv,4));
linkaxes(ax);


linkaxes(ax);

%% Make an antennal grooming video.

expdir = '/home/mayank/Dropbox/ForMayankFromAlice/grooming_GMR_30B01_AE_01_CsChr_RigB_20150903T161828/';
fly = 1;
frames = 5961:6040; 

hofim = cell(1,numel(frames));
hogim = cell(1,numel(frames));
parfor count = 1:numel(frames)
  fnum = frames(count);
  curhofim = VisualizeFlowFeatures(expdir,fly,fnum,0,'trxfilename','registered_trx.mat','method','hs-sup');
  curhogim = VisualizeHogFeatures(expdir,fly,fnum,'trxfilename','registered_trx.mat');
  hofim{count} = curhofim.cdata;
  hogim{count} = curhogim.cdata;
end

%%
fig = figure;
for ndx = 61:numel(frames)
  figure(fig);
  imshow(uint8([hofim{ndx} hogim{ndx}]));
  title(sprintf('%d',ndx))
  pause(0.5);
end
close(fig)

%%
fig = figure;
ndx = 72;
figure(fig);
imshow(uint8([hofim{ndx} hogim{ndx}]));

%%
expdir = '/home/mayank/Dropbox/ForMayankFromAlice/grooming_GMR_30B01_AE_01_CsChr_RigB_20150903T161828/';
fly = 1;
frames = 5961:6000; 

fnum = frames(1);
im = VisualizeHogFeatures(expdir,1,fnum,'trxfilename','registered_trx.mat');
figure; imshow(im.cdata);
