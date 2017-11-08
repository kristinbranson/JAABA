%% read stipdet output

fid = fopen('../../temp/out.txt','r');
fgetl(fid);
fgetl(fid);
line = fgetl(fid);
xcords = [];
ycords =[];
ts = [];
sigmas = [];
taus = [];
while ischar(line)
  x = sscanf(line,'%f');
  xcords(end+1) = x(6);
  ycords(end+1) = x(5);
  ts(end+1) = x(7);
  sigmas(end+1) = x(8);
  taus(end+1) = x(9);
  
  line = fgetl(fid);
end

fclose(fid);

%% Generate an input file for stipdet.

infile = '../../videos/6-6-2013SC2LaserOn4-Basler_piA640.avi';

[din,fin,ein] = fileparts(infile);
outfile = fullfile(din, [fin '_pts.txt']);
if ~exist('get_readframe_fcn','builtin');
  addpath('~/JAABA/misc');
end

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(infile);

height = headerinfo.Height;
width = headerinfo.Width;


sigmas = 2.^(3:5);
taus = [2 4];

tskip = 4;

fout = fopen(outfile,'w');

for ss = sigmas(:)'
  psize = round(2*3*3*sqrt(ss));
  sskip = round(psize/2);
  for fno = 2*tskip+1:tskip:nframes
    for y = psize+1:sskip:height-psize
      for x = psize+1:sskip:width-psize
        for curtau = taus(:)'
          fprintf(fout,'# 200 %d %d %d %d %d 1\n',x,y,fno,ss,curtau);
        end
        
      end
    end
  end
  
  
end


fclose(fout);

%%

addpath(genpath('toolbox'));
addpath('~/bransonlab/mayank/JAABA/filehandling/');
addpath(genpath('OpticalFlow'));
addpath(genpath('pdollarOF'));
addpath('~//bransonlab/mayank/JAABA/filehandling/');
addpath('~//bransonlab/mayank/JAABA/misc/');
addpath('~//bransonlab/mayank/JAABA/perframe/');
%%

infile = '../../videos/6-6-2013SC2LaserOn4-Basler_piA640.avi';

[din,fin,ein] = fileparts(infile);
outfile = fullfile(din, [fin '_pts.txt']);

[readframe,nframes,fid,headerinfo] = get_readframe_fcn(infile);

imMat = zeros(headerinfo.Height,headerinfo.Width,nframes);
for ndx = 1:nframes
  imMat(:,:,ndx) = rgb2gray(readframe(ndx));
end


%%

f1 = 180;
f2 = f1+1;
I1 = imMat(:,:,f1);
I2 = imMat(:,:,f2);
[M,O] = gradientMag(single(I1/255)); H1=gradientHist(M,O,8,4,1); H1(H1>0.03) = 0.03;
[M,O] = gradientMag(single(I2/255)); H2=gradientHist(M,O,8,4,1); H2(H2>0.03) = 0.03;

H1 = H1 + 0.0001*randn(size(H1));
H2 = H2 + 0.0001*randn(size(H2));
% ff = figure; s = montage2(H1);
% ff = figure; s = montage2(H2);


%%

[Vx,Vy,zz] = optFlowLk(I1,I2,[],3);
flow(:,:,1) = Vx; flow(:,:,2) = Vy;
% figure; imshow(flowToColor(flow));
M = sqrt(Vx.^2 + Vy.^2);
O = mod(atan2(Vy,Vx),pi);
H = gradientHist(single(M),single(O),30,8,1);
figure; montage2(H,prm);

% hold on; quiver(Vx,Vy);


%%
ff = figure; s = montage2(H1);
H1I = get(s,'CDATA'); close(ff);
ff = figure; s = montage2(H2);
H2I = get(s,'CDATA'); close(ff);
gifH = cat(4,H1I,H2I);
gifH = repmat(gifH,[1 1 3 1]);
gifH = gifH./max(gifH(:));
mm = immovie(gifH);
implay(mm);

%%


gifI = uint8(cat(4,I1,I2));
gifI = repmat(gifI,[1 1 3 1]);
frame2gif(gifI,fullfile('temp.gif'));

%%
off = 3;
dir = 1;
del = 0.0001;
clip = 1;
figure;
for dir = 1:4
switch dir
  case 1
    D = (H1(off:end,:,:)+H2(1:end-off+1,:,:)+del);
  case 2
    D = (H1(:,off:end,:)+H2(:,1:end-off+1,:)+del);
  case 3
    D = (H1(off:end,off:end,:)+H2(1:end-off+1,1:end-off+1,:)+del);
  case 4
    D = (H1(1:end-off+1,off:end,:)+H2(off:end,1:end-off+1,:)+del);
end
D = abs(D);

D(D>clip) = clip;
subplot(2,2,dir); montage2(D);
end

%% 

tic;
[F,H] = genFeatures(infile);
toc;

%%
prm = struct('nn',2);
for ndx = 1:size(F,4)
  figure(1);
  curH = H(:,:,:,ndx);
  subplot(1,2,1); montage2(F(:,:,:,ndx),prm);
  subplot(1,2,2); montage2(curH,prm);
  title(num2str(ndx));
  pause(1);
end

%% test gradient hist.

O = (-4*pi/8-0.0001)*ones(100,100);
M = ones(100,100);
H = gradientHist(single(M),single(O),psize,nbins,1);
figure; montage2(H,prm);

%%


Vx = zeros(100);
Vy = ones(100);

M = sqrt(Vx.^2 + Vy.^2);
O = mod(atan2(Vy,Vx)/2,pi);
H = gradientHist(single(M),single(O),psize,nbins,1);
figure; montage2(H,prm);


%%

V = load('../../data/grab.jab','-mat');

width = 640;
height = 480;
Fs = {}; Hs = {};

ims = [];
for expi = 1:numel(V.x.expDirNames)
  vidfile = fullfile(V.x.expDirNames{expi},V.x.file.moviefilename);
  [readframe,nframes,fid,headerinfo] = get_readframe_fcn(vidfile);
  readFrameFcns.readframe = readframe;
  readFrameFcns.nframes = nframes;
  readFrameFcns.headerinfo = headerinfo;

  [F,H] = genFeatures('ts',10:11,'readframeFcns',readFrameFcns);
  Fs{expi} = F;
  Hs{expi} = H;
  curimmat = rgb2gray(readframe(10));
  
  if size(curimmat,1)>height
    curimmat(height+1:end,:,:) = [];
  elseif size(curimmat,1)<height
    curimmat(end+1:height,:,:) = 0;
  end
  if size(curimmat,2)>width
    curimmat(:,width+1:end,:) = [];
  elseif size(curimmat,2)<width
    curimmat(:,end+1:width,:) = 0;
  end
  
  ims = cat(3,ims,curimmat);
  
  if ~isempty(vidfile) && fid>0,
    fclose(fid);
  end
  
end

figure; montage2(ims);

figure;
for ndx = 1:numel(V.x.expDirNames)
  subplot(1,numel(V.x.expDirNames),ndx);
  vv = Fs{ndx}(:,:,:,1);
  montage2(vv,prm);
end
figure;
for ndx = 1:numel(V.x.expDirNames)
  subplot(1,numel(V.x.expDirNames),ndx);
  montage2(Hs{ndx}(:,:,:,1),prm);
  set(gca,'CLim',[0 0.01]);
end



%%

filtImg = imfilter(single(ims(:,:,1))/255,fspecial('gaussian',2*3*3+1,3));
[M,O] = gradientMag(filtImg); H=gradientHist(M,O,8,nbins,1);
figure; montage2(H,prm);
  set(gca,'CLim',[0 0.01]);

  
  %%
filtImg = imfilter(single(ims(:,:,1))/255,fspecial('gaussian',2*3*3+1,0.5));
H = hog(filtImg,8,15);
figure; V = hogDraw(H); im(V);
  
%%
[readframe,nframes,fid,headerinfo] = get_readframe_fcn('/groups/branson/home/kabram/Desktop/kabram/matlab/adamTracking/videos/vidBig/movie.avi');
T.trx = genTrack(nframes,[640 480]);
save('../../videos/vidBig/trx.mat','-struct','T');
[readframe,nframes,fid,headerinfo] = get_readframe_fcn('/groups/branson/home/kabram/Desktop/kabram/matlab/adamTracking/videos/vid1/movie.avi');
T.trx = genTrack(nframes,[640 480]);
save('../../videos/vid1/trx.mat','-struct','T');
[readframe,nframes,fid,headerinfo] = get_readframe_fcn('/groups/branson/home/kabram/Desktop/kabram/matlab/adamTracking/videos/vid2/movie.avi');
T.trx = genTrack(nframes,[640 480]);
save('../../videos/vid2/trx.mat','-struct','T');

%%
annotateMice('../../videos/vidBig','../../data/grab.jab');
annotateMice('../../videos/vid1','../../data/grab.jab');
annotateMice('../../videos/vid2','../../data/grab.jab');

%%
[readframe,nframes,fid,headerinfo] = get_readframe_fcn('/groups/branson/home/kabram/Desktop/kabram/matlab/adamTracking/videos/vid2/movie.avi');
I1 = single(rgb2gray(readframe(10)))/255;

[M,O] = gradientMag(I1);
[gx,gy] = gradient(I1);

%%
I = I1;
   tic, [M1,O1]=gradientMag(I); toc
   tic, [Gx,Gy]=gradient2(I); M2=sqrt(Gx.^2+Gy.^2);
   O2=mod(atan2(Gy,Gx),pi); toc, mean2(abs(M1-M2))
   d=abs(O1-O2); d(d>pi/2)=pi-d(d>pi/2); mean2(d)

%% Set scorenorm..
h = findall(0,'type','figure','name','JAABA');
vv = guidata(h);
vv.data.windowdata.scoreNorm = 30;

%% for supination. Difference between two flows. Seems to work.
[Vx,Vy,~] = optFlowLk(I0,I1,[],3);

M = sqrt(Vx.^2 + Vy.^2);
sel = M>2;
Z = find(sel(:));
[Y,X] = ind2sub(size(M),Z);
figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);imshow(I0/255);
hold on;
quiver(X,Y,Vx(sel(:)),Vy(sel(:)));
M = sqrt(Vx.^2 + Vy.^2);
O = mod(atan2(Vy,Vx)/2,pi);
O = min(O,pi-1e-6);
H = gradientHist(single(M),single(O),psize,nbins,1);
subplot(1,2,2); montage2(H);

% 
[Vx,Vy,~] = optFlowLk(I0,I1,[],20);

M = sqrt(Vx.^2 + Vy.^2);
sel = M>1;
Z = find(sel(:));
[Y,X] = ind2sub(size(M),Z);
figure('units','normalized','outerposition',[0 0 1 1]);

subplot(1,2,1);imshow(I0/255);
hold on;
quiver(X,Y,Vx(sel(:)),Vy(sel(:)));
M = sqrt(Vx.^2 + Vy.^2);
O = mod(atan2(Vy,Vx)/2,pi);
O = min(O,pi-1e-6);
H = gradientHist(single(M),single(O),psize,nbins,1);
subplot(1,2,2); montage2(H);

% %

[Vx1,Vy1,~] = optFlowLk(I0,I1,[],3);
[Vx2,Vy2,~] = optFlowLk(I0,I1,[],20);

Vxd = Vx1-Vx2;
Vyd = Vy1-Vy2;

M = sqrt(Vxd.^2 + Vyd.^2);
O = mod(atan2(Vyd,Vxd)/2,pi);
O = min(O,pi-1e-6);
H = gradientHist(single(M),single(O),psize,nbins,1);


figure('units','normalized','outerposition',[0 0 1 1]);
subplot(1,2,1);
imshow(I0/255);
% hold on;
% quiver(Vx1-Vx2,Vy1-Vy2);

subplot(1,2,2);
montage2(H);
J = cat(4,I0,I1);
J = repmat(J,[1 1 3 1]);
J = uint8(J);
frame2gif(J,'temp.gif');

%%


M = sqrt(Vx.^2 + Vy.^2);
O = mod(atan2(Vy,Vx)/2,pi);
O = min(O,pi-1e-6);
H = gradientHist(single(M),single(O),psize,nbins,1);
figure; subplot(1,2,1);

imshow(I0/255);
hold on;
quiver(Vx,Vy);
subplot(1,2,2);
montage2(H);

%% Visualize classifier

Fimg = F(:,:,:,1);
Fimg(:) = 0;
Himg = H(:,:,:,1);
Himg(:) = 0;

for ndx = 1:numel(classifier)
  dim = classifier(ndx).dim;
  alpha = classifier(ndx).alpha;
  if dim>numel(Fimg)
    dim = dim-numel(Fimg);
    Himg(dim) = Himg(dim) + abs(alpha);
  else
    Fimg(dim) = Fimg(dim) + abs(alpha);
  end
  
end

figure; 
subplot(1,2,1); montage2(Fimg);
subplot(1,2,2); montage2(Himg);

%%

cl = trainDetectorST('../../data/grab.jab');
save('grabClassifier','cl');
classifyAllDirs('/groups/branson/home/kabram/Desktop/kabram/matlab/adamTracking/videos/20130628','../../data/grab.jab','grabClassifier');
classifyAllDirs('/groups/branson/home/kabram/Desktop/kabram/matlab/adamTracking/videos/20130701','../../data/grab.jab','grabClassifier');


%%

rootdir = '/groups/branson/home/kabram/bransonlab/mayank/adamData/diana/Pre-DTA';
dd = dir(rootdir);

for dndx = 1:numel(dd);
  if strcmp(dd(dndx).name(1),'.'), continue; end
  dd1 = dir(fullfile(rootdir,dd(dndx).name));
  for ndx = 1:numel(dd1)
    if strcmp(dd(ndx).name(1),'.'), continue; end
    curd = fullfile(rootdir,dd(dndx).name,dd1(ndx).name);
    setUpDir(curd,'doforce',true,'features',true,'annotate',false);
    
  end
end
%%
% trx = load('/groups/hantman/hantmanlab/Diana/dta_expts/Post-DTA/10-June/41/10june_postdta_41_1-Basler piA640/trx.mat');
% ips.pf = trx.trx(1).arena.food;
% ips.pm = trx.trx(1).arena.mouth;
% ips.pl = trx.trx(1).arena.perch;

bdir = '/groups/hantman/hantmanlab/Mayank/JayExpts/';
dd = dir(bdir);

for ndx = 1:numel(dd)
  if dd(ndx).name(1) == '.'; continue; end
  if exist(fullfile(bdir,dd(ndx).name,'41'),'dir');
    setUpDir(fullfile(bdir,dd(ndx).name,'41'),'ips',ips);
  end
end


%%
infile = '/groups/branson/home/kabram/bransonlab/mayank/adamData/diana/Pre-DTA/26-May-2013/39/26may_39_predta_3-Basler piA640/trx.mat';
Q = load(infile);
Q.trx(2) = Q.trx(1);
Q.trx(3) = Q.trx(1);
Q.trx(1).x(:) = Q.trx(1).arena.food(1);
Q.trx(1).y(:) = Q.trx(1).arena.food(2);
Q.trx(2).x(:) = Q.trx(1).arena.mouth(1);
Q.trx(2).y(:) = Q.trx(1).arena.mouth(2);
Q.trx(3).x(:) = Q.trx(1).arena.perch(1);
Q.trx(3).y(:) = Q.trx(1).arena.perch(2);

save(infile,'-struct','Q');


%%

bdir = '/groups/branson/home/kabram/bransonlab/mayank/adamData/diana/';

[a,~] = getAllFiles(bdir,'trx.mat');

for ndx = 1:numel(a)
  addAnnotationTracks(a{ndx});
end

[a,~] = getAllFiles(bdir,'features.mat');

for ndx = 1:numel(a)
  cc = regexprep(a{ndx},' ','\\ ');
  cmd = sprintf('touch %s',cc);
  system(cmd);
end

%%

bdir = '/groups/branson/home/kabram/bransonlab/mayank/adamData/diana/';

[~,a] = getAllFiles(bdir,'perframe');

for ndx = 1:numel(a)
  cc = regexprep(a{ndx},' ','\\ ');
  cmd = sprintf('rm -rf %s',cc);
  system(cmd);
end
