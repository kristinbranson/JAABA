stationary = true;
patchsz = 80;

moviename = '/home/mayank/Work/FlySpaceTime/mated10_20140714T131113/movie.ufmf';
trxfilename = '/home/mayank/Work/FlySpaceTime/mated10_20140714T131113/trx.mat';
% fnum = 14003; fly = 2;
fly = 2;

fall = figure;
%%
mm = 1;
count= 1;
tott = 0;
allI = [];
allF{mm} = [];
for fnum = 39017:1:39026 %39007:10:39106
[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
tracks = load(trxfilename);
tracks = tracks.trx;
im1 = readfcn(fnum);
im2 = readfcn(fnum+1);

trackndx = fnum - tracks(fly).firstframe + 1;
locy = round(tracks(fly).y(trackndx));
locx = round(tracks(fly).x(trackndx));
im1 = extractPatch(im1,...
  locy,locx,tracks(fly).theta(trackndx),patchsz);
if stationary
  locy = round(tracks(fly).y(trackndx+1));
  locx = round(tracks(fly).x(trackndx+1));
end
im2 = extractPatch(im2,...
  locy,locx,tracks(fly).theta(trackndx),patchsz);

 
im1 = repmat(im1,[1 1 3]);
im2 = repmat(im2,[1 1 3]);


% uv = estimate_flow_interface(im1,im2,'classic-c-a'); 
tic;
[p,q,~] = size(im1);
uv = LDOF(im1,im2,get_para_flow(p,q),0);
tott = tott + toc;
imflow = flowToColor(uv); 
f = figure;
hh = imshowpair(im1,im2);
allI = [allI get(hh,'CData')];
close(f);
allF{mm} = [allF{mm} imflow];
count= count+1;
end
figure(fall);
subplot(1,2,1); imshow(allI);
subplot(1,2,2); imshow(allF{mm});
title(sprintf('%s,%.3f',methods{mm},tott/10));
