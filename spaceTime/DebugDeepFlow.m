params = getParams;
npatches = params.npatches;
psize = params.psize;
nbins = params.nbins; 
patchsz = params.patchsz;

optflowwinsig = params.optflowwinsig ;
optflowsig = params.optflowsig ;
optreliability = params.optreliability ;


maxflow = 5;
fly_thres = 90;
flow_thres = sqrt(2)/4;
scale = 6;

makeVideo = false;

% expdir = '../mated10_20140714T131113';
% frames = 10900 +(1:200);
% fly = 1;
% fly = 1;
% frames = 11043; % large middle leg movement
% frames = 10925; % large middle leg movement
%  frames = 10914; 

expdir = '/home/mayank/Work/FlySpaceTime/walkMovies/SS03500';
fly = 1;
% frames = 10200 + (1:100);
% frames = 735:780;
frames = 744;
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

%%

im1curr = im1(:,:,1);
im2curr = im2(:,:,1);
imwrite(repmat(uint8(im2curr),[1 1 3]),'epicflow/deepmatching/im2.png');
imwrite(repmat(uint8(im1curr),[1 1 3]),'epicflow/deepmatching/im1.png');
cmd = sprintf('%s %s %s','epicflow/deepmatching/deepmatching ',...
  ' epicflow/deepmatching/im1.png epicflow/deepmatching/im2.png -resize 320 320 -ngh_rad 120 ',...
  ' -out epicflow/deepmatching/temp.out');
[a,b] = system(cmd);
if ~a, fprintf(b); end
%%

ofile = 'epicflow/deepmatching/temp.out';
sz = size(im1curr);
A = dlmread(ofile);
Vx = zeros(sz); 
Vy = Vx; Vs = Vx;
idx = sub2ind(sz,round(A(:,2)-0.0006)+1,round(A(:,1)-0.0006)+1);
Vx(idx) = A(:,3)-A(:,1);
Vy(idx) = A(:,4)-A(:,2);
Vs(idx) = A(:,5);
uv = cat(3,Vx,Vy);
figure;  ax = [];
ns = 5;
ax(1) = subplot(1,ns,1);
imshowpair(im1curr,im2curr);
% hold on;
% quiver(Vx,Vy,0);

[xx,yy] = meshgrid((-5:5),(-5:5));

uv = uv(2:2:end,2:2:end,:);
Vss = Vs(2:2:end,2:2:end);
selpx = Vss<3.5;
% selpx = Vs<4;
tic;
im1sm = imresize(im1curr,0.5);
im2sm = imresize(im2curr,0.5);

uvold = estimate_flow_interface(im1sm,im2sm,'hs-brightness',...
  {'max_warping_iter',2});
toc;
% uvoldss = uvold(2:2:end,2:2:end,:);
uv(selpx)=uvold(selpx);
uv = imresize(uv,2);
uv(1:11,1:11,:)=cat(3,xx,yy);
ax(2) = subplot(1,ns,2);
imshow(flowToColor(uv));
ax(3) = subplot(1,ns,3);
imagesc(Vs); axis image;
ax(4) = subplot(1,ns,4);
imshow(imresize(flowToColor(uvold),2)); 
ax(5) = subplot(1,ns,5);
imshow(flowToColor(cat(3,Vx,Vy))); 

linkaxes(ax);

