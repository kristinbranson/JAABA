%% Read the images.


% J = load('/home/mayank/Work/FlyBowl/newfeatures_walk_2.jab','-mat');
% 
% count = zeros(1,2);
% for ndxx = 1:numel(J.x.labels.t0s)
%   for bout = 1:numel(J.x.labels.t0s{ndxx})
%     if J.x.labels.t1s{ndxx}(bout)>1000, continue;end
%     fly = J.x.labels.flies(ndxx);
%     frames = J.x.labels.t0s{ndxx}(bout):J.x.labels.t1s{ndxx}(bout)-1;
%     bothim = {};
%     for tt = 1:2
%       if tt ==1
%         expdir = '/home/mayank/Work/FlySpaceTime/walkMovies/SS03500_test';
%         ftrname = 'DSs';
%       else
%         expdir = '/home/mayank/Work/FlySpaceTime/walkMovies/SS03500';
%         ftrname = 'hs_sups';
%       end

%%
clear tt;
expdir = '/home/mayank/Work/FlySpaceTime/walkMovies/SS03500_test';
ftrname = 'DSs';
% fly = 2;
% frames = 51:200;
fly = 4;
frames = 3342:3353;
% frames = 3532:3543;
% selftrs = [6 7 1; 2 4 1; 6 3 1];
selftrs = [3 7 1; 6 5 1; 3 3 1];
savefigure = false;

%% params.
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
flow_thres = sqrt(2)/4;


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


%% Read the features;

ftrVal = [];
for ndx = 1:size(selftrs,1)
  curname = sprintf('%s_%02d_%02d_%d.mat',ftrname,...
    selftrs(ndx,1),selftrs(ndx,2),selftrs(ndx,3));
  Q = load(fullfile(expdir,'perframe',curname));
  startndx = frames(1)-tracks(fly).firstframe+1;
  endndx  = frames(end)-tracks(fly).firstframe+1;
  ftrVal(ndx,:) = Q.data{fly}(startndx:endndx);
end
maxval = 2;
ftrVal(ftrVal>maxval) = maxval;
%% Visualize

minImg = min(im1,[],3);

nftrs = size(selftrs,1);
padsz = 50; margin = 4;
tsz = (padsz-2*margin)/nftrs;
padImg = padgrab(minImg,255,1,nc+padsz,1,nr,1,size(minImg,3));

if exist('tt','var')
  ff = figure(104+tt);
else
  ff = figure;
end
imshow(uint8(imresize(padImg,scale)));hold on;

xidx = linspace(margin,nc-margin,numel(frames));
cc = hsv(nftrs);
for ndx = 1:nftrs
  ybase = nc + (padsz-margin)/nftrs*ndx;
  
  plot(xidx*scale,(ybase-ftrVal(ndx,:)*tsz/maxval)*scale,'Color',cc(ndx,:));
  px = [0 psize psize 0 0] + (selftrs(ndx,2)-1)*psize;
  py = [0 0 psize psize 0] + (selftrs(ndx,1)-1)*psize;
  patch(px*scale,py*scale,cc(ndx,:),'FaceColor',cc(ndx,:),'FaceAlpha',0.075,'LineStyle','none');
end


for ndx = 1:npatches-1
  plot([0.5 nr*scale+0.5],[ndx*psize ndx*psize]*scale,'y');
  plot([ndx*psize ndx*psize]*scale,[0.5 nr*scale+0.5],'y');
  
end

% fd = getframe(ff);
% bothim{tt} = fd.cdata;
%%
%     end
%     
%     if savefigure,
%       outimgname = sprintf('%s_%s_Fly%d_From%dTo%d_%s.png',...
%         J.x.labels.names{ndxx}{bout},...
%         expname,fly,frames(1),frames(end),...
%         datestr(now,'yyyymmdd'));
%       outname = fullfile('..','figures',outimgname);
%       imwrite(cat(2,bothim{1},bothim{2}),outname);
%     end
% 
%   end
% end