function im = VisualizeFlowFeaturesSeq(bdir,fly,fnum,method)
%% inputs.

method = 'hs-sup';

if strcmp(method,'LK')
  fname = 'ff';
elseif strcmp(method,'hs-brightness')
  fname = 'hs_ff';
elseif strcmp(method,'hs-sup')
  fname = 'hs_sup';
else
  error('Unknown method %s',method);
end
fname = [fname 's'];

moviename = fullfile(bdir,'movie.ufmf');
trackfilename = fullfile(bdir,'trx.mat');

%% params

params = getParams;
npatches = params.npatches;
psize = params.psize;
nbins = params.nbins; 
patchsz = params.patchsz;
scale = params.scale;

% this seems to be what the centers correspond to
bincenters = linspace(0,pi,nbins+1);
bincenters = bincenters(1:nbins);

% mayank divides by 2
bincenters2 = bincenters*2;
dt = mean(diff(bincenters2));
binedges2 = [bincenters2(1)-dt/2,(bincenters2(1:end-1)+bincenters2(2:end))/2,bincenters2(end)+dt/2];

%% Sequence

seqlen = 4;
nshow = 10;
fstart = fnum-seqlen;
fend = fnum+seqlen;

%% 

tracks = load(trackfilename);
tracks = tracks.trx;

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);

im1 = [];im2 = [];
count = 1;
for fno = fstart:fend+1
  i1 = readfcn(fno);
  trackndx = fno - tracks(fly).firstframe + 1;
  locy = round(tracks(fly).y(trackndx));
  locx = round(tracks(fly).x(trackndx));
  im1(:,:,count) = extractPatch(i1,...
    locy,locx,tracks(fly).theta(trackndx),patchsz);
  count = count+1;
end

F = zeros(npatches,npatches,nbins,2*seqlen+1);
parfor yy = 1:npatches
  for xx = 1:npatches
    for oo = 1:nbins
      pfname = fullfile(bdir,'perframe',sprintf('%s_%02d_%02d_%d.mat',fname,yy,xx,oo));
      q = load(pfname);
      trackndx = fnum - tracks(fly).firstframe + 1;
      F(yy,xx,oo,:) = q.data{fly}(trackndx-seqlen:trackndx+seqlen);
      
    end
  end
end

Fall = F;
F = mean(Fall,4);

%% plot
icol = hsv(nshow);
idx = round(linspace(1,2*seqlen+2,nshow));
cimg = zeros(size(im1,1),size(im1,2),3);
count = 1;
for ndx = idx(:)'
  for ch = 1:3
    cimg(:,:,ch) = cimg(:,:,ch)+im1(:,:,ndx)*icol(count,ch);
  end
  count = count+1;
end
for ch = 1:3
  cimg(:,:,ch)= cimg(:,:,ch)/sum(icol(:,ch));
end

cimg = min(im1,[],3);
% cimg = zeros(size(im1,1),size(im1,2),3);
% tim1 = im1(:,:,idx);
% [mimg,mndx] = min(tim1,[],3);
% for ndx = 1:numel(idx)
%   for ch=1:3
%     jj = cimg(:,:,ch);
%     kk = tim1(:,:,ndx);
%     selpx = abs(kk-mimg)<65;
%     jj(selpx)= jj(selpx) + kk(selpx)*icol(ndx,ch);
%     cimg(:,:,ch) = jj;
%   end
% end
% for ch = 1:3
%   cimg(:,:,ch)= cimg(:,:,ch)/sum(icol(:,ch));
% end

hfig = figure;
clf;
hax = axes('Position',[0,0,1,1]);
set(hfig,'Units','pixels','Position',get(0,'ScreenSize'));


him = imshow(uint8(imresize(cimg,scale)));
axis image;
truesize;
hold on;
axis off;
%%
colors = hsv(nbins);

[nr,nc,~] = size(im1);
% maxv2 = max(F(:));
maxv2 = 1;

h = [];
for xi = 1:ceil(nc/psize),
  cx = (psize/2 +  (xi-1)*psize)*scale1 +1;
  if cx+psize/2 > nc*scale,
    break;
  end
  for yi = 1:ceil(nr/psize),
    cy = (psize/2  + (yi-1)*psize)*scale+1;
    if cy+psize/2 > nr*scale,
      break;
    end
    
    for bini = 1:nbins,
      tmp = linspace(binedges2(bini),binedges2(bini+1),20);
      xcurr = cx + [0,psize/2*cos(tmp),0]*scale;
      ycurr = cy + [0,psize/2*sin(tmp),0]*scale;
      h(yi,xi,bini) = patch(xcurr,ycurr,colors(bini,:),'LineStyle','none','FaceAlpha',min(1,F(yi,xi,bini)/maxv2));
    end
    
  end
end
im = getframe(hax);
