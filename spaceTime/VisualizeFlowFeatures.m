function im = VisualizeFlowFeatures(bdir,fly,fnum,stationary,varargin)
%% inputs.

[method,moviename,trxfilename] = myparse(varargin,...
  'method','deep-sup',...
  'moviename','movie.ufmf','trxfilename','trx.mat');

if strcmp(method,'LK')
  fname = 'ff';
elseif strcmp(method,'hs-brightness')
  fname = 'hs_ff';
elseif strcmp(method,'hs-sup')
  fname = 'hs_sup';
elseif strcmp(method,'deep-sup')
  fname = 'DS';
else
  error('Unknown method %s',method);
end
if stationary,
  fname = [fname 's']; 
end

% fnum = 10125;
% fly = 1;
% fstart = 10000;
% stationary = true;
% 
% bdir = '../mated10_20140714T131113/';
moviename = fullfile(bdir,moviename);
trackfilename = fullfile(bdir,trxfilename);

%% params
params = getParams;
npatches = params.npatches;
psize = params.psize;
nbins = params.nbins; 
patchsz = params.patchsz;
scale = params.scale;

%% compute the bins

% this seems to be what the centers correspond to
bincenters = linspace(0,pi,nbins+1);
bincenters = bincenters(1:nbins);

bincenters2 = bincenters*2;
dt = mean(diff(bincenters2));
binedges2 = [bincenters2(1)-dt/2,(bincenters2(1:end-1)+bincenters2(2:end))/2,bincenters2(end)+dt/2];

%% 

tracks = load(trackfilename);
tracks = tracks.trx;

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
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

if ~strcmp(method,'hs_sup')
im2 = extractPatch(im2,...
  locy,locx,tracks(fly).theta(trackndx),patchsz);
else
im2 = extractPatch(im2,...
  locy,locx,tracks(fly).theta(trackndx+1),patchsz);  
end

F = zeros(npatches,npatches,nbins);
for yy = 1:npatches
  for xx = 1:npatches
    for oo = 1:nbins
      pfname = fullfile(bdir,'perframe',sprintf('%s_%02d_%02d_%d.mat',fname,yy,xx,oo));
      q = load(pfname);
      trackndx = fnum - tracks(fly).firstframe + 1;
      F(yy,xx,oo) = q.data{fly}(trackndx);
      
    end
  end
end
  
% plot

hfig = figure('Visible','off');
clf;
hax = axes('Position',[0,0,1,1]);
set(hfig,'Units','pixels','Position',get(0,'ScreenSize'));

im1curr = im1;
im2curr = im2;


him = imshowpair(imresize(im1curr,scale),imresize(im2curr,scale),'Scaling','none');
axis image;
truesize;
colormap gray;
hold on;
axis off;

colors = hsv(nbins);

[nr,nc,~] = size(im1);
% maxv2 = max(F(:));
maxv2 = 2;

h = [];
for xi = 1:ceil(nc/psize),
  cx = (psize/2  + (xi-1)*psize)*scale+ 1;
  if cx+psize/2 > nc*scale,
    break;
  end
  for yi = 1:ceil(nr/psize),
    cy = (psize/2 + (yi-1)*psize)*scale+ 1;
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
truesize(hfig);
im = getframe(hax);
pause(0.01);
close(hfig);
