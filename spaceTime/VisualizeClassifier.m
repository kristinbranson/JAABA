function im = VisualizeClassifier(jabfile)


%% params.
params = getParams;
npatches = params.npatches;
psize = params.psize;
nbins = params.nbins; 
patchsz = params.patchsz;
scale = params.scale;

wd = params.wd;

%% load an example fly.

Q = load(jabfile,'-mat');

bdir = Q.x.expDirNames{1};
moviename = fullfile(bdir,'movie.ufmf');
trackfilename = fullfile(bdir,'trx.mat');
tracks = load(trackfilename);
tracks = tracks.trx;

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
trackndx = 1; fly = 1;
fnum = tracks(fly).firstframe;
im1 = readfcn(fnum);

locy = round(tracks(fly).y(trackndx));
locx = round(tracks(fly).x(trackndx));
im1 = extractPatch(im1,...
  locy,locx,tracks(fly).theta(trackndx),patchsz);

im1 =[im1 im1];

%% Find Hog/HOF components.

ncl = numel(Q.x.classifierStuff);

for clnum = 1:ncl
  
  H{clnum} = zeros(npatches,npatches,nbins);
  F{clnum} = H{clnum};
  
  cl = Q.x.classifierStuff(clnum).params;
  fnames = Q.x.classifierStuff(clnum).featureNames;
  for ndx = 1:numel(cl)
    curf = fnames{cl(ndx).dim}{1};
    val = cl(ndx).alpha;
    prts = strsplit(curf,'_');
    loc = str2double(prts(end-2:end));
    if strcmp(curf(1:2),'hf')
      H{clnum}(loc(1),loc(2),loc(3))=H{clnum}(loc(1),loc(2),loc(3))+val;
    else
      F{clnum}(loc(1),loc(2),loc(3))=F{clnum}(loc(1),loc(2),loc(3))+val;
    end
  end
end
%% Draw them

% For HOG
bincenters = linspace(0,pi,nbins+1);
bincenters = bincenters(1:nbins);
dt = mean(diff(bincenters));
binedges = [bincenters(1)-dt/2,(bincenters(1:end-1)+bincenters(2:end))/2,bincenters(end)+dt/2];

% For HOF
bincenters2 = bincenters*2;
dt = mean(diff(bincenters2));
binedges2 = [bincenters2(1)-dt/2,(bincenters2(1:end-1)+bincenters2(2:end))/2,bincenters2(end)+dt/2];

maxv = 3;
im = {};
for clnum = 1:ncl
  hfig = figure;
  clf;
  hax = axes('Position',[0,0,1,1]);
  set(hfig,'Units','pixels','Position',get(0,'ScreenSize'));
  
  im1curr = im1;
  
  him = imshow(imresize(im1curr,scale));
  axis image;
  truesize;
  colormap gray;
  hold on;
  axis off;
  
  colors = hsv(nbins);
  
  [nr,nc,~] = size(im1);
  nc = nc/2;
  
  h = [];
  for xi = 1:ceil(nc/psize),
    cx = (psize/2 + (xi-1)*psize)*scale + 1 ;
    if cx+psize/2 > nc*scale,
      break;
    end
    for yi = 1:ceil(nr/psize),
      cy = (psize/2 + (yi-1)*psize)*scale + 1;
      if cy+psize/2 > nr*scale,
        break;
      end
      
      for bini = 1:nbins,
        tmp = linspace(binedges2(bini),binedges2(bini+1),20);
        xcurr = cx + [0,psize/2*cos(tmp),0]*scale;
        ycurr = cy + [0,psize/2*sin(tmp),0]*scale;
        h(yi,xi,bini) = patch(xcurr,ycurr,colors(bini,:),'LineStyle','none','FaceAlpha',min(1,F{clnum}(yi,xi,bini)/maxv));
      end
      
    end
  end
  
  colors = hsv(nbins);
  colors = colors([ (end/2+1):end 1:end/2],:);
  hogpatch = [wd wd -wd -wd wd;-psize psize psize -psize -psize]/2;
  h = [];
  for xi = 1:ceil(nc/psize),
    cx = (psize/2 + 1 + (xi-1)*psize)*scale;
    if cx+psize/2 > nc*scale,
      break;
    end
    for yi = 1:ceil(nr/psize),
      cy = (psize/2 + 1 + (yi-1)*psize)*scale;
      if cy+psize/2 > nr*scale,
        break;
      end
      
      for bini = 1:nbins,
        tmp = bincenters(bini);
        curpatch = [cos(tmp) -sin(tmp); sin(tmp) cos(tmp)]*hogpatch;
        xcurr = nc*scale + cx + curpatch(1,:)*scale;
        ycurr = cy + curpatch(2,:)*scale;
        h(yi,xi,bini) = patch(xcurr,ycurr,colors(bini,:),'LineStyle','none','FaceAlpha',min(1,H{clnum}(yi,xi,bini)/maxv));
      end
      
    end
  end
  im{clnum} = getframe(hax);
end
