function ConvertCtraxMouseToTrx(indir,intrxfile,moviefilename,outtrxfile,arenaSize,framerate,method,sex,frames2sel)
% function ConvertCtraxMouseToTrx(indir,arenaSize,dt,sex)
% Convert tracking results from ctrx into trx format. Mainly for vivek's
% movies.
% indir: experiment directory
% intrxfile: input trx file (ctrax output)
% moviefilename:
% outtrxfile: output trx file
% arenaSize: The size of the edge of the box in mm
% framerate
% method: 'manual','auto', 'both'. Both will try out auto and if it fails
% will ask for for manual intervention.
% sex: animals sex - optional 
% frames2sel: Frames to use for corner detection.

% intrxfile = 'trx.mat';
% moviefilename = 'movie.mov';
% outtrxfile = 'registered_trx.mat';

if nargin< 8
  sex = 'm';
end

sampleframes = 1;

trx = load_tracks(fullfile(indir, intrxfile));

[~,~,ext] = fileparts(moviefilename);
if strcmp(ext,'.ufmf')
  isufmf = true;
else
  isufmf = false;
end

[readframe,nframes,movie_fid,movieheaderinfo] = ...
  get_readframe_fcn(fullfile(indir,moviefilename));

imgfile = fullfile(indir,'corner_detection.png');

if nargin<9
  frames2read = randsample(nframes,1);
else
  frames2read=randsample(frames2sel,1);
end
tl = []; tr = []; bl = []; br = [];
img = readframe(frames2read);

if strcmp(method,'auto'),
  
  [success,C] = cornerDetection(img,imgfile);
  if success,
    ctl = C(1,:);
    ctr = C(2,:);
    cbl = C(3,:);
    cbr = C(4,:);
  else
    fprintf('Could not detect corners automatically for %s. Please label them manually\n',indir);
    return;
  end
elseif strcmp(method,'manual'),
  [ctl ctr cbl cbr] = getCornersManual(img);
elseif strcmp(method,'both');
  [success,C] = cornerDetection(img,imgfile);
  if success,
    ctl = C(1,:);
    ctr = C(2,:);
    cbl = C(3,:);
    cbr = C(4,:);
  else
    [ctl ctr cbl cbr] = getCornersManual(img);
  end
else
  fprintf('Invalid method: It should be either ''manual'', ''auto'' or ''both''\n');
  return;
end


tl(end+1,:) = ctl;
tr(end+1,:) = ctr;
bl(end+1,:) = cbl;
br(end+1,:) = cbr;

ftl = median(tl,1);
ftr = median(tr,1);
fbl = median(bl,1);
fbr = median(br,1);

image_size = size(img);

d_top_px = sqrt( sum((ftl-ftr).^2));
d_left_px = sqrt( sum((ftl-fbl).^2));
d_bottom_px = sqrt( sum((fbl-fbr).^2));
d_right_px = sqrt( sum((ftr-fbr).^2));
mean_size = mean([d_top_px d_right_px d_left_px d_bottom_px]);
scaleFactor = arenaSize/mean_size;


for ndx = 1:numel(trx)
  trx(ndx).x_mm = trx(ndx).x*scaleFactor;
  if isufmf,
    trx(ndx).y_mm = trx(ndx).y*scaleFactor;
  else
    trx(ndx).y_mm = (image_size(1)-trx(ndx).y)*scaleFactor;
    trx(ndx).y = (image_size(1)-trx(ndx).y);
  end
  trx(ndx).a_mm = trx(ndx).a*scaleFactor;
  trx(ndx).b_mm = trx(ndx).b*scaleFactor;
  if isufmf,
    trx(ndx).theta_mm = trx(ndx).theta;
  else
    trx(ndx).theta_mm = 2*pi - trx(ndx).theta;
    trx(ndx).theta = 2*pi - trx(ndx).theta;
  end
  trx(ndx).dt = repmat(1/framerate,1,numel(trx(ndx).x)-1);
  trx(ndx).sex = repmat({sex},1,numel(trx(ndx).x));
  trx(ndx).arena.tl = ftl*scaleFactor;
  trx(ndx).arena.tr = ftr*scaleFactor;
  trx(ndx).arena.bl = fbl*scaleFactor;
  trx(ndx).arena.br = fbr*scaleFactor;
end

save(fullfile(indir,outtrxfile),'trx');
if movie_fid > 0, fclose(movie_fid); end

function [tl tr bl br] = getCorners(img)
  
pp = prctile(double(img(:)),[10:3:90]);
[~,pndx ]= max(pp(2:end)-pp(1:end-1));
tr = pp(pndx-3);

bb = img>tr;
bb = imfill(bb,'holes');
[y,x] = ind2sub(size(bb),find(bb));

ll = bwlabel(bb);
props = regionprops(ll,'Area');
[~,sel] = max([props.Area]);
bb = ll == sel;
cc = imclose(bb,strel('disk',5));
dd = imopen(cc,strel('disk',5));
corners = corner(dd,4,'FilterCoefficients',fspecial('gaussian',[25 1],5));

[~,xord] = sort(corners(:,1));
lefts = xord(1:2); rights = xord(3:4);
[~,yord] = sort(corners(:,2));
tops = yord(1:2); bottoms = yord(3:4);
tl = corners(intersect(lefts,tops),:);
tr = corners(intersect(rights,tops),:);
bl = corners(intersect(lefts,bottoms),:);
br = corners(intersect(rights,bottoms),:);

if size(tl,1)>1 || size(tr,1)>1 || size(bl,1)>1 || size(br,1)>1,
  tl = []; tr = []; bl = []; br = [];
end

function [tl tr bl br] = getCornersManual(img)
  
f = figure; h = imshow(img); 
title('Click on the top left corner, double click to finalize');
p1 = impoint(gca);
wait(p1);

title('Click on the top right corner, double click to finalize');
p2 = impoint(gca);
wait(p2);

title('Click on the bottom right corner, double click to finalize');
p3 = impoint(gca);
wait(p3);

title('Click on the bottom left corner, double click to finalize');
p4 = impoint(gca);
wait(p4);

tl = p1.getPosition;
tr = p2.getPosition;
br = p3.getPosition;
bl = p4.getPosition;
close(f);
pause(0.1);
