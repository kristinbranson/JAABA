%% Corner Detection for Vivek
function [success, C] = cornerDetection(moviefile,frameRange,showResult)
% function C = cornerDetection(moviefile,frameRange,showResult)
% Detects corners 
% frameRange is the range of the frames to be used to detect corners.
% Show results shows the result of corner detection.

count = 0;
success = false;

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(...
  moviefile);
while count < 10,
%% Read a frame.
count = count+1;

if isempty(frameRange),
  frameRange = 1:nframes;
end

rfr = randsample(frameRange,1);
I = readfcn(rfr);


%%
CM = cornermetric(I(:,:,1));
bw = CM>0.000005;
bw([1:10 end-10:end],:) = 0;
bw(:,[1:10 end-10:end]) = 0;
bwl = bwlabel(bw);
bwcenter = regionprops(bwl,'Centroid','Area');
small2rem = find([bwcenter.Area]<3);

xc = []; yc = [];
for ndx = 1:numel(bwcenter)
xc(ndx) = bwcenter(ndx).Centroid(1);
yc(ndx) = bwcenter(ndx).Centroid(2);
end
xc(small2rem) = size(I,2)/2;
yc(small2rem) = size(I,1)/2;
minx = min(xc);
miny = min(yc);
maxx = max(xc);
maxy = max(yc);

tol = 10;
topleft = find( (abs(xc-minx)<tol) & (abs(yc-miny)<tol));
topright = find( (abs(xc-maxx)<tol) & (abs(yc-miny)<tol));
bottomleft = find( (abs(xc-minx)<tol) & (abs(yc-maxy)<tol));
bottomright = find( (abs(xc-maxx)<tol) & (abs(yc-maxy)<tol));

% C(1,:) = bwcenter(topleft).Centroid; 
% C(2,:) = bwcenter(topright).Centroid; 
% C(3,:) = bwcenter(bottomleft).Centroid; 
% C(4,:) = bwcenter(bottomright).Centroid; 

if isempty(topleft) || isempty(topright) || isempty(bottomleft) || isempty(bottomright),
continue;
end
  
if numel(topleft)>1 || numel(topright) > 1 || numel(bottomleft)>1 || numel(bottomright)>1,
  topleft = topleft(argmax([bwcenter(topleft).Area]));
  topright = topright(argmax([bwcenter(topright).Area]));
  bottomleft = bottomleft(argmax([bwcenter(bottomleft).Area]));
  bottomright = bottomright(argmax([bwcenter(bottomright).Area]));
end
  
success = true; break;
end

if fid> 0,
  fclose(fid);
end

if ~success,
  C = [];
  return;
end

% selImg = (bwl == topleft) | (bwl == bottomleft) | (bwl == topright) | (bwl == bottomright);
% 
% 
% CM(~selImg) = 0;
% corner_peaks = imregionalmax(CM);
% corner_idx = find(corner_peaks == true);
% [xx yy] = ind2sub(size(I(:,:,1)),corner_idx);
% C = [];
C(:,1) = yc([topleft topright bottomleft bottomright]);
C(:,2) = xc([topleft topright bottomleft bottomright]);

if showResult,
  figure;
  imshow(I);
  hold on
  plot(C(:,2), C(:,1), 'r*');
end

%{

%% Detect corners

C = corner(I(:,:,1),'Harris',4,'SensitivityFactor',0.001);

%% Show Corners
imshow(I);
hold on
plot(C(:,1), C(:,2), 'r*');

%% Corner metric

hx = [];
figure;
hx(1) = subplot(1,3,1);
imshow(I);
title('Original Image');
CM = cornermetric(I(:,:,1));

CM_adjusted = imadjust(CM);
hx(2) = subplot(1,3,2);
imshow(CM_adjusted);
title('Corner Metric');


corner_peaks = imregionalmax(CM);
corner_idx = find(corner_peaks == true);
[r g b] = deal(I(:,:,1));
r(corner_idx) = 255;
g(corner_idx) = 255;
b(corner_idx) = 0;
RGB = cat(3,r,g,b);
hx(3) = subplot(1,3,3);
imshow(RGB);
title('Corner Points');
linkaxes(hx);
%}