function imnorm = compute_spacetime_transform(im,x,y,theta,a,b,meana,meanb)

clipsize = 5;
boxwidth2 = round(meanb*8);
boxheight2 = round(meana*4);

if isnan(x) || isnan(y) || isnan(theta) || isnan(a) || isnan(b)
  imnorm=nan(boxheight2*2+1,boxwidth2*2+1);
  return
end

x1 = floor(max(1,x-a*clipsize));
x2 = floor(min(size(im,2),x+a*clipsize));
y1 = ceil(max(1,y-a*clipsize));
y2 = ceil(min(size(im,1),y+a*clipsize));
%[x1,x2,y1,y2] = ellipse_to_bounding_box(pos.x,pos.y,pos.a*2,pos.b*2,pos.theta);
% y1 = floor(y1);
% y2 = ceil(y2);
% x1 = floor(x1);
% x2 = ceil(x2);
imbb = im(y1:y2,x1:x2);

[nrcurr,nccurr,~] = size(imbb);
theta = theta+pi/2;
R = [cos(theta),-sin(theta),0
  sin(theta),cos(theta),0
  0,0,1];
scalex = meanb/b;
scaley = meana/a;
S = [scalex,0,0;0,scaley,0;0,0,1];
tform = maketform('affine',double(R*S));
%tform = maketform('affine',R);
imnorm = imtransform(imbb,tform,...
  'UData',double([x1-x,nccurr-x+x1]),...
  'VData',double([y1-y,nrcurr-y+y1]),...
  'XData',double([-boxwidth2,boxwidth2]),...
  'YData',double([-boxheight2,boxheight2]),...
  'FillValues',0,...
  'XYScale',1);

%   subplot(1,4,3);
%   imagesc(imnorm,[0,1]);
%   axis image;