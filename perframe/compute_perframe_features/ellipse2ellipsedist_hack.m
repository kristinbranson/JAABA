% AR 3/8/2018
function [d,x1i,y1i,theta1i,x2j,y2j,theta2j] = ellipse2ellipsedist_hack(xc1,yc1,a1,b1,theta1,xc2,yc2,a2,b2,theta2,npoints)

if nargin < 11
    npoints = 100;
end

[x1,y1,ntheta1] = ellipsepoints(a1,b1,xc1,yc1,theta1,npoints);
[x2,y2,ntheta2] = ellipsepoints(a2,b2,xc2,yc2,theta2,npoints);


distances = dist2([x1(:),y1(:)],[x2(:),y2(:)]);

% distances(:) turns it into 1 colunm vector, then can take min
[d,idx] = min(distances(:));
[i,j] = ind2sub(size(distances),idx);

d = sqrt(d);
x1i = x1(i);
x2j = x2(j);
y1i = y1(i);
y2j = y2(j);
theta1i = ntheta1(i);
theta2j = ntheta2(j);

d = reshape(d,size(xc1));

end

