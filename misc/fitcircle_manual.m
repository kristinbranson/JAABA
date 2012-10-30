function [xc,yc,radius] = fitcircle_manual(hfig)

if nargin < 1,
  hfig = gcf;
end
hold on;

title('Click points on circle');
[x,y] = getline(hfig,'closed');

% solve for parameters a, b, and c in the least-squares sense by
% using the backslash operator
abc = [x y ones(length(x),1)] \ -(x.^2+y.^2);
a = abc(1); b = abc(2); c = abc(3);

% calculate the location of the center and the radius
xc = -a/2;
yc = -b/2;
radius  =  sqrt((xc^2+yc^2)-c);

% display the calculated center
plot(xc,yc,'mx','LineWidth',2);
text(xc,yc,sprintf('(%.1f, %.1f)',xc,yc),'Color','m','FontWeight','bold');

% plot the entire circle
theta = 0:0.01:2*pi;

% use parametric representation of the circle to obtain coordinates
% of points on the circle
Xfit = radius*cos(theta) + xc;
Yfit = radius*sin(theta) + yc;

plot(Xfit, Yfit,'m-','LineWidth',2);

message = sprintf('The estimated radius is %2.3f pixels', radius);
text(15,15,message,'Color','m','FontWeight','bold','BackgroundColor','k');