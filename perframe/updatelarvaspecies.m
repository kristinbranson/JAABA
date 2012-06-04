function updatelarvaspecies(hfly,hfly_extra,pos)

set(hfly,'XData',pos.xspine([1,6,11]),'YData',pos.yspine([1,6,11]));
%set(hfly,'XData',[pos.xcontour;nan;pos.xspine],'YData',[pos.ycontour;nan;pos.yspine]);
% xhead = pos.x + 2*pos.a*cos(pos.theta);
% yhead = pos.y + 2*pos.a*sin(pos.theta);
xhead = pos.xspine(1);
yhead = pos.yspine(1);
set(hfly_extra,'XData',xhead,'YData',yhead);