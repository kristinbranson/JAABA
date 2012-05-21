function updatelarvaspecies(hfly,hfly_extra,pos)

set(hfly,'XData',pos.xspine,'YData',pos.yspine);
%set(hfly,'XData',[pos.xcontour;nan;pos.xspine],'YData',[pos.ycontour;nan;pos.yspine]);
xhead = pos.x + 2*pos.a*cos(pos.theta);
yhead = pos.y + 2*pos.a*sin(pos.theta);
set(hfly_extra,'XData',xhead,'YData',yhead);