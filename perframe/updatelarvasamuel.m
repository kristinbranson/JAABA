function updatelarvasamuel(hfly,hfly_extra,pos)

set(hfly,'XData',[pos.xcontour,nan,pos.xspine'],'YData',[pos.ycontour,nan,pos.yspine']);
set(hfly_extra,'XData',[pos.xhead,pos.xmid],'YData',[pos.yhead,pos.ymid]);