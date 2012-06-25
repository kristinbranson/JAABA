function updatelarvacontour(hfly,hfly_extra,pos)

set(hfly,'XData',pos.xcontour,'YData',pos.ycontour);
if all(isfield(pos,{'xspine','yspine'})),
  set(hfly_extra,'XData',pos.xspine,'YData',pos.yspine);
end