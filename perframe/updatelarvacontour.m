function updatelarvacontour(hfly,hfly_extra,pos)

if isempty(pos.xcontour),
  set(hfly,'XData',[],'YData',[]);
else
  set(hfly,'XData',[pos.xcontour(:);pos.xcontour(1)],'YData',[pos.ycontour(:);pos.ycontour(1)]);
end
if all(isfield(pos,{'xspine','yspine'})),
  set(hfly_extra(1),'XData',pos.xspine,'YData',pos.yspine);
  set(hfly_extra(2),'XData',pos.xspine(1),'YData',pos.yspine(1));
end