function UpdateTargetPosition(targettype,hfly,hfly_extra,pos)
switch targettype,
  case 'fly',
    updatefly(hfly,pos.x,pos.y,pos.theta,pos.a,pos.b);
  case 'larvavani',
    set(hfly,'XData',pos.skeletonx,'YData',pos.skeletony);
    set(hfly_extra,'XData',pos.skeletonx(1),'YData',pos.skeletony(1));
  case 'wingedfly',
    updatewingedfly(hfly,hfly_extra,pos);
  case 'larvacontour',
    updatelarvacontour(hfly,hfly_extra,pos);
  case 'larvasamuel',
    updatelarvasamuel(hfly,hfly_extra,pos);
  case 'larva',
    updatelarvaspecies(hfly,hfly_extra,pos);
end
