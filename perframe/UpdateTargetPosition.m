function UpdateTargetPosition(targettype,hfly,pos)
switch targettype,
  case 'fly',
    updatefly(hfly,pos.x,pos.y,pos.theta,pos.a,pos.b);
  case 'larva',
    set(hfly,'XData',pos.skeletonx,'YData',pos.skeletony);
end
