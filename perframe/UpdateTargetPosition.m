function UpdateTargetPosition(targettype,hfly,hfly_extra,pos)
switch targettype,
  case 'fly',
    updatefly(hfly,pos.x,pos.y,pos.theta,pos.a,pos.b);
  case 'larvavani',
    set(hfly,'XData',pos.skeletonx,'YData',pos.skeletony);
    set(hfly_extra,'XData',pos.skeletonx(1),'YData',pos.skeletony(1));
  case 'wingedfly',
		updatewingedfly(hfly,hfly_extra,pos);
	case 'center_and_orientation',
		xdata = pos.x + [0,cos(pos.theta)*pos.a];
		ydata = pos.y + [0,sin(pos.theta)*pos.a];
		set(hfly,'XData',xdata,'YData',ydata);
  case 'larvacontour',
    updatelarvacontour(hfly,hfly_extra,pos);
  case 'larvasamuel',
    updatelarvasamuel(hfly,hfly_extra,pos);
  case 'larva',
    updatelarvaspecies(hfly,hfly_extra,pos);
  case 'wingedfly_and_landmarks',
    updatewingedfly(hfly,hfly_extra(:,1,:),pos);
    updatelandmarks(hfly_extra(:,2,:),pos);
end
