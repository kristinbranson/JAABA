function UpdateTargetPosition(targettype,hfly,hfly_extra,pos,hflies_apt,apt_info)
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
    updatewingedfly([],hfly_extra(:,1,:),pos);
    updatelandmarks(hfly,hfly_extra(:,2,:),pos);
end

if isfield(pos,'trk_x')
  % has APT trk 
  set(hflies_apt(1),'XData',pos.trk_x,'YData',pos.trk_y);
  if ~isempty(apt_info.skeletonEdges)
    sk = apt_info.skeletonEdges;
    curdat = nan(3*size(sk,1),2);
    curdat(1:3:end,1) = pos.trk_x(sk(:,1));
    curdat(2:3:end,1) = pos.trk_x(sk(:,2));
    curdat(1:3:end,2) = pos.trk_y(sk(:,1));
    curdat(2:3:end,2) = pos.trk_y(sk(:,2));
    set(hflies_apt(2),'XData',curdat(:,1),'YData',curdat(:,2));
    
  end
end