function updatewingedfly(hbody,hwings,pos)

% body position
if ~isempty(hbody),
  updatefly(hbody,pos.x,pos.y,pos.theta,pos.a,pos.b);
end

% wing positions
set(hwings,'XData',[pos.xwingl,pos.x,pos.xwingr],...
  'YData',[pos.ywingl,pos.y,pos.ywingr]);
