function centerOnParentFigure(fig,parentFig)

% Get our position
pos=get(fig,'position');
offset=pos(1:2);
sz=pos(3:4);

% Get out parent's position
parentPos=get(parentFig,'position');
parentOffset=parentPos(1:2);
parentSz=parentPos(3:4);

% Calculate a new offset that will center us on the parent
newOffset=parentOffset+(parentSz-sz)/2;

% Set our position, using the new offset but the same size as before
set(fig,'position',[newOffset sz]);

end
