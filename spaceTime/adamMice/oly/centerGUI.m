function centerGUI(h)

OldUnits = get(h,'Units');
set(h,'Units','pixels');
OldPos = get(h,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);

ScreenUnits = get(0,'Units');
set(0,'Units','pixels');
ScreenSize = get(0,'MonitorPosition'); % MonitorPosition works better than ScreenSize for multiheaded comps
ScreenSize = ScreenSize(1,:);
set(0,'Units',ScreenUnits);

FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
FigPos(3:4)=[FigWidth FigHeight];
set(h,'Position',FigPos);
set(h,'Units',OldUnits);

end