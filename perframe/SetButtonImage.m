function SetButtonImage(hObj)
  bgColor = get(hObj,'BackgroundColor');
  buttonPos = get(hObj,'Position');
  buttonSz = round(buttonPos([4 3]));
  buttonImg = shiftdim(repmat(bgColor,[buttonSz(2) 1 buttonSz(1)]),2);
  if ismac 
    set(hObj,'CData',buttonImg);
  end
