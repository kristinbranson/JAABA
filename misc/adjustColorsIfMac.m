function adjustColorsIfMac(fig)

if ismac(),
  allpopups = findall(fig,'type','uicontrol','Style','popup');
  set(allpopups,'ForegroundColor',[0 0 0], ...
                'BackgroundColor',[1 1 1]);
  pushbuttons=findall(fig,'type','uicontrol','Style','pushbutton');
  adjustNonLabelButtonColor(pushbuttons)
  togglebuttons=findall(fig,'type','uicontrol','Style','togglebutton');
  adjustNonLabelButtonColor(togglebuttons)
end

end
