function restorePointer(fig,oldPointer)

drawnow('update');  % make sure everything is updated before we restore the pointer
currentPointer=get(fig,'pointer');
if ~strcmpi(oldPointer,currentPointer)
  set(fig,'pointer',oldPointer);
end
drawnow('update');

end
