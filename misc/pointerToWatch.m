function oldPointer=pointerToWatch(fig)

oldPointer=get(fig,'pointer');
if ~strcmpi(oldPointer,'watch')
  set(fig,'pointer','watch');
end
drawnow('update');

end
