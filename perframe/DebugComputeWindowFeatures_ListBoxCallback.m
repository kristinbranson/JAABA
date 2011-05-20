function DebugComputeWindowFeatures_ListBoxCallback(hObject,varargin)

v = get(hObject,'Value');
x = getappdata(hObject,'x');
handles = getappdata(hObject,'handles');

delete(handles.data(ishandle(handles.data)));
n = numel(v);
if n > 7,
  colors = jet(n)*.7;
else
  colors = lines(n);
end
ylim = [min(min(x(v,:))),max(max(x(v,:)))];
dy = diff(ylim);
ylim(1) = ylim(1)-dy*.01;
ylim(2) = ylim(2)+dy*.01;
if ylim(1) >= ylim(2),
  ylim(1) = min(ylim)-.5;
  ylim(2) = max(ylim)+.5;
end
set(handles.ax,'YLim',ylim);
handles.data = plot(handles.ax,x(v,:)','.-');
for i = 1:n,
  set(handles.data(i),'Color',colors(i,:));
end
setappdata(hObject,'handles',handles);