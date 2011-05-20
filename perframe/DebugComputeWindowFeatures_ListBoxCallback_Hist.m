function DebugComputeWindowFeatures_ListBoxCallback_Hist(hObject,event)

v = get(hObject,'Value');
if numel(v) > 1,
  v = v(end);
  set(hObject,'Value',v);
end
frac = getappdata(hObject,'frac');
centers = getappdata(hObject,'centers');
handles = getappdata(hObject,'handles');

delete(handles.data(ishandle(handles.data)));

n = size(frac,3);
if n < 7,
  colors = lines(n);
else
  colors = jet(n)*.6;
end

xlim = [centers(v,1)*2-centers(v,2),centers(v,end)*2-centers(v,end-1)];
ylim = [0,max(frac(v,:))];

handles.data = bar(centers(v,:)',permute(frac(v,:,:),[2,3,1]),.99);
for i = 1:n,
  set(handles.data(i),'FaceColor',colors(i,:),'EdgeColor','none');
end
legend(handles.data,getappdata(hObject,'legends'));
set(handles.ax,'XLim',xlim,'YLim',ylim);
setappdata(hObject,'handles',handles);