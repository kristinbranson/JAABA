function DebugComputeWindowFeatures_SliderCallback_Hist(hObject,varargin)

hfig = get(hObject,'Parent');
handles = getappdata(hfig,'handles');
maxdiff_logprob0 = get(handles.slider0,'Value');
mindiff_logprob1 = get(handles.slider1,'Value');
if mindiff_logprob1 < maxdiff_logprob0,
  maxdiff_logprob0 = (maxdiff_logprob0+mindiff_logprob1)/2;
  mindiff_logprob1 = maxdiff_logprob0;
  set(handles.slider0,'Value',maxdiff_logprob0);
  set(handles.slider1,'Value',mindiff_logprob1);
end
set(handles.text0,'String',sprintf('none<=%.1f',maxdiff_logprob0));
set(handles.text1,'String',sprintf('%.1f<=turn',mindiff_logprob1));

data = getappdata(hfig,'data');

ypred = 2+zeros(1,size(data.dlogprob,1));
ypred(data.dlogprob>=mindiff_logprob1) = 1;
ypred(data.dlogprob<=maxdiff_logprob0) = 0;

frac = nan([size(data.centers),4]);
for i = 1:size(data.centers,1),
  for v = 0:1,
    counts = histc(data.x(i,data.ytrue==v),data.edges(i,:));
    frac(i,:,2*v+1) = counts(1:end-1) / sum(counts(1:end-1));
  end
  for v = 0:1,
    counts = histc(data.x(i,ypred==v),data.edges(i,:));
    frac(i,:,2*v+2) = counts(1:end-1) / sum(counts(1:end-1));
  end
end

v = get(handles.listbox,'Value');
delete(handles.data(ishandle(handles.data)));
handles.data = bar(data.centers(v,:)',permute(frac(v,:,:),[2,3,1]),.99);
for i = 1:numel(handles.data),
  set(handles.data(i),'FaceColor',data.colors(i,:),'EdgeColor','none');
end
legend(handles.data,getappdata(handles.listbox,'legends'));
setappdata(handles.listbox,'frac',frac);
