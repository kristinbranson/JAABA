function DebugComputeWindowFeatures_SliderCallback(hObject,varargin)

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

ytrue = getappdata(hfig,'ytrue');
dlogprob = getappdata(hfig,'dlogprob');

ypred = 2+zeros(1,size(dlogprob,1));
ypred(dlogprob>=mindiff_logprob1) = 1;
ypred(dlogprob<=maxdiff_logprob0) = 0;

ycombine = ytrue*3 + ypred;
vs = setdiff(unique(ycombine),0);
for v = vs,
  [i0,i1] = get_interval_ends(ycombine==v);
  i1 = i1 - 1;
  xtmp = [i0-.5;i1+.5;i1+.5;i0-.5;i0-.5];
  ytmp = repmat(handles.ylim([1,1,2,2,1])',[1,numel(i0)]);
  set(handles.patch(v),'XData',xtmp,'YData',ytmp);
end

n = numel(ycombine);
legend(handles.patch,...
  sprintf('l=none, d=turn (%.3f)',nnz(ycombine==1)/n),...
  sprintf('l=none, d=unknown (%.3f)',nnz(ycombine==2)/n),...
  sprintf('l=turn, d=none (%.3f)',nnz(ycombine==3)/n),...
  sprintf('l=turn, d=turn (%.3f)',nnz(ycombine==4)/n),...
  sprintf('l=turn, d=unknown (%.3f)',nnz(ycombine==5)/n));
