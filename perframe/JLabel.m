function varargout = JLabel(varargin)
% JLABEL MATLAB code for JLabel.fig
%      JLABEL, by itself, creates a new JLABEL or raises the existing
%      singleton*.
%
%      H = JLABEL returns the handle to a new JLABEL or the handle to
%      the existing singleton*.
%
%      JLABEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JLABEL.M with the given input arguments.
%
%      JLABEL('Property','Value',...) creates a new JLABEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JLabel_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JLabel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JLabel

% Last Modified by GUIDE v2.5 23-Jan-2012 10:52:49

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JLabel_OpeningFcn, ...
                   'gui_OutputFcn',  @JLabel_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}) && exist(varargin{1}), %#ok<EXIST>
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before JLabel is made visible.
function JLabel_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JLabel (see VARARGIN)

% parse optional inputs
[handles.classifierfilename,...
  handles.configfilename,...
  handles.defaultpath] = ...
  myparse(varargin,...
  'classifierfilename','',...
  'configfilename','',...
  'defaultpath','');

% initialize statusbar
handles.status_bar_text = sprintf('Status: No experiment loaded');
handles.idlestatuscolor = [0,1,0];
handles.busystatuscolor = [1,0,1];
handles.movie_height = 100;
handles.movie_width = 100;
set(handles.similarFramesButton,'Enable','off');
handles.doFastUpdates = true;
ClearStatus(handles);

% read configuration
[handles,success] = LoadConfig(handles);
if ~success,
  guidata(hObject,handles);
  delete(handles.figure_JLabel);
  return;
end

% get relative locations of stuffs
handles = GetGUIPositions(handles);

% initialize data
handles = InitializeState(handles);

% initialize plot handles
handles = InitializePlots(handles);

% load classifier
if ~isempty(handles.classifierfilename),
  if exist(handles.classifierfilename,'file'),
    [success,msg] = handles.data.SetClassifierFileName(handles.classifierfilename);
    if ~success,
      warning(msg);
      SetStatus(handles,'Error loading classifier from file');
    end
  end
end

if isempty(handles.data.expdirs),
  guidata(hObject,handles);
  menu_file_editfiles_Callback(handles.figure_JLabel, [], handles);
  handles = guidata(hObject);
end

% keypress callback for all non-edit text objects
RecursiveSetKeyPressFcn(handles.figure_JLabel);

% enable gui
EnableGUI(handles);

handles.output = handles.figure_JLabel;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JLabel wait for user response (see UIRESUME)
% UNCOMMENT
%uiwait(handles.figure_JLabel);

% --- Outputs from this function are returned to the command line.
function varargout = JLabel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.figure_JLabel;
% UNCOMMENT
% if isfield(handles,'data'),
%   varargout{1} = handles.data;
% else
%   varargout{1} = [];
% end
% SaveRC(handles);
% delete(handles.figure_JLabel);

function handles = InitializePlots(handles)

handles.axes_preview_curr = 1;
if numel(handles.axes_previews) > numel(handles.ts),
  handles.ts = [handles.ts,repmat(handles.ts(end),[1,numel(handles.axes_previews)-numel(handles.ts)])];
end

% slider callbacks
for i = 1:numel(handles.slider_previews),
  fcn = get(handles.slider_previews(i),'Callback');
  %set(handles.slider_previews(i),'Callback','');
  handles.hslider_listeners(i) = handle.listener(handles.slider_previews(i),...
    'ActionEvent',fcn);
end

% fly current positions
handles.hflies = zeros(handles.nflies_curr,numel(handles.axes_previews));
handles.hflies_extra = zeros(handles.nflies_curr,numel(handles.axes_previews));
handles.hfly_markers = zeros(handles.nflies_curr,numel(handles.axes_previews));
% fly path
handles.htrx = zeros(handles.nflies_label,numel(handles.axes_previews));

% choose colors for flies
% TODO: change hard-coded colormap
handles.fly_colors = jet(handles.nflies_curr)*.7;
handles.fly_colors = handles.fly_colors(randperm(handles.nflies_curr),:);

handles.hlabel_curr = nan(1,numel(handles.axes_previews));
for i = 1:numel(handles.axes_previews),
  cla(handles.axes_previews(i),'reset');

  % image in axes_preview
  handles.himage_previews(i) = imagesc(0,'Parent',handles.axes_previews(i),[0,255]);
  set(handles.himage_previews(i),'HitTest','off');
  axis(handles.axes_previews(i),'equal');
  
  set(handles.axes_previews(i),'ButtonDownFcn',@(hObject,eventdata) JLabel('axes_preview_ButtonDownFcn',hObject,eventdata,guidata(hObject)));
  hold(handles.axes_previews(i),'on');

  % labeled behaviors
  handles.hlabels = nan(1,handles.data.nbehaviors);
  handles.hpredicted = nan(1,handles.data.nbehaviors);
  handles.hlabelstarts = nan(1,handles.data.nbehaviors);
  for j = 1:handles.data.nbehaviors,
    handles.hlabels(j) = plot(handles.axes_previews(i),nan,nan,'-',...
      'color',handles.labelcolors(j,:),'linewidth',5,'HitTest','off');
    handles.hpredicted(j) = plot(handles.axes_previews(i),nan,nan,'-',...
      'color',handles.labelcolors(j,:),'linewidth',5,'HitTest','off');
    % start of label
    handles.hlabelstarts(j) = plot(handles.axes_previews(i),nan,nan,'v',...
      'color',handles.labelcolors(j,:),'markerfacecolor',handles.labelcolors(j,:),...
      'HitTest','off');
    
    set(handles.axes_previews(i),'Color','k','XColor','w','YColor','w');
    
  end
  
  if handles.plot_labels_manual,
    set(handles.hlabels,'Visible','on');
  else
    set(handles.hlabels,'Visible','off');
  end
  if handles.plot_labels_automatic,
    set(handles.hpredicted,'Visible','on');
  else
    set(handles.hpredicted,'Visible','off');
  end
  
  % current label plotted on axes
  handles.hlabel_curr(i) = plot(nan(1,2),nan(1,2),'k-',...
    'Parent',handles.axes_previews(i),...
    'HitTest','off','Linewidth',5);
  
  % trx of flies
  for j = 1:handles.nflies_label,
    handles.htrx(j,i) = plot(handles.axes_previews(i),nan,nan,'.-',...
      'linewidth',1,'HitTest','off');
  end
  
  % fly current positions
  for fly = 1:handles.nflies_curr,
    handles.hflies(fly,i) = plot(handles.axes_previews(i),nan,nan,'-',...
      'color',handles.fly_colors(fly,:),'linewidth',3,...
      'ButtonDownFcn',@(hObject,eventdata) JLabel('fly_ButtonDownFcn',hObject,eventdata,guidata(hObject),fly,i));
    handles.hflies_extra(fly,i) = plot(handles.axes_previews(i),nan,nan,'o',...
      'color',handles.fly_colors(fly,:),'markerfacecolor',handles.fly_colors(fly,:),...
      'MarkerSize',handles.flies_extra_markersize,...
      'ButtonDownFcn',@(hObject,eventdata) JLabel('fly_ButtonDownFcn',hObject,eventdata,guidata(hObject),fly,i));
    handles.hfly_markers(fly,i) = plot(handles.axes_previews(i),nan,nan,'*',...
      'color',handles.fly_colors(fly,:),'linewidth',3,...
      'ButtonDownFcn',@(hObject,eventdata) JLabel('fly_ButtonDownFcn',hObject,eventdata,guidata(hObject),fly,i),...
      'Visible','off');
  end

end

% TODO: allow colormap options
colormap(handles.axes_preview,gray(256));

% timelines

% zoom
handles.hzoom = zoom(handles.figure_JLabel);
handles.hpan = pan(handles.figure_JLabel);
set(handles.hzoom,'ActionPostCallback',@(hObject,eventdata) PostZoomCallback(hObject,eventdata,guidata(eventdata.Axes)));
set(handles.hpan,'ActionPostCallback',@(hObject,eventdata) PostZoomCallback(hObject,eventdata,guidata(eventdata.Axes)));

% manual timeline
timeline_axes_color = get(handles.panel_timelines,'BackgroundColor');
handles.himage_timeline_manual = image(zeros([1,1,3]),'Parent',handles.axes_timeline_manual);
set(handles.himage_timeline_manual,'HitTest','off');
hold(handles.axes_timeline_manual,'on');
handles.htimeline_manual_starts = plot(handles.axes_timeline_manual,nan,nan,'w-','HitTest','off');
ylim = [.5,1.5];
ydata = [ylim(1)+diff(ylim)*.025,ylim(2)-diff(ylim)*.025];
handles.htimeline_label_curr = patch(nan(1,5),ydata([1,2,2,1,1]),'k',...
  'Parent',handles.axes_timeline_manual,'LineStyle','--','EdgeColor','w',...
  'HitTest','off','Linewidth',3,'Clipping','on');
if handles.plot_labels_manual,
  set(handles.timeline_label_manual,'ForegroundColor',handles.emphasiscolor,'FontWeight','bold');
else
  set(handles.timeline_label_manual,'ForegroundColor',handles.unemphasiscolor,'FontWeight','normal');
end

set(handles.axes_timeline_manual,'YTick',[]);
setAxesZoomMotion(handles.hzoom,handles.axes_timeline_manual,'horizontal');
setAllowAxesPan(handles.hpan,handles.axes_timeline_manual,false);

% auto timeline
ydata_im = [2/3,4/3];
handles.himage_timeline_auto = image(zeros([3,1,3]),'Parent',handles.axes_timeline_auto);
set(handles.himage_timeline_auto,'YData',ydata_im);
set(handles.himage_timeline_auto,'HitTest','off');
hold(handles.axes_timeline_auto,'on');
handles.htimeline_auto_starts = plot(handles.axes_timeline_auto,nan,nan,'w-','HitTest','off');
set(handles.axes_timeline_auto,'YTick',[]);
setAxesZoomMotion(handles.hzoom,handles.axes_timeline_auto,'horizontal');
setAllowAxesPan(handles.hpan,handles.axes_timeline_auto,false);
if handles.plot_labels_automatic,
  set(handles.timeline_label_automatic,'ForegroundColor',handles.emphasiscolor,'FontWeight','bold');
else
  set(handles.timeline_label_automatic,'ForegroundColor',handles.unemphasiscolor,'FontWeight','normal');
end

for h = handles.axes_timeline_labels,
  set(h,'YLim',[.5,1.5]);
end

% properties
propi = 1;
handles.htimeline_data(propi) = plot(handles.axes_timeline_props(propi),nan,nan,'w.-','HitTest','off');
hold(handles.axes_timeline_props(propi),'on');

% whether the manual and auto match
handles.htimeline_errors = plot(handles.axes_timeline_manual,nan,nan,'-',...
  'color',handles.incorrectcolor,'HitTest','off','Linewidth',5);
% new suggestions
handles.htimeline_suggestions = plot(handles.axes_timeline_manual,nan,nan,'-',...
  'color',handles.suggestcolor,'HitTest','off','Linewidth',5);

handles.menu_view_zoom_options = setdiff(findall(handles.menu_view_zoom,'Type','uimenu'),...
  handles.menu_view_zoom);

% suggest timeline
% handles.himage_timeline_suggest = image(zeros([1,1,3]),'Parent',handles.axes_timeline_suggest);
% set(handles.himage_timeline_suggest,'HitTest','off');
% hold(handles.axes_timeline_suggest,'on');
% handles.htimeline_suggest_starts = plot(handles.axes_timeline_suggest,nan,nan,'w-','HitTest','off');

% error timeline
% handles.himage_timeline_error = image(zeros([1,1,3]),'Parent',handles.axes_timeline_error);
% set(handles.himage_timeline_error,'HitTest','off');
% hold(handles.axes_timeline_error,'on');
% handles.htimeline_error_starts = plot(handles.axes_timeline_error,nan,nan,'w-','HitTest','off');

for i = 1:numel(handles.axes_timeline_props),
  setAxesZoomMotion(handles.hzoom,handles.axes_timeline_props(i),'vertical');
  setAllowAxesPan(handles.hpan,handles.axes_timeline_props(i),true);
  setAxesPanMotion(handles.hpan,handles.axes_timeline_props(i),'vertical');
end

for i = 1:numel(handles.axes_timelines),
  hold(handles.axes_timelines(i),'on');
  set(handles.axes_timelines(i),'XColor','w','YColor','w','Color',timeline_axes_color);
end

handles.hcurr_timelines = nan(size(handles.axes_timelines));
for i = 1:numel(handles.axes_timelines),
  handles.hcurr_timelines(i) = plot(handles.axes_timelines(i),nan(1,2),[-10^6,10^6],'y-','HitTest','off','linewidth',2);
end
handles.hselection = nan(size(handles.axes_timelines));
for i = 1:numel(handles.axes_timelines),
  ylim = [.5,1.5];
  ydata = [ylim(1)+diff(ylim)*.025,ylim(2)-diff(ylim)*.025];
  handles.hselection(i) = ...
    plot(handles.axes_timelines(i),nan(1,5),ydata([1,2,2,1,1]),'--','color',handles.selection_color,...
    'HitTest','off','Linewidth',3);
end

for i = 2:numel(handles.axes_timelines),
%  if handles.axes_timelines(i) ~= handles.axes_timeline_error,
%  if handles.axes_timelines(i) ~= handles.axes_timeline_auto,
    set(handles.axes_timelines(i),'XTickLabel',{});
%  end
end

linkaxes(handles.axes_timelines,'x');

%for i = 1:numel(handles.axes_timelines),
%  setAxesZoomMotion(handles.hzoom,handles.axes_timelines(i),'horizontal');
%end

% timeline callbacks
% fcn = @(hObject,eventdata) JLabel('axes_timeline_ButtonDownFcn',hObject,eventdata,guidata(hObject));
% for i = 1:numel(handles.axes_timelines),
%   set(handles.axes_timelines(i),'ButtonDownFcn',fcn);
% end


function UpdatePlots(handles,varargin)

% WARNING: we directly access handles.data.trx for speed here -- 
% REMOVED! NOT SO SLOW

[axes,refreshim,refreshflies,refreshtrx,refreshlabels,...
  refresh_timeline_manual,refresh_timeline_auto,refresh_timeline_suggest,refresh_timeline_error,...
  refresh_timeline_xlim,refresh_timeline_hcurr,...
  refresh_timeline_props,refresh_timeline_selection,...
  refresh_curr_prop] = ...
  myparse(varargin,'axes',1:numel(handles.axes_previews),...
  'refreshim',true,'refreshflies',true,'refreshtrx',true,'refreshlabels',true,...
  'refresh_timeline_manual',true,...
  'refresh_timeline_auto',true,...
  'refresh_timeline_suggest',true,...
  'refresh_timeline_error',true,...
  'refresh_timeline_xlim',true,...
  'refresh_timeline_hcurr',true,...
  'refresh_timeline_props',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',true);

% make sure data for this experiment is loaded
% if handles.expi ~= handles.data.expi,
%   SetStatus('Preloading data for experiment %s, flies %s',handles.data.expnames{handles.expi},mat2str(handles.flies));
%   handles.data.PreLoad(handles.expi,handles.flies);
% end

% update timelines
if refresh_timeline_manual,
  set(handles.himage_timeline_manual,'CData',handles.labels_plot.im);
  tmp = find(handles.labels_plot.isstart);
  nstarts = numel(tmp);
  tmpx = reshape(cat(1,repmat(tmp,[2,1]),nan(1,nstarts)),[3*nstarts,1]);
  tmpy = reshape(repmat([.5;1.5;nan],[1,nstarts]),[3*nstarts,1]);
  set(handles.htimeline_manual_starts,'XData',tmpx,'YData',tmpy);  
  if handles.label_state ~= 0,
    ts = sort([handles.label_t0,handles.ts(1)]) + [-.5,.5];
    set(handles.htimeline_label_curr,'XData',ts([1,1,2,2,1]));
  end
end

if refresh_timeline_auto,
  set(handles.himage_timeline_auto,'CData',handles.labels_plot.predicted_im);
end
if refresh_timeline_suggest,
  set(handles.htimeline_suggestions,'XData',handles.labels_plot.suggest_xs,...
    'YData',zeros(size(handles.labels_plot.suggest_xs))+1.5);
  %set(handles.himage_timeline_suggest,'CData',handles.labels_plot.suggested_im);
end
if refresh_timeline_error,
  set(handles.htimeline_errors,'XData',handles.labels_plot.error_xs,...
  'YData',zeros(size(handles.labels_plot.error_xs))+1.5);
  %set(handles.himage_timeline_error,'CData',handles.labels_plot.error_im);
end

if refresh_timeline_xlim,
  xlim = [handles.ts(1)-(handles.timeline_nframes-1)/2,...
    handles.ts(1)+(handles.timeline_nframes-1)/2];
  for i = 1:numel(handles.axes_timelines),
    set(handles.axes_timelines(i),'XLim',xlim);
    %zoom(handles.axes_timelines(i),'reset');
  end
end


if refresh_timeline_hcurr,
  set(handles.hcurr_timelines,'XData',handles.ts([1,1]));
end
if refresh_timeline_selection,
  tmp = handles.selected_ts + .5*[-1,1];
  set(handles.hselection,'XData',tmp([1,1,2,2,1]));
end

if refresh_timeline_props,
  for propi = 1:numel(handles.perframepropis),
    v = handles.perframepropis(propi);
    [perframedata,T0,T1] = handles.data.GetPerFrameData(handles.expi,handles.flies,v);
    set(handles.htimeline_data(propi),'XData',T0:T1,...
      'YData',perframedata);
    %if isnan(handles.timeline_data_ylims(1,v)),
      ylim = [min(perframedata),max(perframedata)];
      set(handles.axes_timeline_props(propi),'YLim',ylim);
      zoom(handles.axes_timeline_props(propi),'reset');
    %end
    if ~isnan(handles.timeline_data_ylims(1,v)),
      ylim = handles.timeline_data_ylims(:,v);
      set(handles.axes_timeline_props(propi),'YLim',ylim);
    end
    ydata = [ylim(1)+diff(ylim)*.025,ylim(2)-diff(ylim)*.025];
    set(handles.hselection(propi),'YData',ydata([1,2,2,1,1]));      
  end
end

%drawnow;

for i = axes,
  
  if refreshim,

    % read in current frame
    %image_cache = getappdata(handles.figure_JLabel,'image_cache');
%     try
%       j = find(handles.ts(i)==image_cache.ts,1);
%       if isempty(j),
        im = handles.readframe(handles.ts(i));
%       else
%         im = image_cache.ims(:,:,:,j);
%       end
%     catch ME
%       uiwait(warndlg(sprintf('Could not read frame %d from current movie: %s',handles.ts(i),getReport(ME))));
%       return;
%     end
  
    % update frame
    set(handles.himage_previews(i),'CData',im);

%     if isempty(j),
%       j = argmin(image_cache.timestamps);
%       image_cache.ims(:,:,:,j) = im;
%       image_cache.ts(j) = handles.ts(i);
%     end
%     image_cache.timestamps(j) = now;
%     setappdata(handles.figure_JLabel,'image_cache',image_cache);
    
  end
  
  % TODO: change hard-coded linewidths, nprev, npost
  
  % update current position
  if refreshflies,
    if handles.ts(i) < handles.t0_curr || handles.ts(i) > handles.t1_curr,
      labelidx = [];
    elseif handles.label_state ~= 0,
      labelidx = handles.label_state;
    elseif handles.plot_labels_manual,
      labelidxStruct = handles.data.GetLabelIdx(handles.expi,handles.flies,handles.ts(i),handles.ts(i));
      labelidx = labelidxStruct.vals;
    elseif handles.plot_labels_automatic,
       prediction = handles.data.GetPredictedIdx(handles.expi,handles.flies,handles.ts(i),handles.ts(i));
       labelidx = prediction.predictedidx;
    end
    inbounds = handles.data.firstframes_per_exp{handles.expi} <= handles.ts(i) & ...
      handles.data.endframes_per_exp{handles.expi} >= handles.ts(i);
    set(handles.hflies(~inbounds,i),'XData',nan,'YData',nan);
    set(handles.hflies_extra(~inbounds,i),'XData',nan,'YData',nan);
    set(handles.hfly_markers(~inbounds,i),'XData',nan,'YData',nan);
    for fly = find(inbounds),
      % WARNING: this accesses handles.data.trx directly -- make sure that
      % handles.data.trx is loaded for the correct movie
      % REMOVED! NOT SO SLOW

      t = handles.ts(i);
      pos = handles.data.GetTrxPos1(handles.expi,fly,t);
      UpdateTargetPosition(handles.data.targettype,handles.hflies(fly,i),handles.hflies_extra(fly,i),pos);
      set(handles.hfly_markers(fly,i),'XData',pos.x,'YData',pos.y);
      sexcurr = handles.data.GetSex1(handles.expi,fly,t);
      if lower(sexcurr(1)) == 'm',
        set(handles.hfly_markers(fly,i),'Visible','on');
      else
        set(handles.hfly_markers(fly,i),'Visible','off');
      end
%       updatefly(handles.hflies(fly,i),...
%         handles.data.GetTrxX1(handles.expi,fly,t),...
%         handles.data.GetTrxY1(handles.expi,fly,t),...
%         handles.data.GetTrxTheta1(handles.expi,fly,t),...
%         handles.data.GetTrxA1(handles.expi,fly,t),...
%         handles.data.GetTrxB1(handles.expi,fly,t));
%       j = handles.ts(i) + handles.data.trx(fly).off;
%       updatefly(handles.hflies(fly,i),handles.data.trx(fly).x(j),...
%         handles.data.trx(fly).y(j),...
%         handles.data.trx(fly).theta(j),...
%         handles.data.trx(fly).a(j),...
%         handles.data.trx(fly).b(j));
      %updatefly(handles.hflies(fly,i),trx(fly).x,trx(fly).y,trx(fly).theta,trx(fly).a,trx(fly).b);
      if ismember(fly,handles.flies),
        set(handles.hflies(fly,i),'LineWidth',3);
        if labelidx <= 0,
          set(handles.hflies(fly,i),'Color',handles.labelunknowncolor);
          set(handles.hflies_extra(fly,i),'Color',handles.labelunknowncolor,...
            'MarkerFaceColor',handles.labelunknowncolor);
        else
          set(handles.hflies(fly,i),'Color',handles.labelcolors(labelidx,:),...
            'MarkerFaceColor',handles.labelcolors(labelidx,:));
        end
      else
        set(handles.hflies(fly,i),'LineWidth',1);
      end
    end
    
    if strcmpi(handles.preview_zoom_mode,'center_on_fly'),
      ZoomInOnFlies(handles,i);
    elseif strcmpi(handles.preview_zoom_mode,'follow_fly'),
      KeepFliesInView(handles,i);
    end    
  end

  % update trx
  % TODO: remove hard-coded nprev, npost
  nprev = handles.traj_nprev;
  npost = handles.traj_npost;
  if refreshtrx,
    for j = 1:numel(handles.flies),
      fly = handles.flies(j);
      tmp = handles.ts(i);
      t0 = handles.data.firstframes_per_exp{handles.expi}(fly);
      t1 = handles.data.endframes_per_exp{handles.expi}(fly);
      ts = max(t0,tmp-nprev):min(t1,tmp+npost);
      set(handles.htrx(j,i),'XData',handles.data.GetTrxValues('X1',handles.expi,fly,ts),...
        'YData',handles.data.GetTrxValues('Y1',handles.expi,fly,ts));
      %j0 = max(1,tmp-nprev);
      %j1 = min(handles.data.trx(fly).nframes,tmp+npost);
      %set(handles.htrx(j,i),'XData',handles.data.trx(fly).x(j0:j1),...
      %  'YData',handles.data.trx(fly).y(j0:j1));
      %trx = handles.data.GetTrx(handles.expi,fly,handles.ts(i)-nprev:handles.ts(i)+npost);
      %set(handles.htrx(j,i),'XData',trx.x,'YData',trx.y);
    end
  end  
  
  % update labels plotted
  if refreshlabels,
    for k = 1:numel(handles.flies),
      fly = handles.flies(k);
      T0 = handles.data.firstframes_per_exp{handles.expi}(fly);
      T1 = handles.data.endframes_per_exp{handles.expi}(fly);
%       T0 = handles.data.GetTrxFirstFrame(handles.expi,fly);
%       T1 = handles.data.GetTrxEndFrame(handles.expi,fly);
      t0 = min(T1,max(T0,handles.ts(i)-nprev));
      t1 = min(T1,max(T0,handles.ts(i)+npost));
      for j = 1:handles.data.nbehaviors,
        xplot = handles.labels_plot.x(:,handles.labels_plot_off+t0:handles.labels_plot_off+t1,j,k);
        yplot = handles.labels_plot.y(:,handles.labels_plot_off+t0:handles.labels_plot_off+t1,j,k);
        set(handles.hlabels(j),'XData',xplot(:),'YData',yplot(:));
        xpred = handles.labels_plot.predx(:,handles.labels_plot_off+t0:handles.labels_plot_off+t1,j,k);
        ypred = handles.labels_plot.predy(:,handles.labels_plot_off+t0:handles.labels_plot_off+t1,j,k);
        set(handles.hpredicted(j),'XData',xpred(:),'YData',ypred(:));
      end
      if handles.label_state ~= 0,
        ts = sort([handles.label_t0,handles.ts(1)]);
        t0 = max(t0,ts(1));
        t1 = min(t1,ts(2)+1);
        xdata = handles.data.GetTrxValues('X1',handles.expi,handles.flies(1),t0:t1);
        ydata = handles.data.GetTrxValues('Y1',handles.expi,handles.flies(1),t0:t1);
        set(handles.hlabel_curr(1),'XData',xdata,'YData',ydata);
        if handles.label_state == -1,
          set(handles.hlabel_curr(1),'Color',handles.labelunknowncolor);
        else
          set(handles.hlabel_curr(1),'Color',handles.labelcolors(handles.label_state,:));
        end
      else
        set(handles.hlabel_curr(1),'XData',nan,'YData',nan);
      end

    end
  end
  
  if refresh_curr_prop,
    for propi = 1:numel(handles.perframepropis),
      v = handles.perframepropis(propi);
      if handles.ts(i) < handles.t0_curr || handles.ts(i) > handles.t1_curr,
        s = '';
      else
        perframedata = handles.data.GetPerFrameData1(handles.expi,handles.flies,v,handles.ts(i));
        s = sprintf('%.3f',perframedata);
      end
      if numel(handles.text_timeline_props) >= propi && ishandle(handles.text_timeline_props(propi)),
        set(handles.text_timeline_props(propi),'String',s);
      end
    end
  end
  
  %drawnow;
  
end


function [handles,success] = SetCurrentMovie(handles,expi)

success = false;

if expi == handles.expi,
  success = true;
  return;
end

% check that the current movie exists
moviefilename = handles.data.GetFile('movie',expi);
if ~exist(moviefilename,'file'),
  uiwait(warndlg(sprintf('Movie file %s does not exist.',moviefilename)),'Error setting movie');
  return;
end

% close previous movie
if isfield(handles,'movie_fid') && ~isempty(fopen(handles.movie_fid)),
  fclose(handles.movie_fid);
end

% open new movie
% try
  SetStatus(handles,'Opening movie...');
  [handles.readframe,handles.nframes,handles.movie_fid,handles.movieheaderinfo] = ...
    get_readframe_fcn(moviefilename,'interruptible',false);
  im = handles.readframe(1);
  handles.movie_width = size(im,2);
  handles.movie_height = size(im,1);
% catch ME,
%   uiwait(warndlg(sprintf('Error opening movie file %s: %s',moviefilename,getReport(ME)),'Error setting movie'));
%   ClearStatus(handles);
%   return;
% end

% number of flies
handles.nflies_curr = handles.data.nflies_per_exp(expi);

% choose flies
if handles.nflies_curr == 0,
  flies = [];
else
  flies = 1;
end

% load trx
[success,msg] = handles.data.PreLoad(expi,flies);
if ~success,
  uiwait(errordlg(sprintf('Error loading data for experiment %d: %s',expi,msg)));
  return;
end

% set zoom radius
if isnan(handles.zoom_fly_radius(1)),
  handles.zoom_fly_radius = nanmean([handles.data.trx.a])*20 + [0,0];
end


handles.expi = expi;

ClearStatus(handles);

% TODO: change hard-coded colormap
% update colors
handles.fly_colors = jet(handles.nflies_curr)*.7;
handles.fly_colors = handles.fly_colors(randperm(handles.nflies_curr),:);

% delete old fly current positions
if isfield(handles,'hflies'),
  delete(handles.hflies(ishandle(handles.hflies)));
  handles.hflies = [];
end
if isfield(handles,'hflies_extra'),
  delete(handles.hflies_extra(ishandle(handles.hflies_extra)));
  handles.hflies_extra = [];
end
if isfield(handles,'hfly_markers'),
  delete(handles.hfly_markers(ishandle(handles.hfly_markers)));
  handles.hfly_markers = [];
end

% update plotted trx handles, as number of flies will change
handles.hflies = zeros(handles.nflies_curr,numel(handles.axes_previews));
handles.hflies_extra = zeros(handles.nflies_curr,numel(handles.axes_previews));
handles.hfly_markers = zeros(handles.nflies_curr,numel(handles.axes_previews));

for i = 1:numel(handles.axes_previews),
  % fly current positions
  for fly = 1:handles.nflies_curr,
    handles.hflies(fly,i) = plot(handles.axes_previews(i),nan,nan,'-',...
      'color',handles.fly_colors(fly,:),'linewidth',3,...
      'ButtonDownFcn',@(hObject,eventdata) JLabel('fly_ButtonDownFcn',hObject,eventdata,guidata(hObject),fly,i));
    handles.hflies_extra(fly,i) = plot(handles.axes_previews(i),nan,nan,'o',...
      'color',handles.fly_colors(fly,:),'MarkerFaceColor',handles.fly_colors(fly,:),...
      'MarkerSize',handles.flies_extra_markersize,...
      'ButtonDownFcn',@(hObject,eventdata) JLabel('fly_ButtonDownFcn',hObject,eventdata,guidata(hObject),fly,i));
    handles.hfly_markers(fly,i) = plot(handles.axes_previews(i),nan,nan,'*',...
      'color',handles.fly_colors(fly,:),'linewidth',3,...
      'ButtonDownFcn',@(hObject,eventdata) JLabel('fly_ButtonDownFcn',hObject,eventdata,guidata(hObject),fly,i),...
      'Visible','off');
  end
end

% set flies
handles = SetCurrentFlies(handles,flies,true,false);

% update slider steps, range
for i = 1:numel(handles.slider_previews),
  set(handles.slider_previews(i),'Min',1,'Max',handles.nframes,...
    'Value',1,...
    'SliderStep',[1/(handles.nframes-1),100/(handles.nframes-1)]);
end

% choose frame
for i = 1:numel(handles.axes_previews),
  handles = SetCurrentFrame(handles,i,handles.t0_curr,nan,true,false);
end

% update zoom
for i = 1:numel(handles.axes_previews),
  axis(handles.axes_previews(i),[.5,handles.movie_width+.5,.5,handles.movie_height+.5]);
end
for i = 1:numel(handles.axes_previews),
  zoom(handles.axes_previews(i),'reset');
end

% update plot
UpdatePlots(handles,'refresh_timeline_props',true,'refresh_timeline_selection',true);

for h = handles.axes_timeline_labels,
  zoom(h,'reset');
end

% enable GUI components
EnableGUI(handles);

success = true;

function handles = UnsetCurrentMovie(handles)

% close previous movie
if isfield(handles,'movie_fid') && ~isempty(fopen(handles.movie_fid)),
  fclose(handles.movie_fid);
end

handles.expi = 0;
handles.flies = nan(1,handles.nflies_label);
handles.ts = zeros(1,numel(handles.axes_previews));
handles.label_state = 0;
handles.label_imp = [];
handles.nflies_curr = 0;
% delete old fly current positions
if isfield(handles,'hflies'),
  delete(handles.hflies(ishandle(handles.hflies)));
  handles.hflies = [];
end
if isfield(handles,'hflies_extra'),
  delete(handles.hflies_extra(ishandle(handles.hflies_extra)));
  handles.hflies_extra = [];
end
if isfield(handles,'hfly_markers'),
  delete(handles.hfly_markers(ishandle(handles.hfly_markers)));
  handles.hfly_markers = [];
end


% enable GUI components
EnableGUI(handles);


function i = GetPreviewPanelNumber(hObject)

i = regexp(get(get(hObject,'Parent'),'Tag'),'^panel_axes(\d+)$','tokens','once');
if isempty(i),
  warning('Could not find index of parent panel');
  i = 1;
else
  i = str2double(i{1});
end


% --- Executes on slider movement.
function slider_preview_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to slider_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% get slider value
t = min(max(1,round(get(hObject,'Value'))),handles.nframes);

% which preview panel is this
i = GetPreviewPanelNumber(hObject);

% set current frame
SetCurrentFrame(handles,i,t,hObject);

function handles = SetCurrentFlies(handles,flies,doforce,doupdateplot)

if ~exist('doforce','var'),
  doforce = false;
end
if ~exist('doupdateplot','var'),
  doupdateplot = true;
end

[success,msg] = handles.data.PreLoad(handles.expi,flies);
if ~success,
  uiwait(waitdlg(sprintf('Error loading data for current set of flies: %s',msg)));
  return;
end

% same flies, return
if ~doforce && isempty(setdiff(flies,handles.flies)) && ...
    isempty(setdiff(handles.flies,flies)),
  return;
end

handles.flies = flies;

% frames these flies are both alive
handles.t0_curr = max(handles.data.GetTrxFirstFrame(handles.expi,handles.flies));
handles.t1_curr = min(handles.data.GetTrxEndFrame(handles.expi,handles.flies));

% form of labels for easier plotting:
% x, y positions of all labels
handles.labels_plot = struct;
n = handles.t1_curr-handles.t0_curr+1;
handles.labels_plot.im = zeros([1,n,3]);
handles.labels_plot.predicted_im = zeros([1,n,3]);
handles.labels_plot.suggest_xs = nan;
handles.labels_plot.error_xs = nan;
%handles.labels_plot.suggested_im = zeros([1,n,3]);
%handles.labels_plot.error_im = zeros([1,n,3]);
handles.labels_plot.x = nan(2,n,handles.data.nbehaviors,numel(handles.flies));
handles.labels_plot.y = nan(2,n,handles.data.nbehaviors,numel(handles.flies));
handles.labels_plot.predx = nan(2,n,handles.data.nbehaviors,numel(handles.flies));
handles.labels_plot.predy = nan(2,n,handles.data.nbehaviors,numel(handles.flies));
handles.labels_plot_off = 1-handles.t0_curr;
% set([handles.himage_timeline_manual,handles.himage_timeline_auto,...
%   handles.himage_timeline_error,handles.himage_timeline_suggest],...
%   'XData',[handles.t0_curr,handles.t1_curr]);
set([handles.himage_timeline_manual,handles.himage_timeline_auto],...
  'XData',[handles.t0_curr,handles.t1_curr]);

labelidxStruct = handles.data.GetLabelIdx(handles.expi,flies);
labelidx = labelidxStruct.vals;

prediction = handles.data.GetPredictedIdx(handles.expi,flies);
predictedidx = prediction.predictedidx;
scores = handles.data.NormalizeScores(prediction.scoresidx);
for flyi = 1:numel(flies),
  fly = flies(flyi);
  x = handles.data.GetTrxValues('X1',handles.expi,fly,handles.t0_curr:handles.t1_curr);
  y = handles.data.GetTrxValues('Y1',handles.expi,fly,handles.t0_curr:handles.t1_curr);
  for behaviori = 1:handles.data.nbehaviors
    idx = find(labelidx == behaviori);
    idx1 = min(idx+1,numel(x));
    handles.labels_plot.x(1,idx,behaviori,flyi) = x(idx);
    handles.labels_plot.x(2,idx,behaviori,flyi) = x(idx1);
    handles.labels_plot.y(1,idx,behaviori,flyi) = y(idx);
    handles.labels_plot.y(2,idx,behaviori,flyi) = y(idx1);
    
%     idx = find(predictedidx == behaviori);
    idx = find((predictedidx == behaviori) & ...
      (abs(scores)>handles.data.GetConfidenceThreshold(behaviori)));
    idx1 = min(idx+1,numel(x));
    handles.labels_plot.predx(1,idx,behaviori,flyi) = x(idx);
    handles.labels_plot.predx(2,idx,behaviori,flyi) = x(idx1);
    handles.labels_plot.predy(1,idx,behaviori,flyi) = y(idx);
    handles.labels_plot.predy(2,idx,behaviori,flyi) = y(idx1);
  end
end
handles = UpdateTimelineIms(handles);

% which interval we're currently within
handles.current_interval = [];

% update timelines
set(handles.himage_timeline_manual,'CData',handles.labels_plot.im);
axis(handles.axes_timeline_manual,[handles.t0_curr-.5,handles.t1_curr+.5,.5,1.5]);
% update zoom
for h = handles.axes_timeline_labels,
  zoom(h,'reset');
end

% update trx colors
for i = 1:numel(handles.axes_previews),
  for j = 1:numel(handles.flies),
    fly = handles.flies(j);
    set(handles.htrx(j,i),'Color',handles.fly_colors(fly,:));
  end
end

% Update colors for all other flies. 
inbounds = handles.data.firstframes_per_exp{handles.expi} <= handles.ts(i) & ...
  handles.data.endframes_per_exp{handles.expi} >= handles.ts(i);

for i = 1:numel(handles.axes_previews),
  for fly = find(inbounds),
    set(handles.hflies(fly,i),'Color',handles.fly_colors(fly,:));
    set(handles.hflies_extra(fly,i),'Color',handles.fly_colors(fly,:),...
      'MarkerFaceColor',handles.fly_colors(fly,:));
  end
end

% status bar text
[~,expname] = myfileparts(handles.data.expdirs{handles.expi});
if numel(handles.flies) == 1,
  handles.status_bar_text = sprintf('%s, fly %d',expname,handles.flies);
else
  handles.status_bar_text = [sprintf('%s, flies',expname),sprintf(' %d',handles.flies)];
end

% make sure frame is within bounds
isset = handles.ts ~= 0;
ts = max(handles.t0_curr,min(handles.t1_curr,handles.ts));
ts(~isset) = 0;
for i = 1:numel(ts),
  if ts(i) ~= handles.ts(i),
    handles = SetCurrentFrame(handles,i,ts(i),nan);
    handles.ts(i) = ts(i);
  end
end

ClearStatus(handles);

% TODO: generalize to multiple flies
s = GetTargetInfo(handles,flies(1));
set(handles.text_selection_info,'String',s);

guidata(handles.figure_JLabel,handles);

if doupdateplot,
  UpdatePlots(handles,'refresh_timeline_props',true,'refresh_timeline_selection',true);
end

function handles = UpdateTimelineIms(handles)

% Note: this function directly accesses handles.data.labelidx,
% handles.data.predictedidx for speed, so make sure we've preloaded the
% right experiment, flies
% REMOVED!

% if handles.expi ~= handles.data.expi || ~all(handles.flies == handles.data.flies),
%   handles.data.Preload(handles.expi,handles.flies);
% end

handles.labels_plot.im(:) = 0;
labelidx = handles.data.GetLabelIdx(handles.expi,handles.flies);

for behaviori = 1:handles.data.nbehaviors
  idx = (labelidx.vals == behaviori) & labelidx.imp;
  curColor = ShiftColor.increaseIntensity(handles.labelcolors(behaviori,:));
  for channel = 1:3,
    handles.labels_plot.im(1,idx,channel) = curColor(channel);
  end
  
  idx = (labelidx.vals == behaviori) & ~labelidx.imp;
  curColor = handles.labelcolors(behaviori,:);
  for channel = 1:3,
    handles.labels_plot.im(1,idx,channel) = curColor(channel);
  end
end

handles.labels_plot.predicted_im(:) = 0;
prediction= handles.data.GetPredictedIdx(handles.expi,handles.flies);
predictedidx = prediction.predictedidx;
scores = handles.data.NormalizeScores(prediction.scoresidx);
scores_old = handles.data.NormalizeScores(prediction.scoresidx_old);
latest = prediction.latest;

for behaviori = 1:handles.data.nbehaviors

  idxScores = (predictedidx == behaviori) & latest;
  idxPredict = idxScores & ...
    (abs(scores)>handles.data.GetConfidenceThreshold(behaviori));
  for channel = 1:3,
    if ~prediction.isValidated % Different color whether the scores are validated or not
      handles.labels_plot.predicted_im(1,idxPredict,channel) = handles.labelcolors(behaviori,channel);
      scoreNdx = ceil(scores(idxScores)*31)+32;
      handles.labels_plot.predicted_im(2,idxScores,channel) = handles.scorecolor(scoreNdx,channel,1);
    else
      idxLabeledPredict = idxPredict & labelidx.vals>0.5;
      idxUnlabeledPredict = idxPredict & labelidx.vals<0.5;
      idxLabeledScores = idxScores & labelidx.vals>0.5;
      idxUnlabeledScores = idxScores & labelidx.vals<0.5;

      shiftedColor = ShiftColor.shiftColorBkwd(handles.labelcolors(behaviori,:));
      handles.labels_plot.predicted_im(1,idxLabeledPredict,channel) = shiftedColor(channel);
      scoreNdx = ceil(scores(idxLabeledScores)*31)+32;
      handles.labels_plot.predicted_im(2,idxLabeledScores,channel) = handles.scorecolor(scoreNdx,channel,3);
    
      % Color doesn't change for unlabeled data.
      handles.labels_plot.predicted_im(1,idxUnlabeledPredict,channel) = handles.labelcolors(behaviori,channel);
      scoreNdx = ceil(scores(idxUnlabeledScores)*31)+32;
      handles.labels_plot.predicted_im(2,idxUnlabeledScores,channel) = handles.scorecolor(scoreNdx,channel,1);

    end
    
    % Old scores.
    oldScoreNdx = ceil(scores_old(idxScores)*31)+32;
    handles.labels_plot.predicted_im(3,idxScores,channel) = handles.scorecolor(oldScoreNdx,channel,1);
    
  end    
    
    
  % Old scores loaded from scores files are shown in different color
  idxScores = (predictedidx == behaviori) & (~latest);
  idxPredict = idxScores & ...
    (abs(scores)>handles.data.GetConfidenceThreshold(behaviori));
  shiftedColor = ShiftColor.shiftColorFwd(handles.labelcolors(behaviori,:));
  for channel = 1:3,
    handles.labels_plot.predicted_im(1,idxPredict,channel) = shiftedColor(channel);
    scoreNdx = ceil(scores(idxScores)*31)+32;
    handles.labels_plot.predicted_im(2,idxScores,channel) = handles.scorecolor(scoreNdx,channel,2);
    handles.labels_plot.predicted_im(3,idxScores,channel) = handles.scorecolor(scoreNdx,channel,2);
  end
  
end

[error_t0s,error_t1s] = get_interval_ends(labelidx.vals ~= 0 & predictedidx ~= 0 & ...
  labelidx.vals ~= predictedidx);
error_t0s = error_t0s + handles.t0_curr - 1.5;
error_t1s = error_t1s + handles.t0_curr - 1.5;
handles.labels_plot.error_xs = reshape([error_t0s;error_t1s;nan(size(error_t0s))],[1,numel(error_t0s)*3]);
set(handles.htimeline_errors,'XData',handles.labels_plot.error_xs,...
  'YData',zeros(size(handles.labels_plot.error_xs))+1.5);
[suggest_t0s,suggest_t1s] = get_interval_ends(labelidx.vals == 0 & predictedidx ~= 0);
suggest_t0s = suggest_t0s + handles.t0_curr - 1.5;
suggest_t1s = suggest_t1s + handles.t0_curr - 1.5;
handles.labels_plot.suggest_xs = reshape([suggest_t0s;suggest_t1s;nan(size(suggest_t0s))],[1,numel(suggest_t0s)*3]);
set(handles.htimeline_suggestions,'XData',handles.labels_plot.suggest_xs,...
  'YData',zeros(size(handles.labels_plot.suggest_xs))+1.5);

%{
%handles.labels_plot.suggested_im(:) = 0;
%for behaviori = 1:handles.data.nbehaviors
%  idx = handles.data.suggestedidx == behaviori;
%  for channel = 1:3,
%    handles.labels_plot.suggested_im(1,idx,channel) = handles.labelcolors(behaviori,channel);
%  end
%end
%handles.labels_plot.error_im(:) = 0;
%idx = handles.data.erroridx == 1;
%for channel = 1:3,
%  handles.labels_plot.error_im(1,idx,channel) = handles.correctcolor(channel);
%end
%idx = handles.data.erroridx == 2;
%for channel = 1:3,
%  handles.labels_plot.error_im(1,idx,channel) = handles.incorrectcolor(channel);
%end
%}

handles.labels_plot.isstart = ...
cat(2,labelidx.vals(1)~=0,...
labelidx.vals(2:end)~=0 & ...
labelidx.vals(1:end-1)~=labelidx.vals(2:end));


% set current frame
function handles = SetCurrentFrame(handles,i,t,hObject,doforce,doupdateplot)

if ~exist('doforce','var'),
  doforce = false;
end
if ~exist('doupdateplot','var'),
  doupdateplot = true;
end

t = round(t);

% check for change
if doforce || handles.ts(i) ~= t,

  handles.ts(i) = t;
  
  % update labels
%   if handles.label_state < 0,
%     handles = SetLabelPlot(handles,min(handles.t1_curr,max(handles.t0_curr,t)),0);
%   elseif handles.label_state > 0,
%     handles = SetLabelPlot(handles,min(handles.t1_curr,max(handles.t0_curr,t)),handles.label_state);
%   end
  
  % update slider
  if hObject ~= handles.slider_previews(i),
    set(handles.slider_previews(i),'Value',t);
  end

  % update frame number edit box
  if hObject ~= handles.edit_framenumbers(i),
    set(handles.edit_framenumbers(i),'String',num2str(t));
  end
  
  % update selection
  if handles.selecting,
    handles.selected_ts(end) = t;
    UpdateSelection(handles);
  end
  
  guidata(handles.figure_JLabel,handles);

  % update plot
  if doupdateplot,
    UpdatePlots(handles,'axes',i);
  end
  
  % TODO: update timeline zoom
  for h = handles.axes_timeline_labels,
    zoom(h,'reset');
  end
  
%   % out of bounds for labeling? then turn off labeling
%   if (t < handles.t0_curr || t > handles.t1_curr),
%     if handles.label_state > 0,
%       set(handles.togglebutton_label_behaviors(handles.label_state),'Value',0);
%     elseif handles.label_state < 0,
%       set(handles.togglebutton_label_unknown,'Value',0);
%     end
%     handles.label_state = 0;
%     set([handles.togglebutton_label_behaviors,handles.togglebutton_label_unknown],'Enable','off');
%   else
%     set([handles.togglebutton_label_behaviors,handles.togglebutton_label_unknown],'Enable','on');
%   end

  
end

% --- Executes during object creation, after setting all properties.
function slider_preview_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to slider_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function pushtool_save_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to pushtool_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_editfiles_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_editfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles,'data'),
%   [success,msg] = handles.data.UpdateStatusTable();
%   if ~success,
%     error(msg);
%   end
  params = {handles.data};
else
  params = {handles.configfilename};
  if isfield(handles,'defaultpath'),
    params(end+1:end+2) = {'defaultpath',handles.defaultpath};
  end
end
if ~isempty(handles.expi) && handles.expi > 0,
  oldexpdir = handles.data.expdirs{handles.expi};
else
  oldexpdir = '';
end
[handles.data,success] = JLabelEditFiles(params{:});
handles.data.SetStatusFn(@(s) SetStatusCallback(s,handles.figure_JLabel));
handles.data.SetClearStatusFn(@() ClearStatusCallback(handles.figure_JLabel));

handles.defaultpath = handles.data.defaultpath;
if ~success,
  guidata(hObject,handles);
  return;
end
  
%   % rearrange labels
%   newlabels = struct('t0s',cell(1,numel(expdirs)),...
%     't1s',cell(1,numel(expdirs)),...
%     'names',cell(1,numel(expdirs)),...
%     'flies',cell(1,numel(expdirs)));
% 
%   % loop through new experiments
%   for i = 1:numel(expdirs),
%     
%     % is this experiment already loaded in?
%     j = find(strcmp(expdirs{i},handles.expdirs),1);
%     
%     % if this is a new experiment, try to load the labels from file
%     if isempty(j),
%       labelmatname = fullfile(expdirs{i},handles.configparams.file.labelfilename);
%       if exist(labelmatname,'file'),
%         try
%           loadedlabels = load(labelmatname,'t0s','t1s','names','flies');
%           newlabels(i).t0s = loadedlabels.t0s;
%           newlabels(i).t1s = loadedlabels.t1s;
%           newlabels(i).names = loadedlabels.names;
%           newlabels(i).flies = loadedlabels.flies;
%           newwindowdata_labeled = [];
%           newwindowdata_labeled_movie = [];
%           newwindowdata_labeled_flies = [];
%         catch ME,
%           uiwait(warndlg('Error loading labels from %s: %s',labelmatname,getReport(ME)));
%         end
%         didloadtrxfile = false;
%         for k = 1:numel(newlabels(i).t0s),
%           if isempty(newlabels(i).t0s{k}),
%             continue;
%           end
%           flies = newlabels(i).flies{k};
%           % TODO: multiple flies
%           windowfilename = GetWindowFileName(handles,outexpdirs{i},flies(1));
%           windowdata = load(windowfilename);
%           islabeled = false(1,size(windowdata.X,2));
% 
%           if ~didloadtrxfile,
%             % TODO: remove this
%             trxfilename = fullfile(expdirs{i},handles.configparams.file.trxfilename);
%             global CACHED_TRX; %#ok<TLEV>
%             if isempty(CACHED_TRX),
%               trx_curr = load_tracks(trxfilename);
%               CACHED_TRX = trx_curr;
%             else
%               trx_curr = CACHED_TRX;
%             end
%             didloadtrxfile = true;
%           end          
%           off = trx_curr(flies(1)).off;
%           
%           for l = 1:numel(newlabels(i).t0s{k}),
%             i0 = newlabels(i).t0s{k}(l)+off;
%             i1 = newlabels(i).t1s{k}(l)+off;
%             islabeled(i0:i1) = true;
%           end
%           if ~isempty(newwindowdata_labeled),
%             nnew = size(windowdata.X,1);
%             nold = size(newwindowdata_labeled,1);
%             if nnew > nold,
%               newwindowdata_labeled(end+1:end+nnew-nold,:) = nan;
%             elseif nold > nnew,
%               windowdata.X(end+1:end+nold-nnew) = nan;
%             end
%           end
%           n = nnz(islabeled);
%           newwindowdata_labeled(:,end+1:end+n) = windowdata.X(:,islabeled);
%           newwindowdata_labeled_movie(end+1:end+n) = i;
%           newwindowdata_labeled_flies(end+1:end+n,:) = repmat(flies,[n,1]);
%         end
%       end
%     else
%       % existing experiment, just copy labels
%       newlabels(i) = handles.labels(j);
%       idx = handles.windowdata_labeled_movie == i;
%       n = nnz(idx);
%       newwindowdata_labeled(:,end+1:end+n) = handles.windowdata_labeled(:,idx);
%       newwindowdata_labeled_movie(end+1:end+n) = i;
%       newwindowdata_labeled_flies(end+1:end+n,:) = handles.windowdata_labeled_flies(idx,:);      
%     end
%   end
%   handles.labels = newlabels;
    
% save needed if list has changed
handles = SetNeedSave(handles);

if ~isempty(oldexpdir) && ismember(oldexpdir,handles.data.expdirs),
  j = find(strcmp(oldexpdir,handles.data.expdirs),1);
  handles.expi = j;
else
  handles = UnsetCurrentMovie(handles);
  if handles.data.nexps > 0 && handles.data.expi == 0,
    handles = SetCurrentMovie(handles,1);
  else
    handles = SetCurrentMovie(handles,handles.data.expi);
  end
end
  
guidata(hObject,handles);
  
function handles = SetNeedSave(handles)

handles.needsave = true;
set(handles.menu_file_save,'Enable','on');
set(handles.menu_file_save_labels,'Enable','on');

function handles = SetSaved(handles)

handles.needsave = false;
set(handles.menu_file_save,'Enable','off');
set(handles.menu_file_save_labels,'Enable','off');

% --------------------------------------------------------------------
function success = menu_file_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uiputfile('*.mat','Save classifier',handles.data.classifierfilename);
if ~ischar(filename),
  success = false;
  return;
end
handles.data.classifierfilename = fullfile(pathname,filename);
SetStatus(handles,sprintf('Saving classifier to %s',handles.data.classifierfilename));
handles.data.SaveLabels();
handles.data.SaveClassifier();
handles = SetSaved(handles);
ClearStatus(handles);
success = true;

% --------------------------------------------------------------------
function menu_file_exit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_JLabel_CloseRequestFcn(hObject, eventdata, handles);

% --------------------------------------------------------------------
function menu_edit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_view_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_go_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function menu_edit_undo_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_undo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function [handles,success] = LoadConfig(handles,forceui)

if ~exist('forceui','var'),
  forceui = false;
end

% initialize output to success = false
success = false;

% get config file
havefilename = false;
if ~forceui && ~isempty(handles.configfilename),
  havefilename = true;
end

% default config file name
if isempty(handles.configfilename),
  defaultconfigfilename = fullfile(handles.defaultpath,'JLabelConfig.xml');
else
  defaultpath = myfileparts(handles.configfilename);
  if exist(defaultpath,'file'),
    handles.defaultpath = defaultpath;
  end
  defaultconfigfilename = handles.configfilename;
end

% loop until we have a valid config file
while true,

  if ~havefilename,
    
    if ~exist(defaultconfigfilename,'file'),
      defaultconfigfilename = handles.defaultpath;
    end
  
    % get from user
    [filename,pathname] = uigetfile('*.xml','Choose XML config file',defaultconfigfilename);
    
    % if cancel clicked, just return
    if ~ischar(filename),
      return;
    end
    
    handles.configfilename = fullfile(pathname,filename);

    % store path as default
    handles.defaultpath = pathname;

  end
  
  % make sure the file exists
  if ~exist(handles.configfilename,'file'),
    havefilename = false;
    uiwait(warndlg(sprintf('File %s does not exist',handles.configfilename),'Error reading config file'));
    continue;
  end
    
%   try
    handles.configparams = ReadXMLParams(handles.configfilename);
    if ~isfield(handles.configparams,'JLabelMode')
      handles.configparams.JLabelMode.mode = 'basic';
    end
%   catch ME,
%     uiwait(warndlg(sprintf('Error reading configuration from file %s: %s',handles.configfilename,getReport(ME)),'Error reading config file'));
%     havefilename = false;
%     continue;
%   end
  
  % success -- break
  break;

end

success = true;

function handles = InitializeState(handles)

handles = LoadRC(handles);

% whether save is necessary
handles.needsave = false;

% initialize data structure
handles.data = JLabelData(handles.configfilename,...
  'defaultpath',handles.defaultpath,...
  'classifierfilename',handles.classifierfilename,...
  'setstatusfn',@(s) SetStatusCallback(s,handles.figure_JLabel),...
  'clearstatusfn',@() ClearStatusCallback(handles.figure_JLabel));

% number of flies to label at a time
handles.nflies_label = 1;

% learned classifier
handles.classifier = [];

% currently shown experiment
handles.expi = 0;
% currently labeled flies
handles.flies = 1:handles.nflies_label;
% currently shown frame
handles.ts = 0;

% current behavior labeling state: nothing down
handles.label_state = 0;
handles.label_imp = [];

% number of flies for the current movie
handles.nflies_curr = 0;

% label colors
if isfield(handles.configparams,'behaviors') && ...
    isfield(handles.configparams.behaviors,'labelcolors'),
  labelcolors = handles.configparams.behaviors.labelcolors;
  if numel(labelcolors) >= 3*handles.data.nbehaviors,
    handles.labelcolors = reshape(labelcolors(1:3*handles.data.nbehaviors),[handles.data.nbehaviors,3]);
  else
    uiwait(warndlg('Error parsing label colors from config file, automatically assigning','Error parsing config label colors'));
    if isfield(handles.configparams,'labels') && ...
        isfield(handles.configparams.labels,'colormap'),
      cm = handles.configparams.labels.colormap;
    else
      cm = 'lines';
    end
    if ~exist(cm,'file'),
      cm = 'lines';
    end
%     try
      handles.labelcolors = eval(sprintf('%s(%d)',cm,handles.data.nbehaviors));
%     catch ME,
%       uiwait(warndlg(sprintf('Error using label colormap from config file: %s',getReport(ME)),'Error parsing config label colors'));
%       handles.labelcolors = lines(handles.data.nbehaviors);
%     end
  end
end
handles.labelunknowncolor = [0,0,0];
if isfield(handles.configparams,'behaviors') && ...
    isfield(handles.configparams.behaviors,'unknowncolor'),
  unknowncolor = handles.configparams.behaviors.unknowncolor;
  if numel(unknowncolor) >= 3,
    handles.labelunknowncolor = reshape(unknowncolor(1:3),[1,3]);
  else
    uiwait(warndlg('Error parsing unknown color from config file, automatically assigning','Error parsing config unknown colors'));
  end
end
handles.flies_extra_markersize = 12;
if isfield(handles.configparams,'plot') && ...
    isfield(handles.configparams.plot,'trx') && ...
    isfield(handles.configparams.plot.trx,'extra_markersize'),
  handles.flies_extra_markersize = handles.configparams.plot.trx.extra_markersize(1);
end

for channel = 1:3
  midValue = handles.labelunknowncolor(channel);
  startValue = handles.labelcolors(2,channel);
  endValue = handles.labelcolors(1,channel);
  handles.scorecolor(1:32,channel,1) = (midValue-startValue)*(0:31)/31+startValue;
  handles.scorecolor(32:63,channel,1) = (endValue-midValue)*(0:31)/31+midValue;
end
for ndx = 1:63
  handles.scorecolor(ndx,:,2) = ShiftColor.shiftColorFwd(handles.scorecolor(ndx,:,1));
  handles.scorecolor(ndx,:,3) = ShiftColor.shiftColorBkwd(handles.scorecolor(ndx,:,1));  
end

handles.correctcolor = [0,.7,0];
handles.incorrectcolor = [.7,.7,0];
handles.suggestcolor = [0,.7,.7];

handles.selection_color = [1,.6,0];
handles.selection_alpha = .5;

% color for showing which labels are being plotted
handles.emphasiscolor = [.7,.7,0];
handles.unemphasiscolor = [1,1,1];

% create buttons for each label
handles = CreateLabelButtons(handles);

% timeline properties
handles.timeline_prop_remove_string = '<html><body><i>Remove</i></body></html>';
handles.timeline_prop_help_string = '<html><body><i>Help</i></body></html>';
handles.timeline_prop_options = ...
  {handles.timeline_prop_remove_string,...
  handles.timeline_prop_help_string};

if ~isempty(handles.data.allperframefns)
  for i = 1:numel(handles.data.allperframefns),
    handles.timeline_prop_options{end+1} = handles.data.allperframefns{i};
  end
  handles.d = handles.data.allperframefns(1);
  handles.perframepropis = 1;
  set(handles.timeline_label_prop1,'String',handles.timeline_prop_options,'Value',3);
  handles.timeline_data_ylims = nan(2,numel(handles.data.allperframefns));
end

% maximum distance squared in fraction of axis to change frames when
% clicking on preview window
handles.max_click_dist_preview = .025^2;

% zoom state
handles.preview_zoom_mode = 'center_on_fly';
handles.zoom_fly_radius = nan(1,2);
handles.menu_view_zoom_options = findall(handles.menu_view_zoom,'Type','menu');
set(handles.menu_view_zoom_options,'Checked','off');
set(handles.menu_view_zoom_in_on_fly,'Checked','on');

% last clicked object
handles.selection_t0 = nan;
handles.selection_t1 = nan;
handles.selected_ts = nan(1,2);
handles.buttondown_t0 = nan;
handles.buttondown_axes = nan;
set([handles.pushbutton_playselection,handles.pushbutton_clearselection],'Enable','off');

% not selecting
handles.selecting = false;
set(handles.togglebutton_select,'Value',0);

% initialize nextjump obj;
handles.NJObj = NextJump();
handles.NJObj.SetSeekBehaviorsGo(1:handles.data.nbehaviors);
handles.NJObj.SetPerframefns(handles.data.allperframefns);
if isfield(handles.rc,'navPreferences')
  handles.NJObj.SetState(handles.rc.navPreferences);
end

% initialize labels for navigation
SetJumpGoMenuLabels(handles)

% label shortcuts
if numel(handles.label_shortcuts) ~= 2*handles.data.nbehaviors + 1,
  handles.label_shortcuts = cellstr(num2str((1:2*handles.data.nbehaviors+1)'))';
end

% play/stop
handles.hplaying = nan;
%handles.play_FPS = 2;

%handles.traj_nprev = 25;
%handles.traj_npost = 25;

% whether to show trajectories
set(handles.menu_view_plottracks,'Checked','on');

% bookmarked clips windows
handles.bookmark_windows = [];

% whether to plot manual labels or automatic labels
handles.plot_labels_manual = true;
handles.plot_labels_automatic = false;
set(handles.menu_view_plot_labels_manual,'Checked','on');
set(handles.menu_view_plot_labels_automatic,'Checked','off');

buttonNames = {'pushbutton_train','pushbutton_predict',...
              'togglebutton_select','pushbutton_clearselection',...
              'pushbutton_playselection','pushbutton_playstop',...
              'similarFramesButton','bagButton'};
  
for buttonNum = 1:numel(buttonNames)
  SetButtonImage(handles.(buttonNames{buttonNum}));
end


function SetJumpGoMenuLabels(handles)

set(handles.menu_go_forward_X_frames,'Label',sprintf('Forward %d frames (down arrow)',handles.nframes_jump_go));
set(handles.menu_go_back_X_frames,'Label',sprintf('Back %d frames (up arrow)',handles.nframes_jump_go));
jumpType = handles.NJObj.GetCurrentType();
set(handles.menu_go_next_automatic_bout_start,'Label',...
  sprintf('Next %s bout start (shift + right arrow)',jumpType));
set(handles.menu_go_previous_automatic_bout_end,'Label',...
  sprintf('Next %s bout end (shift + left arrow)',jumpType));

% create buttons for each label
function handles = CreateLabelButtons(handles)

% get positions of stuff
set(handles.panel_labelbuttons,'Units','pixels');
panel_pos = get(handles.panel_labelbuttons,'Position');
select_pos = get(handles.panel_select,'Position');
set(handles.togglebutton_label_behavior1,'Units','pixels');
button1_pos = get(handles.togglebutton_label_behavior1,'Position');
set(handles.togglebutton_label_unknown,'Units','pixels');
unknown_button_pos = get(handles.togglebutton_label_unknown,'Position');
out_border_y = unknown_button_pos(2);
out_border_x = unknown_button_pos(1);
in_border_y = button1_pos(2) - (unknown_button_pos(2)+unknown_button_pos(4));
button_width = button1_pos(3);
button_height = button1_pos(4);

% calculate new height for the panel

if strcmp(handles.configparams.JLabelMode.mode,'basic');
new_panel_height = 2*out_border_y + (handles.data.nbehaviors+1)*button_height + ...
  handles.data.nbehaviors*in_border_y;
else
new_panel_height = 2*out_border_y + (2*handles.data.nbehaviors+1)*button_height + ...
  2*handles.data.nbehaviors*in_border_y;
end
% update panel position
panel_top = panel_pos(2)+panel_pos(4);
new_panel_pos = [panel_pos(1),panel_top-new_panel_height,panel_pos(3),new_panel_height];
set(handles.panel_labelbuttons,'Position',new_panel_pos);
dy_label_select = panel_pos(2) - select_pos(2) - select_pos(4);
new_select_pos = [select_pos(1),new_panel_pos(2)-select_pos(4)-dy_label_select,select_pos(3:4)];
set(handles.panel_select,'Position',new_select_pos);

% move unknown button to the bottom
new_unknown_button_pos = [unknown_button_pos(1),out_border_y,unknown_button_pos(3),button_height];
set(handles.togglebutton_label_unknown,'Position',new_unknown_button_pos);

% list of buttons
handles.togglebutton_label_behaviors = nan(1,2*handles.data.nbehaviors);

% update first button
if ~strcmp(handles.configparams.JLabelMode.mode,'basic')
  new_button1_pos = [out_border_x,new_panel_height-out_border_y-button_height,button_width,button_height];
  set(handles.togglebutton_label_behavior1,...
    'String',sprintf('Important %s',handles.data.labelnames{1}),...
    'ForegroundColor','w','Units','pixels','FontUnits','pixels','FontSize',14,...
    'FontWeight','bold','BackgroundColor',ShiftColor.increaseIntensity(handles.labelcolors(1,:)),...
    'Position',new_button1_pos,...
    'UserData',1);
  handles.togglebutton_label_behaviors(1) = handles.togglebutton_label_behavior1;
  SetButtonImage(handles.togglebutton_label_behavior1);
  pos = [out_border_x,new_panel_height-out_border_y-2*button_height-in_border_y,button_width,button_height];
  handles.togglebutton_label_behaviors(2) = uicontrol('Style','togglebutton',...
    'String',sprintf('%s',handles.data.labelnames{1}),...
    'ForegroundColor','w','Units','pixels','FontUnits','pixels','FontSize',14,...
    'FontWeight','bold','BackgroundColor',handles.labelcolors(1,:),...
    'Parent',handles.panel_labelbuttons,...
    'Callback',get(handles.togglebutton_label_behavior1,'Callback'),...
    'Position',pos,'Tag',sprintf('togglebutton_label_normbehavior1'),...
    'UserData',2);
  SetButtonImage(handles.togglebutton_label_behaviors(2));
else
  pos = [out_border_x,new_panel_height-out_border_y-button_height,button_width,button_height];
  set(handles.togglebutton_label_behavior1,...
    'String',sprintf('%s',handles.data.labelnames{1}),...
    'ForegroundColor','w','Units','pixels','FontUnits','pixels','FontSize',14,...
    'FontWeight','bold','BackgroundColor',handles.labelcolors(1,:),...
    'Position',pos,...
    'UserData',2);
  handles.togglebutton_label_behaviors(2) = handles.togglebutton_label_behavior1;
end

  

% create the rest of the buttons
for i = 2:handles.data.nbehaviors,
  if ~strcmp(handles.configparams.JLabelMode.mode,'basic')
    pos = [out_border_x,new_panel_height-out_border_y-button_height*(2*i-1)-in_border_y*(2*i-2),...
      button_width,button_height];
    handles.togglebutton_label_behaviors(2*i-1) = ...
      uicontrol('Style','togglebutton','String',sprintf('Important %s',handles.data.labelnames{i}),...
      'ForegroundColor','w','Units','pixels','FontUnits','pixels','FontSize',14,...
      'FontWeight','bold','BackgroundColor',ShiftColor.increaseIntensity(handles.labelcolors(i,:)),...
      'Position',pos,...
      'Callback',get(handles.togglebutton_label_behavior1,'Callback'),...
      'Parent',handles.panel_labelbuttons,...
      'Tag',sprintf('togglebutton_label_behavior%d',i),...
      'UserData',2*i-1);
    SetButtonImage(handles.togglebutton_label_behaviors(2*i-1));
  pos = [out_border_x,new_panel_height-out_border_y-button_height*(2*i)-in_border_y*(2*i-1),...
    button_width,button_height];
  else
  pos = [out_border_x,new_panel_height-out_border_y-button_height*i-in_border_y*(i-1),...
    button_width,button_height];    
  end
  handles.togglebutton_label_behaviors(2*i) = ...
    uicontrol('Style','togglebutton','String',sprintf('%s',handles.data.labelnames{i}),...
    'ForegroundColor','w','Units','pixels','FontUnits','pixels','FontSize',14,...
    'FontWeight','bold','BackgroundColor',handles.labelcolors(i,:),...
    'Position',pos,...
    'Callback',get(handles.togglebutton_label_behavior1,'Callback'),...
    'Parent',handles.panel_labelbuttons,...
    'Tag',sprintf('togglebutton_label_normbehavior%d',i),...
    'UserData',2*i);
  SetButtonImage(handles.togglebutton_label_behaviors(2*i));
end

% set props for unknown button
set(handles.togglebutton_label_unknown,...
  'String','Label Unknown',...
  'ForegroundColor','w','Units','pixels','FontUnits','pixels','FontSize',14,...
  'FontWeight','bold','BackgroundColor',handles.labelunknowncolor,...
  'UserData',-1);
SetButtonImage(handles.togglebutton_label_unknown);

  
function EnableGUI(handles)

% these controls require a movie to currently be open
if strcmp(handles.configparams.JLabelMode.mode,'basic')
  h = [handles.contextmenu_timeline_manual_timeline_options,...
    handles.togglebutton_label_behaviors(2:2:end),...
    handles.togglebutton_label_unknown,...
    handles.menu_view_zoom_options(:)'];
else
  h = [handles.contextmenu_timeline_manual_timeline_options,...
    handles.togglebutton_label_behaviors(:)',...
    handles.togglebutton_label_unknown,...
    handles.menu_view_zoom_options(:)'];
end
hp = [handles.panel_previews(:)',...
  handles.panel_timelines,...
  handles.panel_learn];
if handles.expi >= 1 && handles.expi <= handles.data.nexps,
  set(h,'Enable','on');
  set(hp,'Visible','on');
else
  set(h,'Enable','off');
  set(hp,'Visible','off');
end

% whether we need to save
if handles.needsave,
  set(handles.menu_file_save,'Enable','on');
else
  set(handles.menu_file_save,'Enable','off');
end

function handles = LoadRC(handles)

% rc file name
handles.rcfilename = fullfile(myfileparts(which('JLabel')),'.JLabelrc.mat');
handles.rc = struct;
if exist(handles.rcfilename,'file'),
%   try
    handles.rc = load(handles.rcfilename);
%   catch ME,
%     warning('Error loading rc file %s: %s',handles.rcfilename,getReport(ME));
%   end
end
% try
  if isfield(handles.rc,'defaultpath'),
    handles.defaultpath = handles.rc.defaultpath;
    if isfield(handles,'data'),
      handles.data.SetDefaultPath(handles.defaultpath);
    end
  end
  if isfield(handles.rc,'figure_JLabel_Position_px'),
    pos = handles.rc.figure_JLabel_Position_px;
    set(handles.figure_JLabel,'Units','pixels');
    % TODO: remove this once resizing is implemented
    pos0 = get(handles.figure_JLabel,'Position');
    pos(3:4) = pos0(3:4);
    set(handles.figure_JLabel,'Position',pos);
  end
  if isfield(handles.rc,'timeline_nframes'),
    handles.timeline_nframes = handles.rc.timeline_nframes;
  else
    handles.timeline_nframes = 250;
  end
  if isfield(handles.rc,'nframes_jump_go')
    handles.nframes_jump_go = handles.rc.nframes_jump_go;
  else
    handles.nframes_jump_go = 30;
  end
  
  if isfield(handles.rc,'label_shortcuts'),
    handles.label_shortcuts = handles.rc.label_shortcuts;
  else
    handles.label_shortcuts = [];
  end

  %output avi options
  
  % compression: scheme for compression for output avis
  if isfield(handles.rc,'outavi_compression'),
    handles.outavi_compression = handles.rc.outavi_compression;
  else
    handles.outavi_compression = 'None';
  end
  % outavi_fps: output frames per second
  if isfield(handles.rc,'outavi_fps'),
    handles.outavi_fps = handles.rc.outavi_fps;
  else
    handles.outavi_fps = 15;
  end
  % useVideoWriter: whether to use videowriter class
  if isfield(handles.rc,'useVideoWriter'),
    handles.useVideoWriter = handles.rc.useVideoWriter > 0;
  else
    handles.useVideoWriter = exist('VideoWriter','file') > 0;
  end
  
  % preview options
  
  % playback speed
  if isfield(handles.rc,'play_FPS'),
    handles.play_FPS = handles.rc.play_FPS;
  else
    handles.play_FPS = 2;
  end
  
  if isfield(handles.rc,'traj_nprev'),
    handles.traj_nprev = handles.rc.traj_nprev;
  else
    handles.traj_nprev = 25;
  end
  
  if isfield(handles.rc,'traj_npost'),
    handles.traj_npost = handles.rc.traj_npost;
  else
    handles.traj_npost = 25;
  end
  
% catch ME,
%   warning('Error loading RC file: %s',getReport(ME));  
% end

function handles = SaveRC(handles)

% try
  if ~isfield(handles,'rcfilename'),
    handles.rcfilename = fullfile(myfileparts(which('JLabel')),'.JLabelrc.mat');
  end
  
  if isfield(handles,'rc'),
    rc = handles.rc;
  else
    rc = struct;
  end
  
  if isfield(handles,'data'),
    rc.defaultpath = handles.data.defaultpath;
  elseif isfield(handles,'defaultpath'),
    rc.defaultpath = handles.defaultpath;
  end
  if isfield(handles,'timeline_nframes'),
    rc.timeline_nframes = handles.timeline_nframes;
  end
  
  set(handles.figure_JLabel,'Units','pixels');
  rc.figure_JLabel_Position_px = get(handles.figure_JLabel,'Position');
  
  if isfield(handles,'nframes_jump_go'),
    rc.nframes_jump_go = handles.nframes_jump_go;
  end
  
  % label shortcuts
  if isfield(handles,'label_shortcuts'),
    rc.label_shortcuts = handles.label_shortcuts;
  end
  
  
  %output avi options
  
  % compression: scheme for compression for output avis
  if isfield(handles,'outavi_compression'),
    rc.outavi_compression = handles.outavi_compression;
  end
  % outavi_fps: output frames per second
  if isfield(handles,'outavi_fps'),
    rc.outavi_fps = handles.outavi_fps;
  end
  % useVideoWriter: whether to use videowriter class
  if isfield(handles,'useVideoWriter'),
    rc.useVideoWriter = handles.useVideoWriter;
  end
  
  % preview options
  
  % playback speed
  if isfield(handles,'play_FPS'),
    rc.play_FPS = handles.play_FPS;
  end
  
  if isfield(handles,'traj_nprev'),
    rc.traj_nprev = handles.traj_nprev;
  end
  
  if isfield(handles,'traj_npost'),
    rc.traj_npost = handles.traj_npost;
  end  
  
  % navigation preferences
  rc.navPreferences = handles.NJObj.GetState();
  
  save(handles.rcfilename,'-struct','rc');

% catch ME,
%   warning('Error saving RC file: %s',getReport(ME));
% end


% --- Executes when user attempts to close figure_JLabel.
function figure_JLabel_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_JLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%delete(hObject);

% check if we need to save
if handles.needsave,
  res = questdlg('Save before quitting?','Save?','Yes','No','Cancel','Yes');
  if strcmpi(res,'Yes'),
    success = menu_file_save_Callback(hObject, eventdata, handles);
    if ~success,
      return;
    end
    handles = guidata(hObject);
  elseif strcmpi(res,'Cancel'),
    return;
  end
end

if isfield(handles,'movie_fid') && ~isempty(handles.movie_fid) && ...
    handles.movie_fid > 1 && ~isempty(fopen(handles.movie_fid)),
  fclose(handles.movie_fid);
  handles.movie_fid = [];
end
% try
  % turn off zooming
  zoom(handles.figure_JLabel,'off');
% catch %#ok<CTCH>
% end
% SWITCH THIS
if true,
  SaveRC(handles);
  delete(handles.figure_JLabel);
else
  uiresume(handles.figure_JLabel); %#ok<UNRCH>
end


% --------------------------------------------------------------------
function toggletool_zoomin_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toggletool_zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function toggletool_zoomin_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggletool_zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set([handles.toggletool_zoomout,handles.toggletool_pan],'State','off');
pan(handles.figure_JLabel,'off');
set(handles.hzoom,'Direction','in','Enable','on');


% --------------------------------------------------------------------
function toggletool_zoomin_OffCallback(hObject, eventdata, handles)
% hObject    handle to toggletool_zoomin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.hzoom,'Enable','off');


% --------------------------------------------------------------------
function toggletool_zoomout_OffCallback(hObject, eventdata, handles)
% hObject    handle to toggletool_zoomout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.hzoom,'Enable','off');

% --------------------------------------------------------------------
function toggletool_zoomout_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggletool_zoomout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set([handles.toggletool_zoomin,handles.toggletool_pan],'State','off');
pan(handles.figure_JLabel,'off');
set(handles.hzoom,'Direction','out','Enable','on');


% --------------------------------------------------------------------
function toggletool_zoomout_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to toggletool_zoomout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function toggletool_pan_OnCallback(hObject, eventdata, handles)
% hObject    handle to toggletool_pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set([handles.toggletool_zoomin,handles.toggletool_zoomout],'State','off');
zoom(handles.figure_JLabel,'off');
pan(handles.figure_JLabel,'on');

% --------------------------------------------------------------------
function toggletool_pan_OffCallback(hObject, eventdata, handles)
% hObject    handle to toggletool_pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pan(handles.figure_JLabel,'off');


% --- Executes on button press in togglebutton_label_behavior1.
function togglebutton_label_behavior1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_label_behavior1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_label_behavior1
buttonNum = get(hObject,'UserData');
behaviori = ceil(buttonNum/2);
isImportant = mod(buttonNum,2);

if get(hObject,'Value'),
  % toggle on, label pen is down.
  
  handles.label_state = behaviori;
  handles.label_imp = isImportant;
  handles.label_t0 = handles.ts(1);
  
  % set everything else to off
  for j = 1:2*handles.data.nbehaviors,
    if j == buttonNum || isnan(handles.togglebutton_label_behaviors(j)),
      continue;
    end
    set(handles.togglebutton_label_behaviors(j),'Value',0,'Enable','off');
  end
  set(handles.togglebutton_label_unknown,'Value',0,'Enable','off');

  curColor = handles.labelcolors(behaviori,:);
  if isImportant, curColor = ShiftColor.increaseIntensity(curColor); end
  set(handles.htimeline_label_curr,'XData',handles.label_t0 + [-.5,-.5,.5,.5,-.5],...
    'FaceColor',curColor);
  % set the current frame to be labeled
  %handles.lastframe_labeled = [];
  %handles = SetLabelPlot(handles,min(handles.t1_curr,max(handles.t0_curr,handles.ts(1))),behaviori);
  
  UpdatePlots(handles,...
     'refreshim',false,'refreshflies',true,'refreshtrx',false,'refreshlabels',true,...
     'refresh_timeline_manual',true,...
     'refresh_timeline_auto',false,...
     'refresh_timeline_suggest',false,...
     'refresh_timeline_error',true,...
     'refresh_timeline_xlim',false,...
     'refresh_timeline_hcurr',false,...
     'refresh_timeline_props',false,...
     'refresh_timeline_selection',false,...
     'refresh_curr_prop',false);

   set(handles.menu_file,'enable','off');
   set(handles.menu_edit,'enable','off');
   set(handles.menu_go,'enable','off');
   set(handles.menu_classifier,'enable','off');
   
else % label pen is up.
  
  
  if handles.ts(1) <= handles.label_t0,
    t0 = handles.ts(1);
    t1 = handles.label_t0;
  else
    t0 = handles.label_t0;
    t1 = handles.ts(1);
  end
  t0 = min(handles.t1_curr,max(handles.t0_curr,t0));
  t1 = min(handles.t1_curr,max(handles.t0_curr,t1));
  handles = SetLabelPlot(handles,t0,t1,handles.label_state,handles.label_imp);

  handles.label_state = 0;
  handles.label_imp = [];
  handles.label_t0 = [];
  set(handles.htimeline_label_curr,'XData',nan(1,5));
  UpdatePlots(handles,...
    'refreshim',false,'refreshflies',true,'refreshtrx',false,'refreshlabels',true,...
    'refresh_timeline_manual',true,...
    'refresh_timeline_auto',false,...
    'refresh_timeline_suggest',false,...
    'refresh_timeline_error',true,...
    'refresh_timeline_xlim',false,...
    'refresh_timeline_hcurr',false,...
    'refresh_timeline_props',false,...
    'refresh_timeline_selection',false,...
    'refresh_curr_prop',false);
  
%   handles.data.StoreLabels();
  for j = 1:2*handles.data.nbehaviors,
    if isnan(handles.togglebutton_label_behaviors(j)), continue; end
    set(handles.togglebutton_label_behaviors(j),'Value',0,'Enable','on');
  end
  set(handles.togglebutton_label_unknown,'Value',0,'Enable','on');
  %set(handles.togglebutton_label_behaviors(behaviori),'String',sprintf('Label %s',handles.data.labelnames{behaviori}));

   set(handles.menu_file,'enable','on');
   set(handles.menu_edit,'enable','on');
   set(handles.menu_go,'enable','on');
   set(handles.menu_classifier,'enable','on');

end

guidata(hObject,handles);

function handles = SetLabelPlot(handles,t0,t1,behaviori,important)

% if behaviori == 0,
%   return;
% end

% if t == handles.lastframe_labeled,
%   warning('This should never happen');
%   keyboard;
% end
% 
% if isempty(handles.lastframe_labeled),
%   t0 = t;
%   t1 = t;
%   t2 = min(t+1,handles.t1_curr);
% else
%   if t < handles.lastframe_labeled,
%     t0 = t;
%     t1 = handles.lastframe_labeled-1;
%     t2 = handles.lastframe_labeled;
%   elseif t > handles.lastframe_labeled,
%     t0 = handles.lastframe_labeled+1;
%     t1 = t;
%     t2 = min(t+1,handles.t1_curr);
%   end
% end

% WARNING: this function directly accesses handles.data.labelidx, trx make sure
% that we've preloaded the right experiment and flies. 
% REMOVED!
% if handles.expi ~= handles.data.expi || ~all(handles.flies == handles.data.flies),
%   handles.data.Preload(handles.expi,handles.flies);
% end

if t1 < t0,
  tmp = t1;
  t1 = t0;
  t0 = tmp;
end
handles.labels_plot.x(:,t0+handles.labels_plot_off:t1+handles.labels_plot_off,:,:) = nan;
handles.labels_plot.y(:,t0+handles.labels_plot_off:t1+handles.labels_plot_off,:,:) = nan;
% handles.data.labelidx(t0+handles.labels_plot_off:t1+handles.labels_plot_off) = 0;

for channel = 1:3,
  handles.labels_plot.im(1,t0+handles.labels_plot_off:t1+handles.labels_plot_off,channel) = handles.labelunknowncolor(channel);
end
handles.data.SetLabel(handles.expi,handles.flies,t0:t1,behaviori,important);
if behaviori > 0,
  % handles.data.labelidx(t0+handles.labels_plot_off:t1+handles.labels_plot_off) = behaviori;
  for channel = 1:3,
    handles.labels_plot.im(1,t0+handles.labels_plot_off:t1+handles.labels_plot_off,channel) = handles.labelcolors(behaviori,channel);
  end
  for l = 1:numel(handles.flies),
    %off = handles.data.trx(handles.flies(l)).off;
    %j0 = t0+off;
    %j2 = t2+off;
    k0 = t0+handles.labels_plot_off;
    k2 = t1+handles.labels_plot_off+1;
    xplot = handles.data.GetTrxValues('X1',handles.expi,handles.flies(l),min(t0:t1+1,handles.t1_curr));
    yplot = handles.data.GetTrxValues('Y1',handles.expi,handles.flies(l),min(t0:t1+1,handles.t1_curr));
    handles.labels_plot.x(:,k0:k2-1,behaviori,l) = [xplot(1:end-1);xplot(2:end)];
    handles.labels_plot.y(:,k0:k2-1,behaviori,l) = [yplot(1:end-1);yplot(2:end)];      

%     handles.labels_plot.x(k0:k2,behaviori,l) = ...
%       handles.data.trx(handles.flies(l)).x(j0:j2);
%     handles.labels_plot.y(k0:k2,behaviori,l) = ...
%       handles.data.trx(handles.flies(l)).y(j0:j2);
  end
end

% isstart
if t0 == handles.t0_curr,
  handles.labels_plot.isstart(t0+handles.labels_plot_off) = behaviori ~= 0;
end
t00 = max(handles.t0_curr+1,t0);
off0 = t00+handles.labels_plot_off;
off1 = t1+handles.labels_plot_off;
% handles.labels_plot.isstart(off0:off1) = ...
%   handles.data.labelidx(off0:off1)~=0 & ...
%   handles.data.labelidx(off0-1:off1-1)~=handles.data.labelidx(off0:off1);
handles.labels_plot.isstart(off0:off1) = ...
  handles.data.IsLabelStart(handles.expi,handles.flies,t00:t1);

handles = UpdateErrors(handles);

handles = SetNeedSave(handles);

%handles.lastframe_labeled = t;

guidata(handles.figure_JLabel,handles);


function handles = SetLabelsPlot(handles,t0,t1,behavioris)


if t1 < t0,
  tmp = t1;
  t1 = t0;
  t0 = tmp;
end

handles.labels_plot.x(:,t0+handles.labels_plot_off:t1+handles.labels_plot_off,:,:) = nan;
handles.labels_plot.y(:,t0+handles.labels_plot_off:t1+handles.labels_plot_off,:,:) = nan;

for channel = 1:3,
  handles.labels_plot.im(1,t0+handles.labels_plot_off:t1+handles.labels_plot_off,channel) = handles.labelunknowncolor(channel);
end
handles.data.SetLabel(handles.expi,handles.flies,t0:t1,behavioris,0);

for behaviori = 1:handles.data.nbehaviors,

  bidx = find(behaviori == behavioris);
  if isempty(bidx),
    continue;
  end
  for channel = 1:3,
    handles.labels_plot.im(1,t0-1+bidx+handles.labels_plot_off,channel) = handles.labelcolors(behaviori,channel);
  end
  for l = 1:numel(handles.flies),
    ks = t0-1+handles.labels_plot_off+bidx;
    xplot0 = handles.data.GetTrxValues('X1',handles.expi,handles.flies(l),t0-1+bidx);
    xplot1 = handles.data.GetTrxValues('X1',handles.expi,handles.flies(l),min(t0+bidx,handles.t1_curr));
    handles.labels_plot.x(:,ks,behaviori,l) = [xplot0;xplot1];
    yplot0 = handles.data.GetTrxValues('Y1',handles.expi,handles.flies(l),t0-1+bidx);
    yplot1 = handles.data.GetTrxValues('Y1',handles.expi,handles.flies(l),min(t0+bidx,handles.t1_curr));
    handles.labels_plot.y(:,ks,behaviori,l) = [yplot0;yplot1];
  end
  
end

% isstart
if t0 == handles.t0_curr,
  handles.labels_plot.isstart(t0+handles.labels_plot_off) = behavioris(1) ~= 0;
end
t00 = max(handles.t0_curr+1,t0);
off0 = t00+handles.labels_plot_off;
off1 = t1+handles.labels_plot_off;
% handles.labels_plot.isstart(off0:off1) = ...
%   handles.data.labelidx(off0:off1)~=0 & ...
%   handles.data.labelidx(off0-1:off1-1)~=handles.data.labelidx(off0:off1);
handles.labels_plot.isstart(off0:off1) = ...
  handles.data.IsLabelStart(handles.expi,handles.flies,t00:t1);

handles = UpdateErrors(handles);

handles = SetNeedSave(handles);

%handles.lastframe_labeled = t;

guidata(handles.figure_JLabel,handles);


function handles = SetPredictedPlot(handles,t0,t1,behavioris)

if nargin < 2,
  [prediction,t0,t1] = handles.data.GetPredictedIdx(handles.expi,handles.flies);
elseif nargin < 4,
  prediction = handles.data.GetPredictedIdx(handles.expi,handles.flies,t0,t1);
end
behavioris = prediction.predictedidx;
scores = handles.data.NormalizeScores(prediction.scoresidx);
handles.labels_plot.predx(:,t0+handles.labels_plot_off:t1+handles.labels_plot_off,:,:) = nan;
handles.labels_plot.predy(:,t0+handles.labels_plot_off:t1+handles.labels_plot_off,:,:) = nan;

for behaviori = 1:handles.data.nbehaviors,

  bidx = find( (behaviori == behavioris) & ...
      (abs(scores)>handles.data.GetConfidenceThreshold(behaviori)));
  if isempty(bidx),
    continue;
  end
  for l = 1:numel(handles.flies),
    ks = t0-1+handles.labels_plot_off+bidx;
    xplot0 = handles.data.GetTrxValues('X1',handles.expi,handles.flies(l),t0-1+bidx);
    xplot1 = handles.data.GetTrxValues('X1',handles.expi,handles.flies(l),min(t0+bidx,handles.t1_curr));
    handles.labels_plot.predx(:,ks,behaviori,l) = [xplot0;xplot1];
    yplot0 = handles.data.GetTrxValues('Y1',handles.expi,handles.flies(l),t0-1+bidx);
    yplot1 = handles.data.GetTrxValues('Y1',handles.expi,handles.flies(l),min(t0+bidx,handles.t1_curr));
    handles.labels_plot.predy(:,ks,behaviori,l) = [yplot0;yplot1];
  end
  
end

guidata(handles.figure_JLabel,handles);

% --- Executes on button press in togglebutton_label_unknown.
function togglebutton_label_unknown_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_label_unknown (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_label_unknown
if get(hObject,'Value'),
  % toggle on
  handles.label_state = -1; 
  handles.label_imp = [];
  handles.label_t0 = handles.ts(1);

  % set everything else to off
  for j = 1:2*handles.data.nbehaviors,
    if isnan(handles.togglebutton_label_behaviors(j)), continue; end
    set(handles.togglebutton_label_behaviors(j),'Value',0,'Enable','off');
  end

  set(handles.htimeline_label_curr,'XData',handles.label_t0 + [-.5,-.5,.5,.5,-.5],...
    'FaceColor',handles.labelunknowncolor);
  
  UpdatePlots(handles,...
    'refreshim',false,'refreshflies',true,'refreshtrx',false,'refreshlabels',true,...
    'refresh_timeline_manual',true,...
    'refresh_timeline_auto',false,...
    'refresh_timeline_suggest',false,...
    'refresh_timeline_error',true,...
    'refresh_timeline_xlim',false,...
    'refresh_timeline_hcurr',false,...
    'refresh_timeline_props',false,...
    'refresh_timeline_selection',false,...
    'refresh_curr_prop',false);

  set(handles.menu_file,'enable','off');
  set(handles.menu_edit,'enable','off');
  set(handles.menu_go,'enable','off');
  set(handles.menu_classifier,'enable','off');
  

  
%   % set everything else to off
%   for j = 1:handles.data.nbehaviors,
%     set(handles.togglebutton_label_behaviors(j),'Value',0,'String',sprintf('Start %s',handles.data.labelnames{j}));
%   end
%   % set the current frame to be labeled
%   %handles.lastframe_labeled = [];
%   handles = SetLabelPlot(handles,min(handles.t1_curr,max(handles.t0_curr,handles.ts(1))),0);
%   UpdatePlots(handles,'refreshim',false,'refreshtrx',false,'refreshflies',false);
%   set(handles.togglebutton_label_unknown,'String','*Label Unknown*');

else
  
  if handles.ts(1) <= handles.label_t0,
    t0 = handles.ts(1);
    t1 = handles.label_t0;
  else
    t0 = handles.label_t0;
    t1 = handles.ts(1);
  end
  t0 = min(handles.t1_curr,max(handles.t0_curr,t0));
  t1 = min(handles.t1_curr,max(handles.t0_curr,t1));
  handles = SetLabelPlot(handles,t0,t1,0,0);

  handles.label_state = 0;
  handles.label_imp = [];
  handles.label_t0 = [];
  set(handles.htimeline_label_curr,'XData',nan(1,5));
    
  %handles.data.StoreLabels();
  for j = 1:2*handles.data.nbehaviors,
    if isnan(handles.togglebutton_label_behaviors(j)), continue; end
    buttonStr = sprintf('%s',handles.data.labelnames{ceil(j/2)});
    if mod(j,2); buttonStr = sprintf('Important %s',buttonStr); end
    set(handles.togglebutton_label_behaviors(j),'Value',0,'String',buttonStr,'Enable','on');
  end
  set(handles.togglebutton_label_unknown,'Value',0,'String','Unknown','Enable','on');  
  UpdatePlots(handles,...
    'refreshim',false,'refreshflies',true,'refreshtrx',false,'refreshlabels',true,...
    'refresh_timeline_manual',true,...
    'refresh_timeline_auto',false,...
    'refresh_timeline_suggest',false,...
    'refresh_timeline_error',true,...
    'refresh_timeline_xlim',false,...
    'refresh_timeline_hcurr',false,...
    'refresh_timeline_props',false,...
    'refresh_timeline_selection',false,...
    'refresh_curr_prop',false);
  
%   handles.label_state = 0;
%   %handles.data.StoreLabels();
%   set(handles.togglebutton_label_unknown,'String','Start Unknown');
  set(handles.menu_file,'enable','on');
  set(handles.menu_edit,'enable','on');
  set(handles.menu_go,'enable','on');
  set(handles.menu_classifier,'enable','on');
   
end

guidata(hObject,handles);

function edit_framenumber_Callback(hObject, eventdata, handles)
% hObject    handle to edit_framenumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_framenumber as text
%        str2double(get(hObject,'String')) returns contents of edit_framenumber as a double
v = str2double(get(hObject,'String'));
i = GetPreviewPanelNumber(hObject);
if isnan(v) || isempty(v),
  set(hObject,'String',num2str(handles.ts(i)));
else
  v = round(v);
  handles = SetCurrentFrame(handles,i,v,hObject);
  guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function edit_framenumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_framenumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_view_timeline_options_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_timeline_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompts = {'N. frames shown:'};
answers = {num2str(handles.timeline_nframes)};
res = inputdlg(prompts,'Timeline view options',numel(prompts),answers);
if isempty(res); return; end;
handles.timeline_nframes = str2double(res{1});

xlim = [handles.ts(1)-(handles.timeline_nframes-1)/2,...
  handles.ts(1)+(handles.timeline_nframes-1)/2];
for i = 1:numel(handles.axes_timelines),
  set(handles.axes_timelines(i),'XLim',xlim);
  zoom(handles.axes_timelines(i),'reset');
end
guidata(hObject,handles);


% --- Executes on mouse press over axes background.
function axes_preview_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_preview (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% WARNING: this function directly accesses handles.data.trx make sure
% that we've preloaded the right experiment and flies. 
% REMOVED!
if handles.expi ~= handles.data.expi,
  handles.data.Preload(handles.expi,handles.flies);
end

% which preview panel is this
i = GetPreviewPanelNumber(hObject);
nprev = handles.traj_nprev;
npost = handles.traj_npost;
mind = inf;
pt = get(handles.axes_previews(i),'CurrentPoint');
xclick = pt(1,1);
yclick = pt(1,2);
dx = diff(get(handles.axes_previews(i),'XLim'));
dy = diff(get(handles.axes_previews(i),'YLim'));
for j = 1:numel(handles.flies),
  fly = handles.flies(j);
  T0 = handles.data.firstframes_per_exp{handles.expi}(fly);
  T1 = handles.data.endframes_per_exp{handles.expi}(fly);
  t0 = min(T1,max(T0,handles.ts(i)-nprev));
  t1 = min(T1,max(T0,handles.ts(i)+npost));
  %off = handles.data.trx(fly).off;
%   [mindcurr,k] = min( ((handles.data.trx(fly).x(t0+off:t1+off)-xclick)/dx).^2 + ...
%     ((handles.data.trx(fly).y(t0+off:t1+off)-yclick)/dy).^2 );
  [mindcurr,k] = min( ((handles.data.GetTrxValues('X1',handles.expi,fly,t0:t1)-xclick)/dx).^2 + ...
    ((handles.data.GetTrxValues('Y1',handles.expi,fly,t0:t1)-yclick)/dy).^2 );
  if mindcurr < mind,
    mind = mindcurr;
    mint = k+t0-1;
  end
end
if mind <= handles.max_click_dist_preview
  handles = SetCurrentFrame(handles,i,mint,hObject);
end
guidata(hObject,handles);


function fly_ButtonDownFcn(hObject, eventdata, handles, fly, i)

% TODO: figure out how to do this when multiple flies define a behavior

% check for double click
if ~strcmpi(get(handles.figure_JLabel,'SelectionType'),'open') || ...
    numel(handles.flies) == 1 && handles.flies == fly,
  % call the axes button down fcn
  axes_preview_ButtonDownFcn(handles.axes_preview(i), eventdata, handles);
  return;
end

% Dont switch flies when the label pen is down.
penDown = false;
if strcmp(handles.configparams.JLabelMode.mode,'basic'),
  behaviorVals = get(handles.togglebutton_label_behaviors(2:2:end),'Value');
else
  behaviorVals = get(handles.togglebutton_label_behaviors,'Value');
end

for ndx = 1:length(behaviorVals)
  penDown = penDown | behaviorVals{ndx};
end
penDown = penDown | get(handles.togglebutton_label_unknown,'Value');
if penDown, return; end

% check if the user wants to switch to this fly
% TODO: this directly accesses handles.data.labels -- abstract this
[ism,j] = ismember(fly,handles.data.labels(handles.expi).flies,'rows');
if ism,
  nbouts = nnz(~strcmpi(handles.data.labels(handles.expi).names{j},'None'));
else
  nbouts = 0;
end

endframe = handles.data.endframes_per_exp{handles.expi}(fly);
firstframe = handles.data.firstframes_per_exp{handles.expi}(fly);
prompt = {sprintf('Switch to fly %d?',fly),...
  sprintf('Trajectory length = %d',endframe-firstframe+1),...
  sprintf('First frame = %d',firstframe),...
  sprintf('N. bouts labeled: %d',nbouts)};

if handles.data.hassex,
  if handles.data.hasperframesex,
    sexfrac = handles.data.GetSexFrac(handles.expi,fly);
    prompt{end+1} = sprintf('Sex: %d%%M, %d%%F',round(sexfrac.M*100),round(sexfrac.F*100));
  else
    t = max(handles.t0_curr,handles.ts(1));
    sex = handles.data.GetSex(handles.expi,fly,t);
    if iscell(sex),
      sex = sex{1};
    end
    prompt{i} = sprintf('Sex: %s',sex);
  end
end

res = questdlg(prompt,...
  'Change flies?','Yes','No','Yes');

if strcmpi(res,'No'),
  return;
end

SetCurrentFlies(handles,fly);


% --- Executes on mouse press over axes background.
function axes_timeline_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_timeline_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pt = get(hObject,'CurrentPoint');
t = min(max(handles.t0_curr,round(pt(1,1))),handles.t1_curr);%nframes);
% TODO: which axes?
SetCurrentFrame(handles,1,t,hObject);


% --- Executes on button press in pushbutton_train.
function pushbutton_train_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store the current labels to windowdata_labeled
handles.data.StoreLabels();
handles.data.Train(handles.doFastUpdates);
handles = SetPredictedPlot(handles);
% predict for current window
handles = UpdatePrediction(handles);
handles = SetNeedSave(handles);
guidata(hObject,handles);

function handles = UpdatePrediction(handles)

% update prediction for currently shown timeline
% TODO: make this work for multiple axes
t0 = max(handles.t0_curr,floor(handles.ts(1)-handles.timeline_nframes/2));
t1 = min(handles.t1_curr,ceil(handles.ts(1)+7*handles.timeline_nframes/2));
handles.data.Predict(handles.expi,handles.flies,t0:t1);
handles = SetPredictedPlot(handles,t0,t1);

handles = UpdateTimelineIms(handles);
guidata(handles.figure_JLabel,handles);
UpdatePlots(handles,'refreshim',false,'refreshflies',true,...
  'refreshtrx',true,'refreshlabels',true,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);

function handles = UpdateErrors(handles)

% update prediction for currently shown timeline
% TODO: make this work for multiple axes
handles.data.UpdateErrorIdx();
handles = UpdateTimelineIms(handles);
UpdatePlots(handles,'refreshim',false,'refreshflies',...
  false,'refreshtrx',false,'refreshlabels',false,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);

% --- Executes on button press in pushbutton_predict.
function pushbutton_predict_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_predict (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = UpdatePrediction(handles);
guidata(hObject,handles);

% --- Executes on key press with focus on figure_JLabel or any of its controls.
function figure_JLabel_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure_JLabel (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

%disp(eventdata.Character);
%disp(eventdata.Modifier);

function SetStatusCallback(s,h)

handles = guidata(h);
SetStatus(handles,s);

function ClearStatusCallback(h)

handles = guidata(h);
ClearStatus(handles);

function SetStatus(handles,s,isbusy)

if nargin < 3 || isbusy,
  color = handles.busystatuscolor;
else
  color = handles.idlestatuscolor;
end
set(handles.text_status,'ForegroundColor',color,'String',s);

if strcmpi(get(handles.figure_JLabel,'Visible'),'off'),
  msgbox(s,'JLabel Status','modal');
end

function ClearStatus(handles)

set(handles.text_status,'ForegroundColor',handles.idlestatuscolor,'String',handles.status_bar_text);
h = findall(0,'Type','figure','Name','JLabel Status');
if ~isempty(h), delete(h(ishandle(h))); end


% --------------------------------------------------------------------
function menu_file_load_top_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_load_top (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_go_switch_experiment_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_switch_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = cell(1,handles.data.nexps);
for i = 1:handles.data.nexps,
  expStats = handles.data.GetExpStats(i);
  if i == handles.expi,
    s{i} = sprintf('%s, N flies: %d, OPEN NOW',...
      expStats.name,expStats.nflies);
  else
    s{i} = sprintf('%s, N flies: %d',expStats.name,expStats.nflies);
    if handles.data.hassex,
      nmales = sum([handles.data.frac_sex_per_exp{i}.M]);
      nfemales = sum([handles.data.frac_sex_per_exp{i}.F]);
      if handles.data.hasperframesex,
        s{i} = [s{i},sprintf('(%.1f M, %.1f F)',nmales,nfemales)];
      else
        s{i} = [s{i},sprintf('(%d M, %d F)',nmales,nfemales)];
      end
    end
    s{i} = [s{i},sprintf(', N flies labeled: %d, N bouts labeled: %d, last labeled: %s',...
      expStats.nlabeledflies,...
      expStats.nlabeledbouts,...
      expStats.labeldatestr)];
    if ~isempty(expStats.nscoreframes)
      s{i} = [s{i},sprintf(', Frames Predicted as %s:%d, Total Frames Predicted:%d, Classifier used to predict:%s',...
        handles.data.labelnames{1},...
        expStats.nscorepos,...
        expStats.nscoreframes,...
        expStats.classifierfilename)];
    end
    
  end
end
[expi,ok] = listdlg('ListString',s,'SelectionMode','single',...
  'InitialValue',handles.expi,'Name','Switch experiment',...
  'PromptString','Select experiment:',...
  'ListSize',[640,300]);
if ~ok || expi == handles.expi,
  return;
end
[handles,success] = SetCurrentMovie(handles,expi);
if ~success,
  return;
end
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_go_switch_target_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_switch_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: generalize this to multiple flies labeled

nflies = handles.data.nflies_per_exp(handles.expi);
s = cell(1,nflies);
for fly = 1:nflies,
  if fly == handles.flies(1),
    s{fly} = sprintf('Target %3d, CURRENTLY SELECTED',fly);
  else
    flyStats = handles.data.GetFlyStats(handles.expi,fly);
    s{fly} = sprintf('Target %3d, Trajectory length %5d, First frame %5d, N bouts labeled %2d',...
      fly,flyStats.trajLength,flyStats.firstframe,flyStats.nbouts);
    if flyStats.hassex,
      if ~isempty(flyStats.sexfrac),
        s{fly} = [s{fly},sprintf(', Sex: %3d%%M, %3d%%F',...
          round(flyStats.sexfrac.M*100),round(flyStats.sexfrac.F*100))];
      else
        s{fly} = [s{fly},sprintf(', Sex: %s',flyStats.sex{1})];
      end
    end
    if ~isempty(flyStats.nscoreframes)
      s{fly} = [s{fly},sprintf(', Frames Predicted as %s:%d, Total Frames Predicted:%d',...
        handles.data.labelnames{1},flyStats.nscorepos,flyStats.nscoreframes)];
    end
  end
end


[fly,ok] = listdlg('ListString',s,'SelectionMode','single',...
  'InitialValue',handles.flies(1),'Name','Switch target',...
  'PromptString','Select experiment:',...
  'ListSize',[640,300]);
if ~ok || fly == handles.flies(1),
  return;
end
handles = SetCurrentFlies(handles,fly);
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_view_zoom_in_on_fly_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_zoom_in_on_fly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.preview_zoom_mode = 'center_on_fly';
set(setdiff(handles.menu_view_zoom_options,hObject),'Checked','off');
set(hObject,'Checked','on');
ZoomInOnFlies(handles);
guidata(hObject,handles);

function ZoomInOnFlies(handles,is)

% WARNING: this function accesses handles.data.trx directly -- this requires
% the correct experiment to be loaded
% REMOVED!

if nargin < 2,
  is = 1:numel(handles.axes_previews);
end

xs = nan(1,numel(handles.flies));
ys = nan(1,numel(handles.flies));
for i = is,
  firstframes = handles.data.firstframes_per_exp{handles.expi}(handles.flies);
  endframes = handles.data.endframes_per_exp{handles.expi}(handles.flies);
  %inds = handles.ts(i)-firstframes+1;
  for j = 1:numel(handles.flies),
    if handles.ts(i) < firstframes(j) || handles.ts(i) > endframes(j),
      continue;
    end
    xs(j) = handles.data.GetTrxValues('X1',handles.expi,handles.flies(j),handles.ts(i));
    ys(j) = handles.data.GetTrxValues('Y1',handles.expi,handles.flies(j),handles.ts(i));
    %xs(j) = handles.data.trx(handles.flies(j)).x(inds(j));
    %ys(j) = handles.data.trx(handles.flies(j)).y(inds(j));
  end
  if ~all(isnan(xs)) && ~all(isnan(ys)),
    %xlim = [max([.5,xs-handles.zoom_fly_radius(1)]),min([handles.movie_width+.5,xs+handles.zoom_fly_radius(1)])];
    %ylim = [max([.5,ys-handles.zoom_fly_radius(2)]),min([handles.movie_height+.5,ys+handles.zoom_fly_radius(2)])];
    x0 = min(xs) - handles.zoom_fly_radius(1);
    x1 = max(xs) + handles.zoom_fly_radius(1);
    if x1 - x0 + 1 >= handles.movie_width,
      xlim = [.5,handles.movie_width+.5];
    elseif x0 < .5,
      dx = .5 - x0;
      xlim = [.5,x1 + dx];
    elseif x1 > handles.movie_width+.5,
      dx = x1 - (handles.movie_width+.5);
      xlim = [x0-dx,handles.movie_width+.5];
    else
      xlim = [x0,x1];
    end
    y0 = min(ys) - handles.zoom_fly_radius(2);
    y1 = max(ys) + handles.zoom_fly_radius(2);
    if y1 - y0 + 1 >= handles.movie_height,
      ylim = [.5,handles.movie_height+.5];
    elseif y0 < .5,
      dy = .5 - y0;
      ylim = [.5,y1 + dy];
    elseif y1 > handles.movie_height+.5,
      dy = y1 - (handles.movie_height+.5);
      ylim = [y0-dy,handles.movie_height+.5];
    else
      ylim = [y0,y1];
    end    
    set(handles.axes_previews(i),'XLim',xlim,'YLim',ylim);
  end
end

function KeepFliesInView(handles,is)

% WARNING: this function accesses handles.data.trx directly -- this requires
% the correct experiment to be loaded
% REMOVED!

if nargin < 2,
  is = 1:numel(handles.axes_previews);
end

xs = nan(1,numel(handles.flies));
ys = nan(1,numel(handles.flies));
for i = is,
  firstframes = handles.data.firstframes_per_exp{handles.expi}(handles.flies);
  endframes = handles.data.endframes_per_exp{handles.expi}(handles.flies);
  %inds = handles.ts(i)-firstframes+1;
  for j = 1:numel(handles.flies),
    if handles.ts(i) < firstframes(j) || handles.ts(i) > endframes(j),
      continue;
    end
    xs(j) = handles.data.GetTrxValues('X1',handles.expi,handles.flies(j),handles.ts(i));
    ys(j) = handles.data.GetTrxValues('Y1',handles.expi,handles.flies(j),handles.ts(i));
    %xs(j) = handles.data.trx(handles.flies(j)).x(inds(j));
    %ys(j) = handles.data.trx(handles.flies(j)).y(inds(j));
  end
  xlim = get(handles.axes_previews(i),'XLim');
  % a little border at the edge of the image
  border = .1;
  dx = diff(xlim);
  xlim(1) = xlim(1) + dx*border;
  xlim(2) = xlim(2) - dx*border;
  ylim = get(handles.axes_previews(i),'YLim');
  dy = diff(ylim);
  ylim(1) = ylim(1) + dy*border;
  ylim(2) = ylim(2) - dy*border;
  if min(xs) < xlim(1) || min(ys) < ylim(1) || ...
      max(xs) > xlim(2) || max(ys) > ylim(2),
    % center on flies
    newxlim = [max([.5,xs-handles.zoom_fly_radius(1)]),min([handles.movie_width+.5,xs+handles.zoom_fly_radius(1)])];
    newylim = [max([.5,ys-handles.zoom_fly_radius(2)]),min([handles.movie_height+.5,ys+handles.zoom_fly_radius(2)])];
    set(handles.axes_previews(i),'XLim',newxlim,'YLim',newylim);    
  end
end

% --------------------------------------------------------------------
function menu_file_save_labels_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_save_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.SaveLabels();

% --------------------------------------------------------------------
function menu_edit_clear_all_labels_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_clear_all_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = {};
s{end+1} = 'Experiments with labels: ';
for i = 1:numel(handles.data.labelstats),
  if handles.data.labelstats(i).nbouts_labeled > 0,
    s{end+1} = sprintf('%s: %d bouts',handles.data.expnames{i},handles.data.labelstats(i).nbouts_labeled); %#ok<AGROW>
  end
end

res = questdlg(s,'Really delete all labels?','Yes','No','Cancel','Cancel');
if strcmpi(res,'Yes'),
  handles.data.ClearLabels();
  handles = UpdateTimelineIms(handles);
  UpdatePlots(handles);
end

function RecursiveSetKeyPressFcn(hfig)

hchil = findall(hfig,'-property','KeyPressFcn');
goodidx = true(1,numel(hchil));
for i = 1:numel(hchil),
  if strcmpi(get(hchil(i),'Type'),'uicontrol') && strcmpi(get(hchil(i),'Style'),'edit'),
    goodidx(i) = false;
  end
end
set(hchil(goodidx),'KeyPressFcn',get(hfig,'KeyPressFcn'));

% --- Executes on key press with focus on figure_JLabel and none of its controls.
function figure_JLabel_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure_JLabel (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Key,
  
  case 'leftarrow',
    if strcmpi(eventdata.Modifier,'control'),
      menu_go_previous_bout_end_Callback(hObject,eventdata,handles);
    elseif strcmpi(eventdata.Modifier,'shift'),
      menu_go_previous_automatic_bout_end_Callback(hObject,eventdata,handles);
    else
      menu_go_previous_frame_Callback(hObject, eventdata, handles);
    end
     
  case 'rightarrow',
    if strcmpi(eventdata.Modifier,'control'),
      menu_go_next_bout_start_Callback(hObject,eventdata,handles);
    elseif strcmpi(eventdata.Modifier,'shift'),
      menu_go_next_automatic_bout_start_Callback(hObject,eventdata,handles);
    else
      menu_go_next_frame_Callback(hObject, eventdata, handles);
    end
  
  case 'uparrow',
    menu_go_back_X_frames_Callback(hObject, eventdata, handles);
    
  case 'downarrow',
    menu_go_forward_X_frames_Callback(hObject, eventdata, handles);
    
  case handles.label_shortcuts,
    buttonNum = find(strcmp(eventdata.Key,handles.label_shortcuts),1);
    if buttonNum > 2*handles.data.nbehaviors,
      if handles.label_state > 0,
        buttonNum = 2*handles.label_state - handles.label_imp;
        set(handles.togglebutton_label_behaviors(buttonNum),'Value',false);
        togglebutton_label_behavior1_Callback(handles.togglebutton_label_behaviors(buttonNum), eventdata, handles);
        handles = guidata(hObject);
      end
      set(handles.togglebutton_label_unknown,'Value',get(handles.togglebutton_label_unknown,'Value')==0);
      togglebutton_label_unknown_Callback(handles.togglebutton_label_unknown, eventdata, handles);
      return;
    else
      if strcmp(handles.configparams.JLabelMode.mode,'basic') && mod(buttonNum,2); return; end 
      % Don't do anything when important label keys are pressed in the basic mode
      if handles.label_state == -1,
        set(handles.togglebutton_label_unknown,'Value',false);
        togglebutton_label_unknown_Callback(handles.togglebutton_label_unknown, eventdata, handles);
        handles = guidata(hObject);
      elseif handles.label_state > 0 && (2*handles.label_state -handles.label_imp)~= buttonNum,
        prevButtonNum = 2*handles.label_state - handles.label_imp;
        set(handles.togglebutton_label_behaviors(prevButtonNum),'Value',false);
        togglebutton_label_behavior1_Callback(handles.togglebutton_label_behaviors(prevButtonNum), eventdata, handles);
        handles = guidata(hObject);
      end
      set(handles.togglebutton_label_behaviors(buttonNum),'Value',...
        get(handles.togglebutton_label_behaviors(buttonNum),'Value')==0);
      togglebutton_label_behavior1_Callback(handles.togglebutton_label_behaviors(buttonNum), eventdata, handles);
      return;
    end
  case {'esc','escape'},
    if get(handles.togglebutton_label_unknown,'Value') ~= 0,
      set(handles.togglebutton_label_unknown,'Value',0);
      togglebutton_label_unknown_Callback(handles.togglebutton_label_unknown, eventdata, handles);
    else
      for behaviori = 1:handles.data.nbehaviors,
        if get(handles.togglebutton_label_behaviors(behaviori),'Value') ~= 0,
          set(handles.togglebutton_label_behaviors(behaviori),'Value',0);
          togglebutton_label_behavior1_Callback(handles.togglebutton_label_behaviors(behaviori), eventdata, handles);
        end
      end
    end
    
end


% --------------------------------------------------------------------
function menu_go_next_frame_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_next_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;
t = min(max(1,handles.ts(axesi)+1),handles.t1_curr);%handles.nframes);
% set current frame
SetCurrentFrame(handles,axesi,t,hObject);

% --------------------------------------------------------------------
function menu_go_previous_frame_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_previous_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;
t = min(max(handles.t0_curr,handles.ts(axesi)-1),handles.t1_curr);
% set current frame
SetCurrentFrame(handles,axesi,t,hObject);


% --------------------------------------------------------------------
function menu_go_forward_X_frames_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_forward_X_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;
% TODO: hardcoded in 10 as up/down arrow step
t = min(max(handles.t0_curr,handles.ts(axesi)+handles.nframes_jump_go),handles.t1_curr);%nframes);
% set current frame
SetCurrentFrame(handles,axesi,t,hObject);


% --------------------------------------------------------------------
function menu_go_back_X_frames_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_back_X_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;
% TODO: hardcoded in 10 as up/down arrow step
t = min(max(handles.t0_curr,handles.ts(axesi)-handles.nframes_jump_go),handles.t1_curr);%nframes);
% set current frame
SetCurrentFrame(handles,axesi,t,hObject);


% --------------------------------------------------------------------
function menu_go_next_bout_start_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_next_bout_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;

t = handles.NJObj.Manual_bout_start(handles.data,handles.expi,handles.flies,...
  handles.ts(axesi),handles.t0_curr,handles.t1_curr);
if isempty(t); return; end

SetCurrentFrame(handles,axesi,t,hObject);

% --------------------------------------------------------------------
function menu_go_previous_bout_end_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_previous_bout_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;

t = handles.NJObj.Manual_bout_end(handles.data,handles.expi,handles.flies,...
  handles.ts(axesi),handles.t0_curr,handles.t1_curr);
if isempty(t); return; end

SetCurrentFrame(handles,axesi,t,hObject);


% --------------------------------------------------------------------
function menu_go_navigation_preferences_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_navigation_preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'figure_NavigationPreferences') && ishandle(handles.figure_NavigationPreferences),
  figure(handles.figure_NavigationPreferences);
else
  handles.figure_NavigationPreferences = NavigationPreferences(handles.figure_JLabel,handles.NJObj);
  guidata(hObject,handles);
end


% --------------------------------------------------------------------
function menu_edit_label_shortcuts_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_label_shortcuts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompts  = {};
allShortcuts = handles.label_shortcuts;
curShortcuts = {};
for j = 1:2*handles.data.nbehaviors
  if strcmp(handles.configparams.JLabelMode.mode,'basic') && mod(j,2); continue; end
  labelStr = handles.data.labelnames{ceil(j/2)};
  if mod(j,2), labelStr = ['Important ' labelStr];end
  prompts{end+1} = labelStr;
  curShortcuts{end+1} = allShortcuts{j};
end
prompts{end+1} = 'Unknown';
curShortcuts{end+1} = allShortcuts{end};
sh = inputdlg(prompts,'Label Shortcuts',1,curShortcuts);
if isempty(sh),
  return;
end
if strcmp(handles.configparams.JLabelMode.mode,'basic')
  handles.label_shortcuts(2:2:2*handles.data.nbehaviors)= sh(1:handles.data.nbehaviors);
  handles.label_shortcuts(2*handles.data.nbehaviors+1)= sh(handles.data.nbehaviors+1);
else
  handles.label_shortcuts = sh;
end
guidata(hObject,handles);


% --- Executes when figure_JLabel is resized.
function figure_JLabel_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure_JLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'guipos'),
  return;
end

set(handles.figure_JLabel,'Units','pixels');
figpos = get(handles.figure_JLabel,'Position');

minh = 700;
minw = 500;
if figpos(3) < minw || figpos(4) < minh,
  figpos(3:4) = max(figpos(3:4),[minw,minh]);
  set(handles.figure_JLabel,'Position',figpos);
end

labelbuttons_pos = get(handles.panel_labelbuttons,'Position');
select_pos = get(handles.panel_select,'Position');
learn_pos = get(handles.panel_learn,'Position');
similar_pos = get(handles.panel_similar,'Position');
info_pos = get(handles.panel_selection_info,'Position');

width_leftpanels = figpos(3) - handles.guipos.leftborder_leftpanels - ...
  handles.guipos.leftborder_rightpanels - handles.guipos.width_rightpanels - ...
  handles.guipos.rightborder_rightpanels;
h = figpos(4) - handles.guipos.bottomborder_bottompanels - ...
  handles.guipos.topborder_toppanels - handles.guipos.bottomborder_previewpanels;
height_timelines = h*handles.guipos.frac_height_timelines;
height_previews = h - height_timelines;
timelines_pos = [handles.guipos.leftborder_leftpanels,handles.guipos.bottomborder_bottompanels,...
  width_leftpanels,height_timelines];
set(handles.panel_timelines,'Position',timelines_pos);
% TODO: deal with multiple preview panels
preview_pos = [handles.guipos.leftborder_leftpanels,...
  figpos(4) - handles.guipos.topborder_toppanels - height_previews,...
  width_leftpanels,height_previews];
set(handles.panel_previews(1),'Position',preview_pos);

label_pos = [figpos(3) - labelbuttons_pos(3) - handles.guipos.rightborder_rightpanels,...
  figpos(4) - labelbuttons_pos(4) - handles.guipos.topborder_toppanels,...
  labelbuttons_pos(3:4)];
set(handles.panel_labelbuttons,'Position',label_pos);

dy_label_select = labelbuttons_pos(2) - select_pos(2) - select_pos(4);
new_select_pos = [figpos(3) - select_pos(3) - handles.guipos.rightborder_rightpanels,...
  label_pos(2) - select_pos(4) - dy_label_select,...
  select_pos(3:4)];
set(handles.panel_select,'Position',new_select_pos);

%dy_label_select = labelbuttons_pos(2) - select_pos(2) - select_pos(4);
if strcmp(handles.configparams.JLabelMode.mode,'basic')
  set(handles.panel_similar,'Visible','off');
  new_info_pos = [figpos(3) - info_pos(3) - handles.guipos.rightborder_rightpanels,...
    new_select_pos(2) - info_pos(4) - dy_label_select,...
    info_pos(3:4)];
  set(handles.panel_selection_info,'Position',new_info_pos);
else
  new_similar_pos = [figpos(3) - similar_pos(3) - handles.guipos.rightborder_rightpanels,...
    new_select_pos(2) - similar_pos(4) - dy_label_select,...
    similar_pos(3:4)];
  set(handles.panel_similar,'Position',new_similar_pos);
  new_info_pos = [figpos(3) - info_pos(3) - handles.guipos.rightborder_rightpanels,...
    new_similar_pos(2) - info_pos(4) - dy_label_select,...
    info_pos(3:4)];
  set(handles.panel_selection_info,'Position',new_info_pos);
end


new_learn_pos = [figpos(3) - learn_pos(3) - handles.guipos.rightborder_rightpanels,...
  handles.guipos.bottomborder_bottompanels,...
  learn_pos(3:4)];
set(handles.panel_learn,'Position',...
  new_learn_pos);

function handles = GetGUIPositions(handles)

% all axes panels
handles.panel_previews = findobj(handles.figure_JLabel,'-regexp','Tag','panel_axes\d+');
% all preview axes
handles.axes_previews = findobj(handles.figure_JLabel,'Tag','axes_preview');
% all sliders
handles.slider_previews = findobj(handles.figure_JLabel,'Tag','slider_preview');
% all frame number edit boxes
handles.edit_framenumbers = findobj(handles.figure_JLabel,'Tag','edit_framenumber');
% all play buttons
handles.pushbutton_playstops = findobj(handles.figure_JLabel,'Tag','pushbutton_playstop');
% all timelines
handles.axes_timelines = findobj(handles.figure_JLabel,'-regexp','Tag','^axes_timeline.*')';
% handles.labels_timelines = findobj(handles.figure_JLabel,'-regexp','Tag','^timeline_label.*');
% Regex messes the order which makes it difficult to remove the last data axes.
handles.labels_timelines(1,1) = handles.timeline_label_prop1;
handles.labels_timelines(2,1) = handles.timeline_label_automatic;
handles.labels_timelines(3,1) = handles.timeline_label_manual;

handles.axes_timeline_props = findobj(handles.figure_JLabel,'-regexp','Tag','^axes_timeline_prop.*')';
handles.axes_timeline_labels = setdiff(handles.axes_timelines,handles.axes_timeline_props);

if numel(handles.labels_timelines) ~= numel(handles.labels_timelines),
  error('Number of timeline axes does not match number of timeline labels');
end
% sort by y-position
ys = nan(1,numel(handles.axes_timelines));
for i = 1:numel(handles.axes_timelines),
  pos = get(handles.axes_timelines(i),'Position');
  ys(i) = pos(2);
end
[~,order] = sort(ys);
handles.axes_timelines = handles.axes_timelines(order);
% sort by y-position. 
% Don't touch the last 2 labels that are part of manual and automatic timeline
% because they are inside a panel and so pos(2) is relative to the panel.
ys = nan(1,numel(handles.labels_timelines)-2);
for i = 1:(numel(handles.labels_timelines)-2),
  pos = get(handles.labels_timelines(i),'Position');
  ys(i) = pos(2);
end
[~,order] = sort(ys);
temp = handles.labels_timelines(1:end-2);
handles.labels_timelines(1:end-2) = temp(order);

handles.text_timeline_props = nan(size(handles.axes_timeline_props));
handles.text_timelines = nan(size(handles.axes_timelines));
[~,idx] = ismember(handles.axes_timeline_props,handles.axes_timelines);
for ii = 1:numel(handles.axes_timeline_props),
  i = idx(ii);
  t = get(handles.axes_timeline_props(ii),'Tag');
  m = regexp(t,'^axes_timeline_prop(.*)$','tokens','once');
  t2 = ['text_timeline_prop',m{1}];
  handles.text_timeline_props(ii) = handles.(t2);
  handles.text_timelines(i) = handles.text_timeline_props(ii);
end

figpos = get(handles.figure_JLabel,'Position');
panel_labelbuttons_pos = get(handles.panel_labelbuttons,'Position');
% panel_learn_pos = get(handles.panel_learn,'Position');
panel_timelines_pos = get(handles.panel_timelines,'Position');
panel_previews_pos = cell(size(handles.panel_previews));
for i = 1:numel(handles.panel_previews),
  panel_previews_pos{i} = get(handles.panel_previews(i),'Position');
end
handles.guipos.width_rightpanels = panel_labelbuttons_pos(3);
handles.guipos.rightborder_rightpanels = figpos(3) - (panel_labelbuttons_pos(1) + panel_labelbuttons_pos(3));
handles.guipos.leftborder_leftpanels = panel_timelines_pos(1);
handles.guipos.leftborder_rightpanels = panel_labelbuttons_pos(1) - (panel_timelines_pos(1) + panel_timelines_pos(3));
handles.guipos.topborder_toppanels = figpos(4) - (panel_labelbuttons_pos(2) + panel_labelbuttons_pos(4));
if handles.guipos.topborder_toppanels < 0
  handles.guipos.topborder_toppanels = 15;
end
handles.guipos.bottomborder_bottompanels = panel_timelines_pos(2);
handles.guipos.bottomborder_previewpanels = panel_previews_pos{end}(2) - (panel_timelines_pos(2)+panel_timelines_pos(4));
handles.guipos.frac_height_timelines = panel_timelines_pos(4) / (panel_timelines_pos(4) + panel_previews_pos{1}(4));

handles.guipos.timeline_bottom_borders = nan(1,numel(handles.axes_timelines));
handles.guipos.timeline_left_borders = nan(1,numel(handles.axes_timelines));
handles.guipos.timeline_label_middle_offsets = nan(1,numel(handles.axes_timelines));
pos0 = get(handles.axes_timelines(1),'Position');
handles.guipos.timeline_bottom_borders(1) = pos0(2);
handles.guipos.timeline_heights(1) = pos0(4);
handles.guipos.timeline_xpos = pos0(1);
handles.guipos.timeline_rightborder = panel_timelines_pos(3) - pos0(1) - pos0(3);
for i = 2:numel(handles.axes_timelines),
  pos1 = get(handles.axes_timelines(i),'Position');
  handles.guipos.timeline_bottom_borders(i) = pos1(2) - pos0(2) - pos0(4);
  handles.guipos.timeline_heights(i) = pos1(4);
  pos0 = pos1;
end
handles.guipos.timeline_top_border = panel_timelines_pos(4) - pos1(2) - pos1(4);
handles.guipos.timeline_heights = handles.guipos.timeline_heights / sum(handles.guipos.timeline_heights);
for i = 1:numel(handles.axes_timelines),
  ax_pos = get(handles.axes_timelines(i),'Position');
  label_pos = get(handles.labels_timelines(i),'Position');
  handles.guipos.timeline_left_borders(i) = label_pos(1);
  m = ax_pos(2) + ax_pos(4)/2;
  handles.guipos.timeline_label_middle_offsets(i) = label_pos(2)-m;
end
ax_pos = get(handles.axes_timeline_prop1,'Position');
handles.guipos.timeline_prop_height = ax_pos(4);
pos = get(handles.timeline_label_prop1,'Position');
handles.guipos.timeline_prop_label_left_border = pos(1);
handles.guipos.timeline_prop_label_size = pos(3:4);
handles.guipos.timeline_prop_label_callback = get(handles.timeline_label_prop1,'Callback');
handles.guipos.timeline_prop_fontsize = get(handles.timeline_label_prop1,'FontSize');
m = ax_pos(2) + ax_pos(4)/2;
handles.guipos.timeline_prop_label_middle_offset = pos(2)-m;

pos = get(handles.text_timeline_prop1,'Position');
handles.guipos.text_timeline_prop_right_border = ax_pos(1) - pos(1) - pos(3);
handles.guipos.text_timeline_prop_size = pos(3:4);
handles.guipos.text_timeline_prop_middle_offset = pos(2)-m;
handles.guipos.text_timeline_prop_fontsize = get(handles.text_timeline_prop1,'FontSize');
handles.guipos.text_timeline_prop_bgcolor = get(handles.text_timeline_prop1,'BackgroundColor');
handles.guipos.text_timeline_prop_fgcolor = get(handles.text_timeline_prop1,'ForegroundColor');

axes_pos = get(handles.axes_preview,'Position');
slider_pos = get(handles.slider_preview,'Position');
edit_pos = get(handles.edit_framenumber,'Position');
play_pos = get(handles.pushbutton_playstop,'Position');
handles.guipos.preview_axes_top_border = panel_previews_pos{end}(4) - axes_pos(4) - axes_pos(2);
handles.guipos.preview_axes_bottom_border = axes_pos(2);
handles.guipos.preview_axes_left_border = axes_pos(1);
handles.guipos.preview_axes_right_border = panel_previews_pos{end}(3) - axes_pos(1) - axes_pos(3);
handles.guipos.preview_slider_left_border = slider_pos(1);
handles.guipos.preview_slider_right_border = panel_previews_pos{end}(3) - slider_pos(1) - slider_pos(3);
handles.guipos.preview_slider_bottom_border = slider_pos(2);
handles.guipos.preview_play_left_border = play_pos(1) - slider_pos(1) - slider_pos(3);
handles.guipos.preview_play_bottom_border = play_pos(2);
handles.guipos.preview_edit_left_border = edit_pos(1) - play_pos(1) - play_pos(3);
handles.guipos.preview_edit_bottom_border = edit_pos(2);

% --- Executes when panel_timelines is resized.
function panel_timelines_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to panel_timelines (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'axes_timelines'),
  return;
end
panel_pos = get(handles.panel_timelines,'Position');

ntimelines = numel(handles.axes_timelines);
h = panel_pos(4) - sum(handles.guipos.timeline_bottom_borders) - handles.guipos.timeline_top_border;
w = panel_pos(3) - handles.guipos.timeline_rightborder - handles.guipos.timeline_xpos;

y0 = 0;
for i = 1:ntimelines,
  y0 = y0 + handles.guipos.timeline_bottom_borders(i);
  axes_pos = [handles.guipos.timeline_xpos,y0,w,h*handles.guipos.timeline_heights(i)];
  m = axes_pos(2) + axes_pos(4)/2;
  set(handles.axes_timelines(i),'Position',axes_pos);
  label_pos = get(handles.labels_timelines(i),'Position');
  new_label_pos = [handles.guipos.timeline_left_borders(i),...
    m+handles.guipos.timeline_label_middle_offsets(i),...
    label_pos(3:4)];
  set(handles.labels_timelines(i),'Position',new_label_pos);
  if ishandle(handles.text_timelines(i)),
    new_text_pos = [axes_pos(1)-handles.guipos.text_timeline_prop_right_border-handles.guipos.text_timeline_prop_size(1),...
      m+handles.guipos.text_timeline_prop_middle_offset,...
      handles.guipos.text_timeline_prop_size];
    set(handles.text_timelines(i),'Position',new_text_pos);
  end
  y0 = y0 + axes_pos(4);
end

timeline_select_pos = get(handles.panel_timeline_select,'Position');
timeline_manual_pos = get(handles.axes_timelines(end),'Position');
timeline_auto_pos = get(handles.axes_timelines(end-1),'Position');
timeline_select_pos(2) = timeline_auto_pos(2);
timeline_select_pos(4) = timeline_manual_pos(2)-timeline_auto_pos(2)+...
                            timeline_manual_pos(4);
set(handles.panel_timeline_select,'Position',timeline_select_pos);

auto_radio_pos = get(handles.timeline_label_automatic,'Position');
manual_radio_pos = get(handles.timeline_label_manual,'Position');
% Position for the auto and manual radio buttons.
auto_radio_pos(2) = timeline_auto_pos(4)/2-auto_radio_pos(4)/2;
set(handles.timeline_label_automatic,'Position',auto_radio_pos);
manual_radio_pos(2) = timeline_select_pos(4)-auto_radio_pos(2)...
  -manual_radio_pos(4);
set(handles.timeline_label_manual,'Position',manual_radio_pos);

%{
% axes_manual_pos = [handles.guipos.timeline_xpos,...
%   panel_pos(4)-handles.guipos.timeline_bordery-h,w,h];
% set(handles.axes_timeline_manual,'Position',axes_manual_pos);  
% 
% axes_auto_pos = [handles.guipos.timeline_xpos,...
%   axes_manual_pos(2)-handles.guipos.timeline_bordery-h,w,h];
% set(handles.axes_timeline_auto,'Position',axes_auto_pos);  
% 
% text_manual_pos = get(handles.timeline_label_manual,'Position');
% m = axes_manual_pos(2) + axes_manual_pos(4)/2;
% new_text_manual_pos = [text_manual_pos(1),m - text_manual_pos(4)/2,...
%   text_manual_pos(3:4)];
% set(handles.timeline_label_manual,'Position',new_text_manual_pos);
% 
% text_auto_pos = get(handles.timeline_label_automatic,'Position');
% m = axes_auto_pos(2) + axes_auto_pos(4)/2;
% new_text_auto_pos = [text_auto_pos(1),m - text_auto_pos(4)/2,...
%   text_auto_pos(3:4)];
% set(handles.timeline_label_automatic,'Position',new_text_auto_pos);
%}

% --- Executes when panel_axes1 is resized.
function panel_axes1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to panel_axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'panel_previews'),
  return;
end
previewi = find(handles.panel_previews==hObject,1);
if isempty(previewi), 
  return;
end

panel_pos = get(handles.panel_previews(previewi),'Position');

axes_pos = [handles.guipos.preview_axes_left_border,...
  handles.guipos.preview_axes_bottom_border,...
  panel_pos(3) - handles.guipos.preview_axes_left_border - handles.guipos.preview_axes_right_border,...
  panel_pos(4) - handles.guipos.preview_axes_top_border - handles.guipos.preview_axes_bottom_border];
set(handles.axes_previews(previewi),'Position',axes_pos);

max_ry = min(handles.zoom_fly_radius(1) * axes_pos(4) / axes_pos(3),handles.movie_height/2);
max_rx = min(handles.zoom_fly_radius(2) * axes_pos(3) / axes_pos(4),handles.movie_width/2);
ischange = false;
if max_ry > handles.zoom_fly_radius(2),
  handles.zoom_fly_radius(2) = max_ry;
  ischange = true;
elseif max_rx > handles.zoom_fly_radius(1),
  handles.zoom_fly_radius(1) = max_rx;
  ischange = true;
end
if ischange,
  guidata(hObject,handles);
  if strcmpi(handles.preview_zoom_mode,'center_on_fly'),
    ZoomInOnFlies(handles,previewi);
  elseif strcmpi(handles.preview_zoom_mode,'follow_fly'),
    KeepFliesInView(handles,previewi);
  end
end

slider_pos = get(handles.slider_previews(previewi),'Position');
new_slider_pos = [handles.guipos.preview_slider_left_border,...
  handles.guipos.preview_slider_bottom_border,...
  panel_pos(3) - handles.guipos.preview_slider_left_border - handles.guipos.preview_slider_right_border,...
  slider_pos(4)];
set(handles.slider_previews(previewi),'Position',new_slider_pos);

play_pos = get(handles.pushbutton_playstops(previewi),'Position');
new_play_pos = [new_slider_pos(1) + new_slider_pos(3) + handles.guipos.preview_play_left_border,...
  handles.guipos.preview_play_bottom_border,play_pos(3:4)];
set(handles.pushbutton_playstops(previewi),'Position',new_play_pos);


edit_pos = get(handles.edit_framenumbers(previewi),'Position');
new_edit_pos = [new_play_pos(1) + new_play_pos(3) + handles.guipos.preview_edit_left_border,...
  handles.guipos.preview_edit_bottom_border,edit_pos(3:4)];
set(handles.edit_framenumbers(previewi),'Position',new_edit_pos);


% --------------------------------------------------------------------
function menu_view_preview_options_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_preview_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompts = {'Playback Speed (fps):','N. previous positions plotted:',...
  'N. future positions plotted:'};

while true,
  defaults = {num2str(handles.play_FPS),num2str(handles.traj_nprev),...
    num2str(handles.traj_npost)};
  res = inputdlg(prompts,'Preview Options',1,defaults);
  errs = {};
  play_FPS = str2double(res{1});
  if isnan(play_FPS) || play_FPS <= 0,
    errs{end+1} = 'Playback speed must be a positive number'; %#ok<AGROW>
  else
    handles.play_FPS = play_FPS;
  end
  
  traj_nprev = str2double(res{2});
  if isnan(traj_nprev) || traj_nprev < 0 || rem(traj_nprev,1) ~= 0,
    errs{end+1} = 'N. previous positions plotted must be a postive integer'; %#ok<AGROW>
  else
    handles.traj_nprev = traj_nprev;
  end
  
  traj_npost = str2double(res{3});
  if isnan(traj_npost) || traj_npost < 0 || rem(traj_npost,1) ~= 0,
    errs{end+1} = 'N. future positions plotted must be a postive integer'; %#ok<AGROW>
  else
    handles.traj_npost = traj_npost;
  end
  
  if isempty(errs),
    break;
  else
    uiwait(warndlg(errs,'Bad preview options'));
  end
  
end
guidata(hObject,handles);
UpdatePlots(handles,...
  'refreshim',false,'refreshflies',false,'refreshtrx',true,'refreshlabels',true,...
  'refresh_timeline_manual',true,...
     'refresh_timeline_auto',false,...
     'refresh_timeline_suggest',false,...
     'refresh_timeline_error',true,...
     'refresh_timeline_xlim',false,...
     'refresh_timeline_hcurr',false,...
     'refresh_timeline_props',false,...
     'refresh_timeline_selection',false,...
     'refresh_curr_prop',false);


% --- Executes on button press in pushbutton_add_timeline.
function pushbutton_add_timeline_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add_timeline (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

AddPropAxes(handles);

% --- Executes on selection change in timeline_label_prop1.
function timeline_label_prop1_Callback(hObject, eventdata, handles)
% hObject    handle to timeline_label_prop1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns timeline_label_prop1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from timeline_label_prop1
propi = GetTimelinePropNumber(hObject,handles);
v = get(hObject,'Value');
s = handles.timeline_prop_options{v};
if strcmpi(s,handles.timeline_prop_remove_string),
  RemovePropAxes(handles,propi);
elseif strcmpi(s,handles.timeline_prop_help_string),
  
  
else
  prop = find(strcmpi(s,handles.data.allperframefns),1);
  handles.perframepropis(propi) = prop;
  handles.perframeprops{propi} = s;
  [perframedata,T0,T1] = handles.data.GetPerFrameData(handles.expi,handles.flies,prop);
  set(handles.htimeline_data(propi),'XData',T0:T1,'YData',perframedata);
  ylim = [min(perframedata),max(perframedata)];
  set(handles.axes_timeline_props(propi),'YLim',ylim);
  zoom(handles.axes_timeline_props(propi),'reset');
  if ~isnan(handles.timeline_data_ylims(1,prop)),
    ylim = handles.timeline_data_ylims(:,prop);
    set(handles.axes_timeline_props(propi),'YLim',ylim);
  end
  ydata = [ylim(1)+diff(ylim)*.025,ylim(2)-diff(ylim)*.025];
  set(handles.hselection(propi),'YData',ydata([1,2,2,1,1]));
  guidata(hObject,handles);
end

function i = GetTimelinePropNumber(hObject,handles)

t = get(hObject,'Type');
if strcmpi(t,'axes'),
  i = find(hObject == handles.axes_timeline_props,1);
elseif strcmpi(t,'uicontrol'),
  s = get(hObject,'Style');
  if strcmpi(s,'popupmenu'),
    j = find(hObject == handles.labels_timelines,1);
    if isempty(j),
      i = [];
    else
      i = find(handles.axes_timelines(j) == handles.axes_timeline_props,1);
    end
  else
    i = [];
  end
else
  i = [];
end
if isempty(i),
  warning('Could not find index of parent panel');
  i = 1;
end

function handles = RemovePropAxes(handles,propi)

% which axes
axi = find(handles.axes_timelines == handles.axes_timeline_props(propi));

% how much height we will remove
axes_pos = get(handles.axes_timeline_props(propi),'Position');
hremove = axes_pos(4) + handles.guipos.timeline_bottom_borders(axi+1);

% set the sizes of the other axes to stretch
Z0 = sum(handles.guipos.timeline_heights);
Z1 = Z0 - handles.guipos.timeline_heights(axi);
handles.guipos.timeline_heights = handles.guipos.timeline_heights * Z0 / Z1;

% delete the axes
delete(handles.axes_timeline_props(propi));
delete(handles.labels_timelines(axi));
if ishandle(handles.text_timelines(axi)),
  delete(handles.text_timelines(axi));
end
handles.axes_timeline_props(propi) = [];
handles.axes_timelines(axi) = [];
handles.labels_timelines(axi) = [];
handles.text_timelines(axi) = [];
handles.htimeline_data(propi) = [];
handles.hcurr_timelines(axi) = [];
handles.hselection(axi) = [];
handles.guipos.timeline_bottom_borders(axi+1) = [];
handles.guipos.timeline_heights(axi) = [];
handles.guipos.timeline_left_borders(axi) = [];
handles.guipos.timeline_label_middle_offsets(axi) = [];
% handles.perframepropfns(propi) = [];
handles.perframepropis(propi) = [];

% show the xticks
set(handles.axes_timelines(1),'XTickLabelMode','auto');

guidata(handles.figure_JLabel,handles);

% make the panel smaller
panel_timelines_pos = get(handles.panel_timelines,'Position');
panel_timelines_pos(4) = panel_timelines_pos(4) - hremove;
set(handles.panel_timelines,'Position',panel_timelines_pos);

% make the preview panel bigger
panel_previews_pos = get(handles.panel_previews,'Position');
panel_previews_pos(2) = panel_previews_pos(2) - hremove;
panel_previews_pos(4) = panel_previews_pos(4) + hremove;
set(handles.panel_previews,'Position',panel_previews_pos);

handles = guidata(handles.figure_JLabel);

handles.guipos.frac_height_timelines = panel_timelines_pos(4) / (panel_timelines_pos(4) + panel_previews_pos(4));
guidata(handles.figure_JLabel,handles);


function handles = AddPropAxes(handles,prop)

% choose a property
if nargin < 2,
  prop = find(~ismember(1:numel(handles.data.allperframefns),handles.perframepropis),1);
  if isempty(prop),
    prop = 1;
  end
end
propi = numel(handles.axes_timeline_props)+1;
% how much height we will add
hadd = handles.guipos.timeline_prop_height + handles.guipos.timeline_bottom_borders(2);

% set the sizes of the other axes to shrink
panel_pos = get(handles.panel_timelines,'Position');
Z0 = panel_pos(4) - sum(handles.guipos.timeline_bottom_borders) - handles.guipos.timeline_top_border;
Z1 = Z0 + hadd;
handles.guipos.timeline_heights = handles.guipos.timeline_heights * Z0 / Z1;

% add the axes
w = panel_pos(3) - handles.guipos.timeline_rightborder - handles.guipos.timeline_xpos;
ax_pos = [handles.guipos.timeline_xpos,handles.guipos.timeline_bottom_borders(1),...
  w,handles.guipos.timeline_prop_height];
hax = axes('Parent',handles.panel_timelines,'Units','pixels',...
  'Position',ax_pos,'XColor','w','YColor','w',...
  'Color',get(handles.panel_timelines,'BackgroundColor'),...
  'Tag',sprintf('timeline_axes_prop%d',propi));
handles.axes_timeline_props = [hax,handles.axes_timeline_props];
handles.axes_timelines = [hax,handles.axes_timelines];
% fcn = get(handles.axes_timelines(1),'ButtonDownFcn');
% set(hax,'ButtonDownFcn',fcn);
setAxesZoomMotion(handles.hzoom,hax,'vertical');
hold(hax,'on');
[perframedata,T0,T1] = handles.data.GetPerFrameData(handles.expi,handles.flies,prop);
maxylim = [min(perframedata),max(perframedata)];
hdata = plot(T0:T1,perframedata,'w.-');
handles.htimeline_data = [hdata,handles.htimeline_data];
xlim = get(handles.axes_timelines(2),'XLim');
if isnan(handles.timeline_data_ylims(1,prop)),
  ylim = maxylim;
else
  ylim = handles.timeline_data_ylims(:,prop)';
end
set(hax,'XLim',xlim,'YLim',ylim);
zoom(hax,'reset');
hcurr = plot(hax,[0,0]+handles.ts(1),[-10^6,10^6],'y-','HitTest','off','linewidth',2);
handles.hcurr_timelines = [hcurr,handles.hcurr_timelines];
ydata = [ylim(1)+diff(ylim)*.025,ylim(2)-diff(ylim)*.025];
hselection = plot(hax,handles.selected_ts([1,1,2,2,1]),ydata([1,2,2,1,1]),'--',...
  'color',handles.selection_color,...
  'HitTest','off',...
  'LineWidth',3);
handles.hselection = [hselection,handles.hselection];
linkaxes(handles.axes_timelines,'x');

% add the label
m = ax_pos(2)+ax_pos(4)/2; 
pos = [handles.guipos.timeline_prop_label_left_border,...
  m+handles.guipos.timeline_prop_label_middle_offset,...
  handles.guipos.timeline_prop_label_size];
hlabel = uicontrol(handles.panel_timelines,...
  'Style','popupmenu',...
  'Units','pixels',...
  'BackgroundColor',get(handles.labels_timelines(1),'BackgroundColor'),...
  'ForegroundColor',get(handles.labels_timelines(1),'ForegroundColor'),...
  'String',handles.timeline_prop_options,...
  'Value',prop+2,...
  'Position',pos,...
  'FontUnits','pixels',...
  'FontSize',handles.guipos.timeline_prop_fontsize,...
  'Tag',sprintf('timeline_label_prop%d',propi));
set(hlabel,'Callback',@(hObject,eventdata) timeline_label_prop1_Callback(hObject,eventdata,guidata(hObject)));

handles.labels_timelines = [hlabel;handles.labels_timelines];

% add new axes sizes
handles.guipos.timeline_heights = [ax_pos(4) / Z1,handles.guipos.timeline_heights];
handles.guipos.timeline_bottom_borders = handles.guipos.timeline_bottom_borders([1,2,2:numel(handles.guipos.timeline_bottom_borders)]);
handles.guipos.timeline_left_borders = [pos(1),handles.guipos.timeline_left_borders];
handles.guipos.timeline_label_middle_offsets = [handles.guipos.timeline_prop_label_middle_offset,handles.guipos.timeline_label_middle_offsets];
% handles.perframepropfns = [handles.data.allperframefns(prop),handles.perframepropfns];
handles.perframepropis = [prop,handles.perframepropis];

% add the text box
pos = [ax_pos(1)-handles.guipos.text_timeline_prop_right_border-handles.guipos.text_timeline_prop_size(1),...
  m+handles.guipos.text_timeline_prop_middle_offset,...
  handles.guipos.text_timeline_prop_size];
htext = uicontrol(handles.panel_timelines,...
  'Style','text',...
  'Units','pixels',...
  'BackgroundColor',handles.guipos.text_timeline_prop_bgcolor,...
  'ForegroundColor',handles.guipos.text_timeline_prop_fgcolor,...
  'String','????????',...
  'Position',pos,...
  'FontUnits','pixels',...
  'FontSize',handles.guipos.text_timeline_prop_fontsize,...
  'Tag',sprintf('text_timeline_prop%d',propi),...
  'HorizontalAlignment','right');

handles.text_timeline_props = [htext;handles.text_timeline_props];
handles.text_timelines = [htext,handles.text_timelines];


% hide the xtick labels
set(handles.axes_timelines(2),'XTickLabel',{});


guidata(handles.figure_JLabel,handles);

% make the panel bigger
panel_timelines_pos = get(handles.panel_timelines,'Position');
panel_timelines_pos(4) = panel_timelines_pos(4) + hadd;
set(handles.panel_timelines,'Position',panel_timelines_pos);

% make the preview panel smaller
panel_previews_pos = get(handles.panel_previews,'Position');
panel_previews_pos(2) = panel_previews_pos(2) + hadd;
panel_previews_pos(4) = panel_previews_pos(4) - hadd;
set(handles.panel_previews,'Position',panel_previews_pos);

handles = guidata(handles.figure_JLabel);

handles.guipos.frac_height_timelines = panel_timelines_pos(4) / (panel_timelines_pos(4) + panel_previews_pos(4));

guidata(handles.figure_JLabel,handles);

UpdatePlots(handles,...
  'refreshim',false,'refreshflies',false,'refreshtrx',false,'refreshlabels',false,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_auto',false,...
  'refresh_timeline_suggest',false,...
  'refresh_timeline_error',false,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_props',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',true);


function PostZoomCallback(hObject,eventdata,handles)

timelinei = find(eventdata.Axes == handles.axes_timelines,1);
previewi = find(eventdata.Axes == handles.axes_previews,1);
if ~isempty(timelinei),
  for propj = 1:numel(handles.perframepropis),
    prop = handles.perframepropis(propj);
    ylim = get(eventdata.Axes,'YLim');
    handles.timeline_data_ylims(:,prop) = ylim;
    ydata = [ylim(1)+diff(ylim)*.025,ylim(2)-diff(ylim)*.025];
    set(handles.hselection(propj),'YData',ydata([1,2,2,1,1]));
  end
  guidata(eventdata.Axes,handles);
elseif ismember(eventdata.Axes,handles.axes_timeline_labels),
  xlim = get(eventdata.Axes,'XLim');
  handles.timeline_nframes = max(1,round(diff(xlim)-1)/2);
  guidata(eventdata.Axes,handles);
elseif ~isempty(previewi),
  xlim = get(eventdata.Axes,'XLim');
  ylim = get(eventdata.Axes,'YLim');
  rx = round((diff(xlim)-1)/2);
  ry = round((diff(ylim)-1)/2);
  axes_pos = get(eventdata.Axes,'Position');
  max_ry = min(rx * axes_pos(4) / axes_pos(3),handles.movie_height/2);
  max_rx = min(ry * axes_pos(3) / axes_pos(4),handles.movie_width/2);
  ischange = false;
  if max_ry > ry,
    ry = max_ry;
    ischange = true;
  elseif max_rx > rx,
    rx = max_rx;
    ischange = true;
  end
  if rx ~= handles.zoom_fly_radius(1) || ...
      ry ~= handles.zoom_fly_radius(2),
    handles.zoom_fly_radius = [rx,ry];
    guidata(eventdata.Axes,handles);
  end  
  if ischange,
    if strcmpi(handles.preview_zoom_mode,'center_on_fly'),
      ZoomInOnFlies(handles,previewi);
    elseif strcmpi(handles.preview_zoom_mode,'follow_fly'),
      KeepFliesInView(handles,previewi);
    end
  end
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure_JLabel_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure_JLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hchil = gco;
if ismember(hchil,handles.axes_timelines),
  seltype = get(hObject,'SelectionType');
  switch lower(seltype),
    case 'normal', %left
      pt = get(hchil,'CurrentPoint');
      handles.buttondown_t0 = round(pt(1,1));
      handles.buttondown_axes = hchil;
      
      handles.didclearselection = ~any(isnan(handles.selected_ts));
      if handles.didclearselection,
        pushbutton_clearselection_Callback(hObject, eventdata, handles);
      else
        guidata(hObject,handles);
      end
      
      %fprintf('buttondown at %d\n',handles.buttondown_t0);
      %handles.selection_t0 = nan;
      %handles.selection_t1 = nan;

    case {'alternate','extend'}, %right,middle
      pt = get(hchil,'CurrentPoint');
      t = pt(1,1);
      if t >= handles.selected_ts(1) && t <= handles.selected_ts(2),
      end
    case 'open', % double click
  end
end


% --- Executes on mouse motion over figure - except title and menu.
function figure_JLabel_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure_JLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~ishandle(handles.buttondown_axes),
  return;
end
if ~isnan(handles.buttondown_t0) && isnan(handles.selection_t0) && ...
    isnan(handles.selection_t1),    
  handles.selection_t0 = handles.buttondown_t0;
  handles.buttondown_t0 = nan;
  if handles.selecting,
    set(handles.togglebutton_select,'Value',0);
    handles.selecting = false;
  end
  guidata(hObject,handles);
end
if ~isnan(handles.selection_t0),
  pt = get(handles.buttondown_axes,'CurrentPoint');
  handles.selection_t1 = round(pt(1,1));
  handles.selected_ts = [handles.selection_t0,handles.selection_t1];
  %fprintf('Selecting %d to %d\n',handles.selection_t0,handles.selection_t1);
  guidata(hObject,handles);
  UpdateSelection(handles);
end
  


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure_JLabel_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure_JLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~ishandle(handles.buttondown_axes),
  return;
end
if isnan(handles.selection_t0),
  h = handles.buttondown_axes;
  handles.buttondown_axes = nan;
  handles.selection_t0 = nan;
  handles.selection_t1 = nan;
  if ~handles.didclearselection,
    axes_timeline_ButtonDownFcn(h, eventdata, handles);
  end
  return;
end
if ~isnan(handles.selection_t0),
  pt = get(handles.buttondown_axes,'CurrentPoint');
  handles.selection_t1 = round(pt(1,1));
  ts = sort([handles.selection_t0,handles.selection_t1]);
  ts(1) = min(max(ts(1),handles.t0_curr),handles.t1_curr);
  ts(2) = min(max(ts(2),handles.t0_curr),handles.t1_curr);
  if ts(1) == ts(2); % outside the range.
    handles.selected_ts = nan(1,2);
  else
    handles.selected_ts = ts;
  end
  %fprintf('Selected %d to %d\n',handles.selected_ts);
  UpdateSelection(handles);
end
handles.buttondown_axes = nan;
handles.selection_t0 = nan;
handles.selection_t1 = nan;
guidata(hObject,handles);

function UpdateSelection(handles)

tmp = handles.selected_ts + .5*[-1,1];
set(handles.hselection,'XData',tmp([1,1,2,2,1]));
buttons = [handles.pushbutton_playselection,handles.pushbutton_clearselection];
if any(isnan(handles.selected_ts)),
  set(buttons,'Enable','off');
else
  set(buttons,'Enable','on');
end

% --- Executes on button press in togglebutton_select.
function togglebutton_select_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_select
if get(hObject,'Value'),
  handles.selecting = true;
  handles.selected_ts = handles.ts(1)+[0,0];
  handles.buttondown_axes = nan;
  UpdateSelection(handles);
else
  handles.selecting = false;
end
guidata(hObject,handles);

% --- Executes on button press in pushbutton_clearselection.
function pushbutton_clearselection_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clearselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.hplaying == handles.pushbutton_playselection,
  handles = stopPlaying(handles);
end

handles.selected_ts = nan(1,2);
handles.buttondown_axes = nan;
handles.selection_t0 = nan;
handles.selection_t1 = nan;
guidata(hObject,handles);
UpdateSelection(handles);


% --- Executes on button press in pushbutton_playstop.
function pushbutton_playstop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playstop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.hplaying == hObject,
  stopPlaying(handles);
else
  if ~isnan(handles.hplaying),
    stopPlaying(handles);
  end
  play(hObject,handles);
end

% --- Executes on button press in pushbutton_playselection.
function pushbutton_playselection_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playselection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.hplaying == hObject,
  stopPlaying(handles);
else
  if ~isnan(handles.hplaying),
    stopPlaying(handles);
  end
  play(hObject,handles,handles.selected_ts(1),handles.selected_ts(2),true);
end

function predictTimerCallback(obj,event,hObject,framesPerTick)
  global PLAY_TIMER_DONE CALC_FEATURES;
  CALC_FEATURES = true;
  handles = guidata(hObject);
  t0 = handles.ts(1)+framesPerTick;%+handles.timeline_nframes;
  handles.data.Predict(handles.expi,handles.flies,t0:(t0+framesPerTick));
  PLAY_TIMER_DONE = true;
  
function handles = play(hObject,handles,t0,t1,doloop)

clear global PLAY_TIME_DONE CALC_FEATURES
global PLAY_TIMER_DONE CALC_FEATURES;
PLAY_TIMER_DONE = false;
CALC_FEATURES = false;

axi = 1;
set(hObject,'String','Stop','BackgroundColor',[.5,0,0]);
SetButtonImage(handles.pushbutton_playstop);

handles.hplaying = hObject;
guidata(hObject,handles);
ticker = tic;
if nargin < 3,
  t0 = handles.ts(axi);
  t1 = handles.t1_curr;%nframes;
  doloop = false;
end

if ~doloop
  framesPerTick = round(handles.timeline_nframes/4);
  T = timer('TimerFcn',{@predictTimerCallback,hObject,framesPerTick},...
        'Period',framesPerTick/handles.play_FPS,...
        'ExecutionMode','fixedRate',...
        'Tag','predictTimer');
  start(T);
end

while true,
  handles = guidata(hObject);
  if handles.hplaying ~= hObject,
    return;
  end
  
  if CALC_FEATURES
    t0 = handles.ts(axi);
    CALC_FEATURES = false;
  end
  
  if PLAY_TIMER_DONE
    ticker = tic;
    PLAY_TIMER_DONE = false;
    predictStart = max(handles.t0_curr,floor(handles.ts(1)-handles.timeline_nframes/2));
    predictEnd = min(handles.t1_curr,ceil(handles.ts(1)+handles.timeline_nframes/2));
    handles = SetPredictedPlot(handles,predictStart,predictEnd);
    handles = UpdateTimelineIms(handles);

    guidata(hObject,handles);
    UpdatePlots(handles,'refreshim',false,'refreshflies',true,...
      'refreshtrx',true,'refreshlabels',true,...
      'refresh_timeline_manual',false,...
      'refresh_timeline_xlim',false,...
      'refresh_timeline_hcurr',false,...
      'refresh_timeline_selection',false,...
      'refresh_curr_prop',false);

  end
  % how long has it been
  dt_sec = toc(ticker);
  % wait until the next frame should be played
  dt = dt_sec*handles.play_FPS;
  t = ceil(dt)+t0;
  if t > t1,
    if doloop,
      ticker = tic;
      continue;
    else
      handles.hplaying = nan;
      guidata(hObject,handles);
      break;
    end
  end
  SetCurrentFrame(handles,axi,t,hObject);
  handles = UpdateTimelineIms(handles);
  dt_sec = toc(ticker);
  pause_time = (t-t0)/handles.play_FPS - dt_sec;
  if pause_time <= 0,
    drawnow;
  else
    pause(pause_time);
  end
end

stopPlaying(handles);

function handles = stopPlaying(handles)

clear global PLAY_TIMER_DONE;
set(handles.hplaying,'String','Play','BackgroundColor',[.2,.4,0]);
SetButtonImage(handles.hplaying);
T = timerfind('Tag','predictTimer');
if ~isempty(T),  stop(T(:)); delete(T(:)); end
  
hObject = handles.hplaying;
handles.hplaying = nan;
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_view_plottracks_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_plottracks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v = get(handles.menu_view_plottracks,'Checked');

if strcmpi(v,'on'),
  h = findall(handles.axes_previews,'Type','line','Visible','on');
  handles.tracks_visible = h;
  set(h,'Visible','off');
  set(hObject,'Checked','off');
else
  handles.tracks_visible = handles.tracks_visible(ishandle(handles.tracks_visible));
  set(handles.tracks_visible(:),'Visible','on');
  set(hObject,'Checked','on');
end
guidata(hObject,handles);


% --------------------------------------------------------------------
function contextmenu_timeline_manual_go_next_bout_start_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_manual_go_next_bout_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_go_next_bout_start_Callback(hObject,eventdata,handles);


% --------------------------------------------------------------------
function contextmenu_timeline_manual_go_previous_bout_end_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_manual_go_previous_bout_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_go_previous_bout_end_Callback(hObject,eventdata,handles);

function [t0,t1,labelidx,label] = GetBoutProperties(handles,t,labeltype)

if nargin < 3,
  labeltype = 'manual';
end

if t < handles.t0_curr && t > handles.t1_curr,
  t0 = nan;
  t1 = nan;
  labelidx = nan;
  label = '';
  return;
end

if strcmpi(labeltype,'manual'),
  [labelidxStruct,T0,T1] = handles.data.GetLabelIdx(handles.expi,handles.flies);
  labelidx = labelidxStruct.vals;
else
  [prediction,T0,T1] = handles.data.GetPredictedIdx(handles.expi,handles.flies);
  labelidx = prediction.predictedidx;
end
i = t - T0 + 1;
i0 = find(labelidx(1:i-1) ~= labelidx(i),1,'last');
if isempty(i0),
  t0 = T0;
else
  t0 = i0 + T0;
end
i1 = find(labelidx(i+1:end) ~= labelidx(i),1,'first');
if isempty(i1),
  t1 = T1;
else
  t1 = i1 + t - 1;
end
labelidx = labelidx(i);
if nargout >= 4,
  if labelidx == 0,
    label = 'Unknown';
  else
    label = handles.data.labelnames{labelidx};
  end
  if ~strcmpi(labeltype,'manual'),
    label = ['Predicted ',label];
  end
end

% --------------------------------------------------------------------
function contextmenu_timeline_manual_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%disp('hi');
pt = get(handles.axes_timeline_manual,'CurrentPoint');
t = pt(1,1);

% inside a bout?
if t >= handles.t0_curr && t <= handles.t1_curr,
  [handles.bookmark_info.t0,handles.bookmark_info.t1,...
    handles.bookmark_info.labelidx,handles.bookmark_info.label] = ...
    GetBoutProperties(handles,round(t));
  s = sprintf('Bookmark %s bout (%d:%d)',handles.bookmark_info.label,...
    handles.bookmark_info.t0,handles.bookmark_info.t1);  
  set(handles.contextmenu_timeline_manual_bookmark_bout,'Visible','on',...
    'Label',s);
else
  set(handles.contextmenu_timeline_manual_bookmark_bout,'Visible','off');
end
  
% inside the current selection?
if t >= handles.selected_ts(1) && t <= handles.selected_ts(2),
  s = sprintf('Bookmark selection (%d:%d)',handles.selected_ts);
  handles.bookmark_info.t0 = min(handles.selected_ts);
  handles.bookmark_info.t1 = max(handles.selected_ts);
  handles.bookmark_info.labelidx = nan;
  handles.bookmark_info.label = 'Selection';
  set(handles.contextmenu_timeline_manual_bookmark_selection,'Visible','on','Label',s);
else
  set(handles.contextmenu_timeline_manual_bookmark_selection,'Visible','off');
end

guidata(hObject,handles);

% --------------------------------------------------------------------
function contextmenu_timeline_manual_bookmark_bout_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_manual_bookmark_bout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clip = handles.bookmark_info;
clip.t0 = max(handles.t0_curr,clip.t0-1);
clip.t1 = min(clip.t1+1,handles.t1_curr);%nframes);
AddBookmark(handles,handles.bookmark_info);

% --------------------------------------------------------------------
function contextmenu_timeline_manual_bookmark_selection_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_manual_bookmark_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

labelidxStruct = handles.data.GetLabelIdx(handles.expi,handles.flies,handles.bookmark_info.t0,handles.bookmark_info.t1);
labelidx = labelidxStruct.vals;
handles.bookmark_info.labelidx = unique(labelidx);
tmp = [{'Unknown'},handles.data.labelnames];
if numel(handles.bookmark_info.labelidx) == 1,
  handles.bookmark_info.label = tmp{handles.bookmark_info.labelidx+1};
else
  counts = hist(labelidx,handles.bookmark_info.labelidx);
  pct = round(counts / numel(labelidx) * 100);
  s = '';
  for i = 1:numel(handles.bookmark_info.labelidx),
    s = [s,sprintf('%s(%d%%), ',tmp{handles.bookmark_info.labelidx(i)+1},pct(i))]; %#ok<AGROW>
  end
  s = s(1:end-2);
  handles.bookmark_info.label = s;
end
guidata(hObject,handles);

AddBookmark(handles,handles.bookmark_info);

function handles = AddBookmark(handles,clip)

fprintf('TODO: Create bookmark for %d:%d\n',clip.t0,clip.t1);
flystr = sprintf('%d, ',handles.flies);
flystr = flystr(1:end-2);
SetStatus(handles,sprintf('Saving AVI for experiment %s, flies %s, frames %d to %d...',...
  handles.data.expnames{handles.expi},flystr,clip.t0,clip.t1));

handles = make_jlabel_results_movie(handles,clip.t0,clip.t1);
ClearStatus(handles);

%{
% clip.expi = handles.expi;
% clip.flies = handles.flies;
% clip.preview_zoom_mode = handles.preview_zoom_mode;
% axesi = 1;
% clip.xlim = get(handles.axes_previews(axesi),'XLim');
% clip.ylim = get(handles.axes_previews(axesi),'YLim');
% clip.zoom_fly_radius = handles.zoom_fly_radius;
% for i = 1:numel(handles.flies),
%   fly = handles.flies(i);
%   t0 = min(max(clip.t0,handles.t0_curr),handles.t1_curr);
%   t1 = min(max(clip.t1,handles.t0_curr),handles.t1_curr);
%   if t0 <= t1,
%     [xcurr,ycurr,thetacurr,acurr,bcurr] = ...
%       handles.data.GetTrxPos1(handles.expi,fly,t0:t1);
%   end
%   clip.trx(i).x = [nan(1,t0-clip.t0),xcurr,nan(1,clip.t1-t1)];
%   clip.trx(i).y = [nan(1,t0-clip.t0),ycurr,nan(1,clip.t1-t1)];
%   clip.trx(i).a = [nan(1,t0-clip.t0),acurr,nan(1,clip.t1-t1)];
%   clip.trx(i).b = [nan(1,t0-clip.t0),bcurr,nan(1,clip.t1-t1)];
%   clip.trx(i).theta = [nan(1,t0-clip.t0),thetacurr,nan(1,clip.t1-t1)];
% end
%  
% BookmarkedClips(handles.figure_JLabel,handles.data,'clips',clip);
%}


% --------------------------------------------------------------------
function contextmenu_timeline_manual_timeline_options_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_manual_timeline_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_view_timeline_options_Callback(hObject, eventdata, handles);


% --- Executes on button press in similarFramesButton.
function similarFramesButton_Callback(hObject, eventdata, handles)
% hObject    handle to similarFramesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of similarFramesButton

curTime = handles.ts(1);
handles.data.SimilarFrames(curTime);

function s = GetTargetInfo(handles,fly)

  s = {};
  i = 1;
  s{i} = sprintf('Target: %d',fly);
  i = i + 1;
  if handles.data.hassex,
    t = max(handles.t0_curr,handles.ts(1));
    sex = handles.data.GetSex(handles.expi,fly,t);
    if iscell(sex),
      sex = sex{1};
    end
    if handles.data.hasperframesex,
      sexfrac = handles.data.GetSexFrac(handles.expi,fly);
      s{i} = sprintf('Sex: %s (%d%%M, %d%%F)',sex,round(sexfrac.M*100),round(sexfrac.F*100));
    else
      s{i} = sprintf('Sex: %s',sex);
    end
    i = i + 1;
  end
  endframe = handles.data.endframes_per_exp{handles.expi}(fly);
  firstframe = handles.data.firstframes_per_exp{handles.expi}(fly);
  s{i} = sprintf('Frames: %d-%d',firstframe,endframe);
  i = i + 1;
  s{i} = sprintf(handles.data.expnames{handles.expi});


% --------------------------------------------------------------------
function menu_go_next_automatic_bout_start_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_next_automatic_bout_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;

t = handles.NJObj.JumpToStart(handles.data,handles.expi,handles.flies,...
  handles.ts(axesi),handles.t0_curr,handles.t1_curr);
if isempty(t),  return; end

SetCurrentFrame(handles,axesi,t,hObject);

% --------------------------------------------------------------------
function menu_go_previous_automatic_bout_end_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_previous_automatic_bout_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;

t = handles.NJObj.JumpToEnd(handles.data,handles.expi,handles.flies,...
  handles.ts(axesi),handles.t0_curr,handles.t1_curr);
if isempty(t); return; end

SetCurrentFrame(handles,axesi,t,hObject);


% --------------------------------------------------------------------
function menu_view_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_view_zoom_keep_target_in_view_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_zoom_keep_target_in_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.preview_zoom_mode = 'follow_fly';
set(setdiff(handles.menu_view_zoom_options,hObject),'Checked','off');
set(hObject,'Checked','on');
KeepFliesInView(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_view_zoom_static_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_zoom_static (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.preview_zoom_mode = 'static';
set(setdiff(handles.menu_view_zoom_options,hObject),'Checked','off');
set(hObject,'Checked','on');
guidata(hObject,handles);


% --------------------------------------------------------------------
function contextmenu_timeline_automatic_go_next_bout_start_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_automatic_go_next_bout_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_go_next_automatic_bout_start_Callback(hObject,eventdata,handles);

% --------------------------------------------------------------------
function contextmenu_timeline_automatic_go_previous_bout_end_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_automatic_go_previous_bout_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_go_previous_automatic_bout_end_Callback(hObject,eventdata,handles);

% --------------------------------------------------------------------
function contextmenu_timeline_automatic_bookmark_bout_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_automatic_bookmark_bout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clip = handles.bookmark_info;
clip.t0 = max(handles.t0_curr,clip.t0-1);
clip.t1 = min(clip.t1+1,handles.t1_curr);%nframes);
AddBookmark(handles,handles.bookmark_info);

% --------------------------------------------------------------------
function contextmenu_timeline_automatic_bookmark_selection_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_automatic_bookmark_selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prediction = handles.data.GetPredictedIdx(handles.expi,handles.flies,handles.bookmark_info.t0,handles.bookmark_info.t1);
labelidx = prediction.predictedidx;
handles.bookmark_info.labelidx = unique(labelidx);
tmp = [{'Unknown'},handles.data.labelnames];
if numel(handles.bookmark_info.labelidx) == 1,
  handles.bookmark_info.label = ['Predicted ',tmp{handles.bookmark_info.labelidx+1}];
else
  counts = hist(labelidx,handles.bookmark_info.labelidx);
  pct = round(counts / numel(labelidx) * 100);
  s = 'Predicted ';
  for i = 1:numel(handles.bookmark_info.labelidx),
    s = [s,sprintf('%s(%d%%), ',tmp{handles.bookmark_info.labelidx(i)+1},pct(i))]; %#ok<AGROW>
  end
  s = s(1:end-2);
  handles.bookmark_info.label = s;
end
guidata(hObject,handles);

AddBookmark(handles,handles.bookmark_info);

% --------------------------------------------------------------------
function contextmenu_timeline_automatic_timeline_options_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_automatic_timeline_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_view_timeline_options_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function contextmenu_timeline_automatic_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_automatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pt = get(handles.axes_timeline_auto,'CurrentPoint');
t = pt(1,1);

% inside a bout?
if t >= handles.t0_curr && t <= handles.t1_curr,
  [handles.bookmark_info.t0,handles.bookmark_info.t1,...
    handles.bookmark_info.labelidx,handles.bookmark_info.label] = ...
    GetBoutProperties(handles,round(t),'automatic');
  s = sprintf('Bookmark %s bout (%d:%d)',handles.bookmark_info.label,...
    handles.bookmark_info.t0,handles.bookmark_info.t1);  
  set(handles.contextmenu_timeline_automatic_bookmark_bout,'Visible','on',...
    'Label',s);
  s = sprintf('Accept %s bout (%d:%d)',handles.bookmark_info.label,...
    handles.bookmark_info.t0,handles.bookmark_info.t1);
  set(handles.contextmenu_timeline_automatic_accept_bout,'Visible','on',...
    'Label',s);
else
  set(handles.contextmenu_timeline_automatic_bookmark_bout,'Visible','off');
  set(handles.contextmenu_timeline_automatic_accept_bout,'Visible','off');
end
  
% inside the current selection?
if t >= handles.selected_ts(1) && t <= handles.selected_ts(2),
  s = sprintf('Bookmark selection (%d:%d)',handles.selected_ts);
  handles.bookmark_info.t0 = min(handles.selected_ts);
  handles.bookmark_info.t1 = max(handles.selected_ts);
  handles.bookmark_info.labelidx = nan;
  handles.bookmark_info.label = 'Selection';
  set(handles.contextmenu_timeline_automatic_bookmark_selection,'Visible','on','Label',s);
  s = sprintf('Accept selected suggested labels (%d:%d)',handles.selected_ts);
  set(handles.contextmenu_timeline_automatic_accept_selected,'Visible','on','Label',s);  
else
  set(handles.contextmenu_timeline_automatic_bookmark_selection,'Visible','off');
  set(handles.contextmenu_timeline_automatic_accept_selected,'Visible','off');
end

guidata(hObject,handles);


% --------------------------------------------------------------------
function contextmenu_timeline_automatic_accept_selected_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_automatic_accept_selected (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

t0 = handles.bookmark_info.t0;
t1 = handles.bookmark_info.t1;

t0 = min(handles.t1_curr,max(handles.t0_curr,t0));
t1 = min(handles.t1_curr,max(handles.t0_curr,t1));
prediction = handles.data.GetPredictedIdx(handles.expi,handles.flies,t0,t1);
predictedidx = prediction.predictedidx;
handles = SetLabelsPlot(handles,t0,t1,predictedidx);

UpdatePlots(handles,...
  'refreshim',false,'refreshflies',true,'refreshtrx',false,'refreshlabels',true,...
  'refresh_timeline_manual',true,...
  'refresh_timeline_auto',false,...
  'refresh_timeline_suggest',false,...
  'refresh_timeline_error',true,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_props',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);

guidata(hObject,handles);

% --------------------------------------------------------------------
function contextmenu_timeline_automatic_accept_bout_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_automatic_accept_bout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_view_plot_labels_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_plot_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_view_plot_labels_manual_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_plot_labels_manual (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.plot_labels_manual = true;
handles.plot_labels_automatic = false;
% set(handles.panel_timeline_select,'SelectedObject',handles.timeline_label_manual);
set(handles.timeline_label_manual,'Value',1);
UpdatePlotLabels(handles);
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_view_plot_labels_automatic_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_plot_labels_automatic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.plot_labels_manual = false;
handles.plot_labels_automatic = true;
% set(handles.panel_timeline_select,'SelectedObject',handles.timeline_label_automatic);
set(handles.timeline_label_automatic,'Value',1);
UpdatePlotLabels(handles);
guidata(hObject,handles);

function UpdatePlotLabels(handles)

if handles.plot_labels_manual,
  set(handles.hlabels,'Visible','on');
  set(handles.menu_view_plot_labels_manual,'Checked','on');
  set(handles.timeline_label_manual,'ForegroundColor',handles.emphasiscolor,'FontWeight','bold');
else
  set(handles.hlabels,'Visible','off');
  set(handles.menu_view_plot_labels_manual,'Checked','off');
  set(handles.menu_view_plot_labels_manual,'Checked','off');
  set(handles.timeline_label_manual,'ForegroundColor',handles.unemphasiscolor,'FontWeight','normal');
end
if handles.plot_labels_automatic,
  set(handles.hpredicted,'Visible','on');
  set(handles.menu_view_plot_labels_automatic,'Checked','on');
  set(handles.timeline_label_automatic,'ForegroundColor',handles.emphasiscolor,'FontWeight','bold');
else
  set(handles.hpredicted,'Visible','off');
  set(handles.menu_view_plot_labels_automatic,'Checked','off');
  set(handles.timeline_label_automatic,'ForegroundColor',handles.unemphasiscolor,'FontWeight','normal');
end

UpdatePlots(handles,...
  'refreshim',false,'refreshflies',true,'refreshtrx',false,'refreshlabels',true,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_auto',false,...
  'refresh_timeline_suggest',false,...
  'refresh_timeline_error',false,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_props',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);


% --------------------------------------------------------------------
function contextmenu_timeline_automatic_overlay_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_automatic_overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_view_plot_labels_automatic_Callback(hObject, eventdata, handles);

% --------------------------------------------------------------------
function contextmenu_timeline_manual_overlay_Callback(hObject, eventdata, handles)
% hObject    handle to contextmenu_timeline_manual_overlay (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

menu_view_plot_labels_manual_Callback(hObject, eventdata, handles);


% --------------------------------------------------------------------
function menu_view_show_bookmarked_clips_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_show_bookmarked_clips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clipdir = handles.data.GetFile('clipsdir',handles.expi);
if ispc,
  winopen(clipdir);
else
  web(clipdir,'-browser');
end


% --------------------------------------------------------------------
function menu_edit_compression_preferences_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_compression_preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CompressionPreferences(handles.figure_JLabel);

% --------------------------------------------------------------------
function menu_classifier_confThresholds_Callback(hObject, eventdata, handles)
% hObject    handle to menu_classifier_confThresholds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.confFig = ConfidenceThresholds;
ConfidenceThresholds('SetJLabelData',handles.confFig,handles.data,handles);
for ndx = 1:2
  curVal = handles.data.GetConfidenceThreshold(ndx);
  ConfidenceThresholds('SetConfidenceThreshold',handles.confFig,ndx,curVal);
  ConfidenceThresholds('SetSliderColor',handles.confFig,ndx,handles.labelcolors(ndx,:));
end


% --------------------------------------------------------------------
function menu_classifier_Callback(hObject, eventdata, handles)
% hObject    handle to menu_classifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_classifier_classifyall_Callback(hObject, eventdata, handles)
% hObject    handle to menu_classifier_classifyall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

for ndx = 1:handles.data.nexps,
  handles.data.PredictWholeMovie(ndx);
end


% --- Executes on button press in bagButton.
function bagButton_Callback(hObject, eventdata, handles)
% hObject    handle to bagButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.DoBagging();
set(handles.similarFramesButton,'Enable','on');


% --------------------------------------------------------------------
function menu_classifier_doFastUpdates_Callback(hObject, eventdata, handles)
% hObject    handle to menu_classifier_doFastUpdates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curVal = get(hObject,'Checked');
if strcmp(curVal,'on')
  set(hObject,'Checked','off');
  handles.doFastUpdates = false;
else
  set(hObject,'Checked','on');
  handles.doFastUpdates = true;
end
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_classifier_selFeatures_Callback(hObject, eventdata, handles)
% hObject    handle to menu_classifier_selFeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
handles.data.ShowSelectFeatures();


  
% --- Executes when selected object is changed in panel_timeline_select.
function panel_timeline_select_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in panel_timeline_select 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if  strcmp(get(eventdata.NewValue,'tag'),'timeline_label_manual')
  menu_view_plot_labels_manual_Callback(hObject, eventdata, handles);
else
  menu_view_plot_labels_automatic_Callback(hObject, eventdata, handles);
end


% --------------------------------------------------------------------
function crossValidate_Callback(hObject, eventdata, handles)
% hObject    handle to crossValidate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.StoreLabels();
crossError = handles.data.CrossValidate();
handles = SetPredictedPlot(handles);
handles = UpdatePrediction(handles);
guidata(hObject,handles);
dialogStr = {};
dialogStr{end+1} = sprintf('%28s Predicted %-10s    Predicted %-10s \n',...
  '',handles.data.labelnames{2},handles.data.labelnames{1});
dialogStr{end+1} = sprintf('Labeled Important %-12s      %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{1},...
  crossError.numbers(1,1),crossError.frac(1,1),...
  crossError.numbers(1,2),crossError.frac(1,2));
dialogStr{end+1} = sprintf('Labeled %-12s                %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{1},...
  crossError.numbers(2,1),crossError.frac(2,1),...
  crossError.numbers(2,2),crossError.frac(2,2));
dialogStr{end+1} = sprintf('Labeled Important %-12s      %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{2},...
  crossError.numbers(3,1),crossError.frac(3,1),...
  crossError.numbers(3,2),crossError.frac(3,2));
dialogStr{end+1} = sprintf('Labeled %-12s                %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{2},...
  crossError.numbers(4,1),crossError.frac(4,1),...
  crossError.numbers(4,2),crossError.frac(4,2));
dialogStr{end+1} = sprintf('');
dialogStr{end+1} = sprintf('----------------- On Old Labels --------------------');
dialogStr{end+1} = sprintf('');
dialogStr{end+1} = sprintf('%28s Predicted %-10s    Predicted %-10s \n',...
  '',handles.data.labelnames{2},handles.data.labelnames{1});
dialogStr{end+1} = sprintf('Labeled Important %-12s       %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{1},...
  crossError.oldNumbers(1,1),crossError.oldFrac(1,1),...
  crossError.oldNumbers(1,2),crossError.oldFrac(1,2));
dialogStr{end+1} = sprintf('Labeled %-12s                 %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{1},...
  crossError.oldNumbers(2,1),crossError.oldFrac(2,1),...
  crossError.oldNumbers(2,2),crossError.oldFrac(2,2));
dialogStr{end+1} = sprintf('Labeled Important %-12s       %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{2},...
  crossError.oldNumbers(3,1),crossError.oldFrac(3,1),...
  crossError.oldNumbers(3,2),crossError.oldFrac(3,2));
dialogStr{end+1} = sprintf('Labeled %-12s                 %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{2},...
  crossError.oldNumbers(4,1),crossError.oldFrac(4,1),...
  crossError.oldNumbers(4,2),crossError.oldFrac(4,2));

helpdlg(dialogStr,'Cross Validation error');
% helpdlg(sprintf('Cross Validation error is %.2f%%',crossError*100),'Cross Validation error');


% --------------------------------------------------------------------
function menu_classifier_classifyCurrentFly_Callback(hObject, eventdata, handles)
% hObject    handle to menu_classifier_classifyCurrentFly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
t0 = max(handles.t0_curr,floor(handles.ts(1)-handles.timeline_nframes/2));
t1 = min(handles.t1_curr,ceil(handles.ts(1)+7*handles.timeline_nframes/2));
handles.data.Predict(handles.expi,handles.flies,handles.t0_curr:handles.t1_curr);
handles = SetPredictedPlot(handles,t0,t1);

handles = UpdateTimelineIms(handles);
guidata(handles.figure_JLabel,handles);
UpdatePlots(handles,'refreshim',false,'refreshflies',true,...
  'refreshtrx',true,'refreshlabels',true,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);


% --------------------------------------------------------------------
function menu_file_loadScores_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_loadScores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_loadscorescurrentexpdefault_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_loadscorescurrentexpdefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.LoadScoresDefault(handles.data.expi);
handles = UpdateTimelineIms(handles);
guidata(handles.figure_JLabel,handles);
UpdatePlots(handles,'refreshim',false,'refreshflies',true,...
  'refreshtrx',true,'refreshlabels',true,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);


% --------------------------------------------------------------------
function menu_file_loadscorescurrentexpselect_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_loadscorescurrentexpselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tstring = sprintf('Scores file for %s',handles.data.expnames{handles.data.expi});
[fname,pname,~] = uigetfile('*.mat',tstring);
if ~fname; return; end;
sfn = fullfile(pname,fname);
handles.data.LoadScores(handles.data.expi,sfn);
handles = UpdateTimelineIms(handles);
guidata(handles.figure_JLabel,handles);
UpdatePlots(handles,'refreshim',false,'refreshflies',true,...
  'refreshtrx',true,'refreshlabels',true,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);


% --------------------------------------------------------------------
function menu_file_loadscoresAlldefault_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_loadscoresAlldefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
for ndx = 1:handles.data.nexps,
  handles.data.LoadScoresDefault(ndx);
end
handles = UpdateTimelineIms(handles);
guidata(handles.figure_JLabel,handles);
UpdatePlots(handles,'refreshim',false,'refreshflies',true,...
  'refreshtrx',true,'refreshlabels',true,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);


% --------------------------------------------------------------------
function menu_file_loadscoresAllselect_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_loadscoresAllselect (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tstring = sprintf('Root dir to load scores for all experiments');
fname = uigetdir('*.mat',tstring);
if ~fname; return; end;
scoreFileName = sprintf('scores_%s.mat',handles.data.labelnames{1});

for ndx = 1:handles.data.nexps,
  sfn = fullfile(fname,handles.data.expnames{ndx},scoreFileName);
  if ~exist(sfn,'file')
    warndlg(sprintf('Scores file %s does not exist for exp:%s',...
      scoreFileName,handles.data.expnames{ndx}));
    continue; 
  end
  handles.data.LoadScores(ndx,sfn);
end
handles = UpdateTimelineIms(handles);
guidata(handles.figure_JLabel,handles);
UpdatePlots(handles,'refreshim',false,'refreshflies',true,...
  'refreshtrx',true,'refreshlabels',true,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);


% --------------------------------------------------------------------
function menu_classifier_testnewlabels_Callback(hObject, eventdata, handles)
% hObject    handle to menu_classifier_testnewlabels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

newError = handles.data.TestOnNewLabels();
handles = UpdatePrediction(handles);
if ~isfield(newError,'numbers'), return; end;
dialogStr = {};
if isfield(newError,'classifierfilename'),
  dialogStr{end+1} = sprintf('Classifier used to generate scores:%s',newError.classifierfilename);
end
dialogStr{end+1} = sprintf('%28s Predicted %10s    Predicted %10s \n',...
  '',handles.data.labelnames{2},handles.data.labelnames{1});
dialogStr{end+1} = sprintf('Labeled Important %12s      %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{1},...
  newError.numbers(1,1),newError.frac(1,1),...
  newError.numbers(1,2),newError.frac(1,2));
dialogStr{end+1} = sprintf('Labeled            %12s      %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{1},...
  newError.numbers(2,1),newError.frac(2,1),...
  newError.numbers(2,2),newError.frac(2,2));
dialogStr{end+1} = sprintf('Labeled Important  %12s      %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{2},...
  newError.numbers(3,1),newError.frac(3,1),...
  newError.numbers(3,2),newError.frac(3,2));
dialogStr{end+1} = sprintf('Labeled            %12s      %d(%.2f)          %d(%.2f)\n',...
  handles.data.labelnames{2},...
  newError.numbers(4,1),newError.frac(4,1),...
  newError.numbers(4,2),newError.frac(4,2));

helpdlg(dialogStr,'Performance on new labeled data');


% --------------------------------------------------------------------
function menu_file_load_exps_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_load_exps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filename = handles.data.classifierfilename;
[filename,pathname] = uigetfile('*.mat','Load classifier',filename);
if ~ischar(filename),
  return;
end
classifiername = fullfile(pathname,filename);
handles.data.SetClassifierFileName(classifiername);


% --------------------------------------------------------------------
function menu_file_load_woexps_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_load_woexps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

filename = handles.data.classifierfilename;
[filename,pathname] = uigetfile('*.mat','Load classifier',filename);
if ~ischar(filename),
  return;
end
classifiername = fullfile(pathname,filename);
handles.data.SetClassifierFileNameWoExp(classifiername);


% --------------------------------------------------------------------
function menu_classifier_setclassifierparameters_Callback(hObject, eventdata, handles)
% hObject    handle to menu_classifier_setclassifierparameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ClassifierOptions(handles.data);


% --------------------------------------------------------------------
function menu_classifier_switch_target_Callback(hObject, eventdata, handles)
% hObject    handle to menu_classifier_switch_target (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
changeTargetHandle = ClassifierChange(hObject);
ClassifierChange('initTable',changeTargetHandle);