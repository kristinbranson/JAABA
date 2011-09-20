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

% Last Modified by GUIDE v2.5 20-Sep-2011 00:32:46

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
ClearStatus(handles);

% read configuration
[handles,success] = LoadConfig(handles);
if ~success,
  guidata(hObject,handles);
  delete(handles.figure_JLabel);
  return;
end

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

% all axes panels
handles.panel_previews = findobj(handles.figure_JLabel,'-regexp','Tag','panel_axes\d+');
% all preview axes
handles.axes_previews = findobj(handles.figure_JLabel,'Tag','axes_preview');
handles.axes_preview_curr = 1;
% all sliders
handles.slider_previews = findobj(handles.figure_JLabel,'Tag','slider_preview');
% all frame number edit boxes
handles.edit_framenumbers = findobj(handles.figure_JLabel,'Tag','edit_framenumber');
if numel(handles.axes_previews) > numel(handles.ts),
  handles.ts = [handles.ts,repmat(handles.ts(end),[1,numel(handles.axes_previews)-numel(handles.ts)])];
end
% all timelines
handles.axes_timelines = findobj(handles.figure_JLabel','-regexp','Tag','^axes_timeline.*');

% slider callbacks
for i = 1:numel(handles.slider_previews),
  fcn = get(handles.slider_previews(i),'Callback');
  %set(handles.slider_previews(i),'Callback','');
  handles.hslider_listeners(i) = handle.listener(handles.slider_previews(i),...
    'ActionEvent',fcn);
end

% fly current positions
handles.hflies = zeros(handles.nflies_curr,numel(handles.axes_previews));
% fly path
handles.htrx = zeros(handles.nflies_label,numel(handles.axes_previews));

% choose colors for flies
% TODO: change hard-coded colormap
handles.fly_colors = jet(handles.nflies_curr)*.7;
handles.fly_colors = handles.fly_colors(randperm(handles.nflies_curr),:);

for i = 1:numel(handles.axes_previews),
  cla(handles.axes_previews(i),'reset');

  % image in axes_preview
  handles.himage_previews(i) = imagesc(0,'Parent',handles.axes_previews(i),[0,255]);
  set(handles.himage_previews(i),'HitTest','off');
  axis(handles.axes_previews(i),'image');
  
  set(handles.axes_previews(i),'ButtonDownFcn',@(hObject,eventdata) JLabel('axes_preview_ButtonDownFcn',hObject,eventdata,guidata(hObject)));
  hold(handles.axes_previews(i),'on');

  % labeled behaviors
  handles.hlabels = nan(1,handles.data.nbehaviors);
  handles.hlabelstarts = nan(1,handles.data.nbehaviors);
  for j = 1:handles.data.nbehaviors,
    handles.hlabels(j) = plot(handles.axes_previews(i),nan,nan,'-',...
      'color',handles.labelcolors(j,:),'linewidth',5,'HitTest','off');
    % start of label
    handles.hlabelstarts(j) = plot(handles.axes_previews(i),nan,nan,'v',...
      'color',handles.labelcolors(j,:),'markerfacecolor',handles.labelcolors(j,:),...
      'HitTest','off');
    
    set(handles.axes_previews(i),'Color','k','XColor','w','YColor','w');
    
  end
  
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
  end

end

% TODO: allow colormap options
colormap(handles.axes_preview,gray(256));

% timelines

% manual timeline
timeline_axes_color = get(handles.panel_timelines,'BackgroundColor');
handles.himage_timeline_manual = image(zeros([1,1,3]),'Parent',handles.axes_timeline_manual);
set(handles.himage_timeline_manual,'HitTest','off');
hold(handles.axes_timeline_manual,'on');
handles.htimeline_manual_starts = plot(handles.axes_timeline_manual,nan,nan,'w-','HitTest','off');

% auto timeline
handles.himage_timeline_auto = image(zeros([1,1,3]),'Parent',handles.axes_timeline_auto);
set(handles.himage_timeline_auto,'HitTest','off');
hold(handles.axes_timeline_auto,'on');
handles.htimeline_auto_starts = plot(handles.axes_timeline_auto,nan,nan,'w-','HitTest','off');

% suggest timeline
handles.himage_timeline_suggest = image(zeros([1,1,3]),'Parent',handles.axes_timeline_suggest);
set(handles.himage_timeline_suggest,'HitTest','off');
hold(handles.axes_timeline_suggest,'on');
handles.htimeline_suggest_starts = plot(handles.axes_timeline_suggest,nan,nan,'w-','HitTest','off');

% error timeline
handles.himage_timeline_error = image(zeros([1,1,3]),'Parent',handles.axes_timeline_error);
set(handles.himage_timeline_error,'HitTest','off');
hold(handles.axes_timeline_error,'on');
handles.htimeline_error_starts = plot(handles.axes_timeline_error,nan,nan,'w-','HitTest','off');

for i = 1:numel(handles.axes_timelines),
  set(handles.axes_timelines(i),'YTick',[],'XColor','w','YColor','w','Color',timeline_axes_color);
end

handles.hcurr_timelines = nan(size(handles.axes_timelines));
for i = 1:numel(handles.axes_timelines),
  handles.hcurr_timelines(i) = plot(handles.axes_timelines(i),nan(1,2),[.5,1.5],'y-','HitTest','off','linewidth',2);
end

for i = 1:numel(handles.axes_timelines),
  if handles.axes_timelines(i) ~= handles.axes_timeline_error,
    set(handles.axes_timelines(i),'XTickLabel',{});
  end
end

linkaxes(handles.axes_timelines);

% zoom
handles.hzoom = zoom(handles.figure_JLabel);

for i = 1:numel(handles.axes_timelines),
  setAxesZoomMotion(handles.hzoom,handles.axes_timelines(i),'horizontal');
end

% timeline callbacks
fcn = @(hObject,eventdata) JLabel('axes_timeline_ButtonDownFcn',hObject,eventdata,guidata(hObject));
for i = 1:numel(handles.axes_timelines),
  set(handles.axes_timelines(i),'ButtonDownFcn',fcn);
end


function UpdatePlots(handles,varargin)

% WARNING: we directly access handles.data.trx for speed here -- 
% REMOVED! NOT SO SLOW

[axes,refreshim,refreshflies,refreshtrx,refreshlabels,...
  refresh_timeline_manual,refresh_timeline_auto,refresh_timeline_suggest,refresh_timeline_error,...
  refresh_timeline_xlim,refresh_timeline_hcurr] = ...
  myparse(varargin,'axes',1:numel(handles.axes_previews),...
  'refreshim',true,'refreshflies',true,'refreshtrx',true,'refreshlabels',true,...
  'refresh_timeline_manual',true,...
  'refresh_timeline_auto',true,...
  'refresh_timeline_suggest',true,...
  'refresh_timeline_error',true,...
  'refresh_timeline_xlim',true,...
  'refresh_timeline_hcurr',true);

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
end

if refresh_timeline_manual,
  set(handles.himage_timeline_auto,'CData',handles.labels_plot.im);
end
if refresh_timeline_auto,
  set(handles.himage_timeline_auto,'CData',handles.labels_plot.predicted_im);
end
if refresh_timeline_suggest,
  set(handles.himage_timeline_suggest,'CData',handles.labels_plot.suggested_im);
end
if refresh_timeline_error,
  set(handles.himage_timeline_error,'CData',handles.labels_plot.error_im);
end

if refresh_timeline_xlim,
  xlim = [handles.ts(1)-(handles.timeline_nframes-1)/2,...
    handles.ts(1)+(handles.timeline_nframes-1)/2];
  set(handles.axes_timelines,'XLim',xlim);
end

%drawnow;

for i = axes,
  
  if refreshim,

    % read in current frame
    %image_cache = getappdata(handles.figure_JLabel,'image_cache');
    try
%       j = find(handles.ts(i)==image_cache.ts,1);
%       if isempty(j),
        im = handles.readframe(handles.ts(i));
%       else
%         im = image_cache.ims(:,:,:,j);
%       end
    catch ME
      uiwait(warndlg(sprintf('Could not read frame %d from current movie: %s',handles.ts(i),getReport(ME))));
      return;
    end
  
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
    inbounds = handles.data.firstframes_per_exp{handles.expi} <= handles.ts(i) & ...
      handles.data.endframes_per_exp{handles.expi} >= handles.ts(i);
    set(handles.hflies(~inbounds,i),'XData',nan,'YData',nan);
    for fly = find(inbounds),
      % WARNING: this accesses handles.data.trx directly -- make sure that
      % handles.data.trx is loaded for the correct movie
      % REMOVED! NOT SO SLOW

      t = handles.ts(i);
      [xcurr,ycurr,thetacurr,acurr,bcurr] = ...
        handles.data.GetTrxPos1(handles.expi,fly,t);
      updatefly(handles.hflies(fly,i),...
        xcurr,ycurr,thetacurr,acurr,bcurr);
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
        set(handles.hflies(fly,i),'LineWidth',5);
      else
        set(handles.hflies(fly,i),'LineWidth',2);
      end
    end
    
    if handles.zoom_fly,
      ZoomInOnFlies(handles,i);
    end    
  end

  % update trx
  % TODO: remove hard-coded nprev, npost
  nprev = 25;
  npost = 25;
  if refreshtrx,
    for j = 1:numel(handles.flies),
      fly = handles.flies(j);
      tmp = handles.ts(i) + handles.data.trx(fly).off;
      t0 = handles.data.firstframes_per_exp{handles.expi}(fly);
      t1 = handles.data.endframes_per_exp{handles.expi}(fly);
      ts = max(t0,min(t1,tmp-nprev:tmp+npost));
      set(handles.htrx(j,i),'XData',handles.data.GetTrxX1(handles.expi,fly,ts),...
        'YData',handles.data.GetTrxY1(handles.expi,fly,ts));
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
        set(handles.hlabels(j),'XData',handles.labels_plot.x(handles.labels_plot_off+t0:handles.labels_plot_off+t1,j,k),...
          'YData',handles.labels_plot.y(handles.labels_plot_off+t0:handles.labels_plot_off+t1,j,k));
      end
    end
  end
  
  %drawnow;
  
end

if refresh_timeline_hcurr,
  set(handles.hcurr_timelines,'XData',handles.ts([1,1]));
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
try
  SetStatus(handles,'Opening movie...');
  [handles.readframe,handles.nframes,handles.movie_fid,handles.movieheaderinfo] = ...
    get_readframe_fcn(moviefilename,'interruptible',false);
  im = handles.readframe(1);
  handles.movie_width = size(im,2);
  handles.movie_height = size(im,1);
catch ME,
  uiwait(warndlg(sprintf('Error opening movie file %s: %s',moviefilename,getReport(ME)),'Error setting movie'));
  ClearStatus(handles);
  return;
end

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
if isnan(handles.zoom_fly_radius),
  handles.zoom_fly_radius = nanmean([handles.data.trx.a])*20;
end


handles.expi = expi;

ClearStatus(handles);

% TODO: change hard-coded colormap
% update colors
handles.fly_colors = jet(handles.nflies_curr)*.7;
handles.fly_colors = handles.fly_colors(randperm(handles.nflies_curr),:);

% set flies
handles = SetCurrentFlies(handles,flies,true,false);

% delete old fly current positions
if isfield(handles,'hflies'),
  delete(handles.hflies(ishandle(handles.hflies)));
  handles.hflies = [];
end

% update plotted trx handles, as number of flies will change
handles.hflies = zeros(handles.nflies_curr,numel(handles.axes_previews));
for i = 1:numel(handles.axes_previews),
  % fly current positions
  for fly = 1:handles.nflies_curr,
    handles.hflies(fly,i) = plot(handles.axes_previews(i),nan,nan,'-',...
      'color',handles.fly_colors(fly,:),'linewidth',3,...
      'ButtonDownFcn',@(hObject,eventdata) JLabel('fly_ButtonDownFcn',hObject,eventdata,guidata(hObject),fly,i));
  end
end

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
UpdatePlots(handles);

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
handles.nflies_curr = 0;
% delete old fly current positions
if isfield(handles,'hflies'),
  delete(handles.hflies(ishandle(handles.hflies)));
  handles.hflies = [];
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
handles.labels_plot.suggested_im = zeros([1,n,3]);
handles.labels_plot.error_im = zeros([1,n,3]);
handles.labels_plot.x = nan(n,handles.data.nbehaviors,numel(handles.flies));
handles.labels_plot.y = nan(n,handles.data.nbehaviors,numel(handles.flies));
handles.labels_plot_off = 1-handles.t0_curr;
set([handles.himage_timeline_manual,handles.himage_timeline_auto,...
  handles.himage_timeline_error,handles.himage_timeline_suggest],...
  'XData',[handles.t0_curr,handles.t1_curr]);

labelidx = handles.data.GetLabelIdx(handles.expi,flies,handles.t0_curr,handles.t1_curr);
for flyi = 1:numel(flies),
  fly = flies(flyi);
  x = handles.data.GetTrxX(handles.expi,fly,handles.t0_curr:handles.t1_curr);
  y = handles.data.GetTrxY(handles.expi,fly,handles.t0_curr:handles.t1_curr);
  for behaviori = 1:handles.data.nbehaviors
    % WARNING: accesses labelidx
    % REMOVED!
    idx = labelidx == behaviori;
    handles.labels_plot.x(idx,behaviori,flyi) = x{1}(idx);
    handles.labels_plot.y(idx,behaviori,flyi) = y{1}(idx);
  end
end
handles = UpdateTimelineIms(handles);

% which interval we're currently within
handles.current_interval = [];

% update timelines
set(handles.himage_timeline_manual,'CData',handles.labels_plot.im);
axis(handles.axes_timeline_manual,[handles.t0_curr-.5,handles.t1_curr+.5,.5,1.5]);
% update zoom
for i = 1:numel(handles.axes_timelines),
  zoom(handles.axes_timelines(i),'reset');
end

% update trx colors
for i = 1:numel(handles.axes_previews),
  for j = 1:numel(handles.flies),
    fly = handles.flies(j);
    set(handles.htrx(j,i),'Color',handles.fly_colors(fly,:));
  end
end

% status bar text
[~,expname] = fileparts(handles.data.expdirs{handles.expi});
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

guidata(handles.figure_JLabel,handles);

if doupdateplot,
  UpdatePlots(handles);
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
labelidx = handles.data.GetLabelIdx(handles.expi,handles.flies,handles.t0_curr,handles.t1_curr);
for behaviori = 1:handles.data.nbehaviors
  idx = labelidx == behaviori;
  for channel = 1:3,
    handles.labels_plot.im(1,idx,channel) = handles.labelcolors(behaviori,channel);
  end
end
handles.labels_plot.predicted_im(:) = 0;
for behaviori = 1:handles.data.nbehaviors
  idx = handles.data.predictedidx == behaviori;
  for channel = 1:3,
    handles.labels_plot.predicted_im(1,idx,channel) = handles.labelcolors(behaviori,channel);
  end
end
handles.labels_plot.suggested_im(:) = 0;
for behaviori = 1:handles.data.nbehaviors
  idx = handles.data.suggestedidx == behaviori;
  for channel = 1:3,
    handles.labels_plot.suggested_im(1,idx,channel) = handles.labelcolors(behaviori,channel);
  end
end
handles.labels_plot.error_im(:) = 0;
idx = handles.data.erroridx == 1;
for channel = 1:3,
  handles.labels_plot.error_im(1,idx,channel) = handles.correctcolor(channel);
end
idx = handles.data.erroridx == 2;
for channel = 1:3,
  handles.labels_plot.error_im(1,idx,channel) = handles.incorrectcolor(channel);
end
handles.labels_plot.isstart = ...
  cat(2,labelidx(1)~=0,...
  labelidx(2:end)~=0 & ...
  labelidx(1:end-1)~=labelidx(2:end));


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
  if handles.label_state ~= 0,
    handles = SetLabelPlot(handles,min(handles.t1_curr,max(handles.t0_curr,t)),handles.label_state);
  end
  
  % update slider
  if hObject ~= handles.slider_previews(i),
    set(handles.slider_previews(i),'Value',t);
  end

  % update frame number edit box
  if hObject ~= handles.edit_framenumbers(i),
    set(handles.edit_framenumbers(i),'String',num2str(t));
  end
  
  guidata(handles.figure_JLabel,handles);

  % update plot
  if doupdateplot,
    UpdatePlots(handles,'axes',i);
  end
  
  % TODO: update timeline zoom
  for i = 1:numel(handles.axes_timelines),
    zoom(handles.axes_timelines(i),'reset');
  end
  
  % out of bounds for labeling? then turn off labeling
  if (t < handles.t0_curr || t > handles.t1_curr),
    if handles.label_state > 0,
      set(handles.togglebutton_label_behaviors(handles.label_state),'Value',0);
    elseif handles.label_state < 0,
      set(handles.togglebutton_label_unknown,'Value',0);
    end
    handles.label_state = 0;
    set([handles.togglebutton_label_behaviors,handles.togglebutton_label_unknown],'Enable','off');
  else
    set([handles.togglebutton_label_behaviors,handles.togglebutton_label_unknown],'Enable','on');
  end

  
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

% --------------------------------------------------------------------
function menu_file_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uiputfile('*.mat','Save classifier',handles.data.classifierfilename);
if ~ischar(filename),
  return;
end
handles.data.classifierfilename = fullfile(pathname,filename);
SetStatus(handles,sprintf('Saving classifier to %s',handles.data.classifierfilename));
handles.data.SaveClassifier();
handles.data.SaveLabels();
ClearStatus(handles);

% --------------------------------------------------------------------
function menu_file_exit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
  defaultpath = fileparts(handles.configfilename);
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
    
  try
    handles.configparams = ReadXMLParams(handles.configfilename);
  catch ME,
    uiwait(warndlg(sprintf('Error reading configuration from file %s: %s',handles.configfilename,getReport(ME)),'Error reading config file'));
    havefilename = false;
    continue;
  end
  
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
    try
      handles.labelcolors = eval(sprintf('%s(%d)',cm,handles.data.nbehaviors));
    catch ME,
      uiwait(warndlg(sprintf('Error using label colormap from config file: %s',getReport(ME)),'Error parsing config label colors'));
      handles.labelcolors = lines(handles.data.nbehaviors);
    end
  end
end
handles.labelunknowncolor = [0,0,0];

handles.correctcolor = [0,.7,0];
handles.incorrectcolor = [.7,0,0];

% create buttons for each label
handles = CreateLabelButtons(handles);

% maximum distance squared in fraction of axis to change frames when
% clicking on preview window
handles.max_click_dist_preview = .025^2;

% zoom state
handles.zoom_fly = true;
handles.zoom_fly_radius = nan;
set(handles.menu_view_zoom_in_on_fly,'Checked','on');

% initialize labels for navigation
SetJumpGoMenuLabels(handles)

% which behavior we seek to, initialize to all except unknown
handles.seek_behaviors_go = 1:handles.data.nbehaviors;

% label shortcuts
if numel(handles.label_shortcuts) ~= handles.data.nbehaviors + 1,
  handles.label_shortcuts = cellstr(num2str((0:handles.data.nbehaviors)'))';
end

function SetJumpGoMenuLabels(handles)

set(handles.menu_go_forward_X_frames,'Label',sprintf('Forward %d frames',handles.nframes_jump_go));
set(handles.menu_go_back_X_frames,'Label',sprintf('Back %d frames',handles.nframes_jump_go));

% create buttons for each label
function handles = CreateLabelButtons(handles)

% get positions of stuff
set(handles.panel_labelbuttons,'Units','pixels');
panel_pos = get(handles.panel_labelbuttons,'Position');
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
new_panel_height = 2*out_border_y + (handles.data.nbehaviors+1)*button_height + ...
  handles.data.nbehaviors*in_border_y;
% update panel position
panel_top = panel_pos(2)+panel_pos(4);
new_panel_pos = [panel_pos(1),panel_top-new_panel_height,panel_pos(3),new_panel_height];
set(handles.panel_labelbuttons,'Position',new_panel_pos);

% move unknown button to the bottom
new_unknown_button_pos = [unknown_button_pos(1),out_border_y,unknown_button_pos(3),button_height];
set(handles.togglebutton_label_unknown,'Position',new_unknown_button_pos);

% list of buttons
handles.togglebutton_label_behaviors = nan(1,handles.data.nbehaviors);

% update first button
new_button1_pos = [out_border_x,new_panel_height-out_border_y-button_height,button_width,button_height];
set(handles.togglebutton_label_behavior1,...
  'String',sprintf('Label %s',handles.data.labelnames{1}),...
  'ForegroundColor','w','Units','pixels','FontUnits','pixels','FontSize',14,...
  'FontWeight','bold','BackgroundColor',handles.labelcolors(1,:),...
  'Position',new_button1_pos,...
  'UserData',1);
handles.togglebutton_label_behaviors(1) = handles.togglebutton_label_behavior1;

% create the rest of the buttons
for i = 2:handles.data.nbehaviors,
  pos = [out_border_x,new_panel_height-out_border_y-button_height*i-in_border_y*(i-1),...
    button_width,button_height];
  handles.togglebutton_label_behaviors(i) = ...
    uicontrol('Style','togglebutton','String',sprintf('Label %s',handles.data.labelnames{i}),...
    'ForegroundColor','w','Units','pixels','FontUnits','pixels','FontSize',14,...
    'FontWeight','bold','BackgroundColor',handles.labelcolors(i,:),...
    'Position',pos,...
    'Callback',get(handles.togglebutton_label_behavior1,'Callback'),...
    'Parent',handles.panel_labelbuttons,...
    'Tag',sprintf('togglebutton_label_behavior%d',i),...
    'UserData',i);
end

% set props for unknown button
set(handles.togglebutton_label_unknown,...
  'String','Label Unknown',...
  'ForegroundColor','w','Units','pixels','FontUnits','pixels','FontSize',14,...
  'FontWeight','bold','BackgroundColor',handles.labelunknowncolor,...
  'UserData',-1);

function EnableGUI(handles)

% these controls require a movie to currently be open
h = [handles.menu_view_timeline_options,...
  handles.togglebutton_label_behaviors(:)',...
  handles.togglebutton_label_unknown,...
  handles.menu_view_zoom_in_on_fly];
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
handles.rcfilename = fullfile(fileparts(which('JLabel')),'.JLabelrc.mat');
handles.rc = struct;
if exist(handles.rcfilename,'file'),
  try
    handles.rc = load(handles.rcfilename);
  catch ME,
    warning('Error loading rc file %s: %s',handles.rcfilename,getReport(ME));
  end
end
try
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
    
catch ME,
  warning('Error loading RC file: %s',getReport(ME));  
end

function handles = SaveRC(handles)

try
  if ~isfield(handles,'rcfilename'),
    handles.rcfilename = fullfile(fileparts(which('JLabel')),'.JLabelrc.mat');
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
  
  save(handles.rcfilename,'-struct','rc');
catch ME,
  warning('Error saving RC file: %s',getReport(ME));
end


% --- Executes when user attempts to close figure_JLabel.
function figure_JLabel_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_JLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%delete(hObject);
if isfield(handles,'movie_fid') && ~isempty(handles.movie_fid) && ...
    ~isempty(fopen(handles.movie_fid)),
  fclose(handles.movie_fid);
  handles.movie_fid = [];
end
try
  % turn off zooming
  zoom(handles.figure_JLabel,'off');
catch %#ok<CTCH>
end
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
behaviori = get(hObject,'UserData');
if get(hObject,'Value'),
  % toggle on
  handles.label_state = behaviori;
  set(handles.togglebutton_label_behaviors(behaviori),'String',sprintf('*Label %s*',handles.data.labelnames{behaviori}));
  
  % set everything else to off
  for j = 1:handles.data.nbehaviors,
    if j == behaviori,
      continue;
    end
    set(handles.togglebutton_label_behaviors(j),'Value',0,'String',sprintf('Label %s',handles.data.labelnames{j}));
  end
  set(handles.togglebutton_label_unknown,'Value',0,'String','Label Unknown');
  
  % set the current frame to be labeled
  handles.lastframe_labeled = [];
  handles = SetLabelPlot(handles,min(handles.t1_curr,max(handles.t0_curr,handles.ts(1))),behaviori);
  UpdatePlots(handles,'refreshim',false,'refreshtrx',false,'refreshflies',false);
else
  handles.label_state = 0;
  %handles.data.StoreLabels();
  set(handles.togglebutton_label_behaviors(behaviori),'String',sprintf('Label %s',handles.data.labelnames{behaviori}));
end

guidata(hObject,handles);

function handles = SetLabelPlot(handles,t,behaviori)

if behaviori == 0,
  return;
end

if t == handles.lastframe_labeled,
  warning('This should never happen');
  keyboard;
end

if isempty(handles.lastframe_labeled),
  t0 = t;
  t1 = t;
  t2 = min(t+1,handles.t1_curr);
else
  if t < handles.lastframe_labeled,
    t0 = t;
    t1 = handles.lastframe_labeled-1;
    t2 = handles.lastframe_labeled;
  elseif t > handles.lastframe_labeled,
    t0 = handles.lastframe_labeled+1;
    t1 = t;
    t2 = min(t+1,handles.t1_curr);
  end
end

% WARNING: this function directly accesses handles.data.labelidx, trx make sure
% that we've preloaded the right experiment and flies. 
% REMOVED!
% if handles.expi ~= handles.data.expi || ~all(handles.flies == handles.data.flies),
%   handles.data.Preload(handles.expi,handles.flies);
% end

handles.labels_plot.x(t0+handles.labels_plot_off:t1+handles.labels_plot_off,:,:) = nan;
handles.labels_plot.y(t0+handles.labels_plot_off:t1+handles.labels_plot_off,:,:) = nan;
% handles.data.labelidx(t0+handles.labels_plot_off:t1+handles.labels_plot_off) = 0;

for channel = 1:3,
  handles.labels_plot.im(1,t0+handles.labels_plot_off:t1+handles.labels_plot_off,channel) = handles.labelunknowncolor(channel);
end
handles.data.SetLabel(handles.expi,handles.flies,t0:t1,behaviori);
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
    k2 = t2+handles.labels_plot_off;
    handles.labels_plot.x(k0:k2,behaviori,l) = ...
      handles.data.GetTrxX1(handles.expi,handles.flies(l),t0:t2);
    handles.labels_plot.y(k0:k2,behaviori,l) = ...
      handles.data.GetTrxY1(handles.expi,handles.flies(l),t0:t2);

%     handles.labels_plot.x(k0:k2,behaviori,l) = ...
%       handles.data.trx(handles.flies(l)).x(j0:j2);
%     handles.labels_plot.y(k0:k2,behaviori,l) = ...
%       handles.data.trx(handles.flies(l)).y(j0:j2);
  end
end

% isstart
if t0 == handles.t0_curr,
  handles.labels_plot.isstart(t+handles.labels_plot_off) = behaviori ~= 0;
end
t00 = max(handles.t0_curr+1,t0);
off0 = t00+handles.labels_plot_off;
off1 = t2+handles.labels_plot_off;
% handles.labels_plot.isstart(off0:off1) = ...
%   handles.data.labelidx(off0:off1)~=0 & ...
%   handles.data.labelidx(off0-1:off1-1)~=handles.data.labelidx(off0:off1);
handles.labels_plot.isstart(off0:off1) = ...
  handles.data.IsLabelStart(handles.expi,handles.flies,t00:t2);

handles = UpdateErrors(handles);

handles.lastframe_labeled = t;

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
  % set everything else to off
  for j = 1:handles.data.nbehaviors,
    set(handles.togglebutton_label_behaviors(j),'Value',0,'String',sprintf('Label %s',handles.data.labelnames{j}));
  end
  % set the current frame to be labeled
  handles.lastframe_labeled = [];
  handles = SetLabelPlot(handles,min(handles.t1_curr,max(handles.t0_curr,handles.ts(1))),-1);
  UpdatePlots(handles,'refreshim',false,'refreshtrx',false,'refreshflies',false);
  set(handles.togglebutton_label_unknown,'String','*Label Unknown*');
else
  handles.label_state = 0;
  %handles.data.StoreLabels();
  set(handles.togglebutton_label_unknown,'String','Label Unknown');
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
nprev = 25;
npost = 25;
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
  [mindcurr,k] = min( ((handles.data.GetTrxX1(handles.expi,fly,t0:t1)-xclick)/dx).^2 + ...
    ((handles.data.GetTrxY1(handles.expi,fly,t0:t1)-yclick)/dy).^2 );
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
res = questdlg({sprintf('Switch to fly %d?',fly),...
  sprintf('Trajectory length = %d',endframe-firstframe+1),...
  sprintf('First frame = %d',firstframe),...
  sprintf('N. bouts labeled: %d',nbouts)},...
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
t = min(handles.nframes,max(1,round(pt(1,1))));
% TODO: which axes?
SetCurrentFrame(handles,1,t,hObject);


% --- Executes on button press in pushbutton_train.
function pushbutton_train_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_train (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store the current labels to windowdata_labeled
handles.data.StoreLabels();
handles.data.Train();
% predict for current window
handles = UpdatePrediction(handles);
guidata(hObject,handles);

function handles = UpdatePrediction(handles)

% update prediction for currently shown timeline
% TODO: make this work for multiple axes
t0 = max(handles.t0_curr,floor(handles.ts(1)-handles.timeline_nframes/2));
t1 = min(handles.t1_curr,ceil(handles.ts(1)+handles.timeline_nframes/2));
handles.data.Predict(handles.expi,handles.flies,t0:t1);
handles = UpdateTimelineIms(handles);
UpdatePlots(handles,'refreshim',false,'refreshflies',...
  false,'refreshtrx',false,'refreshlabels',false,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_xlim',false);


function handles = UpdateErrors(handles)

% update prediction for currently shown timeline
% TODO: make this work for multiple axes
handles.data.UpdateErrorIdx();
handles = UpdateTimelineIms(handles);
UpdatePlots(handles,'refreshim',false,'refreshflies',...
  false,'refreshtrx',false,'refreshlabels',false,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_xlim',false);


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
function menu_file_load_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_load (see GCBO)
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
function menu_go_switch_experiment_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_switch_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = cell(1,handles.data.nexps);
for i = 1:handles.data.nexps,
  if i == handles.expi,
    s{i} = sprintf('%s, N flies: %d, OPEN NOW',...
      handles.data.expnames{i},handles.data.nflies_per_exp(i));
  else
    s{i} = sprintf('%s, N flies: %d, N flies labeled: %d, N bouts labeled: %d, last labeled: %s',...
      handles.data.expnames{i},handles.data.nflies_per_exp(i),...
      handles.data.labelstats(i).nflies_labeled,...
      handles.data.labelstats(i).nbouts_labeled,...
      handles.data.labelstats(i).datestr);
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
  
  % TODO: this function directly accesses handles.data.labels, abstract
  % this
  [ism,j] = ismember(fly,handles.data.labels(handles.expi).flies,'rows');
  if ism,
    nbouts = numel(handles.data.labels(handles.expi).t0s{j});
  else
    nbouts = 0;
  end
  
  if fly == handles.flies(1),
    s{fly} = sprintf('Target %d, CURRENTLY SELECTED',fly);
  else
    endframe = handles.data.endframes_per_exp{handles.expi}(fly);
    firstframe = handles.data.firstframes_per_exp{handles.expi}(fly);
    s{fly} = sprintf('Target %d, Trajectory length %d, First frame %d, N bouts labeled %d',...
      fly,endframe-firstframe+1,...
      firstframe,...
      nbouts);
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

handles.zoom_fly = ~handles.zoom_fly;
if handles.zoom_fly,
  set(hObject,'Checked','on');
else
  set(hObject,'Checked','off');
end
if handles.zoom_fly,
  ZoomInOnFlies(handles);
end

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
    xs(j) = handles.data.GetTrxX1(handles.expi,handles.flies(j),handles.ts(i));
    ys(j) = handles.data.GetTrxY1(handles.expi,handles.flies(j),handles.ts(i));
    %xs(j) = handles.data.trx(handles.flies(j)).x(inds(j));
    %ys(j) = handles.data.trx(handles.flies(j)).y(inds(j));
  end
  if ~all(isnan(xs)) && ~all(isnan(ys)),
    xlim = [max([.5,xs-handles.zoom_fly_radius]),min([handles.movie_width+.5,xs+handles.zoom_fly_radius])];
    ylim = [max([.5,ys-handles.zoom_fly_radius]),min([handles.movie_height+.5,ys+handles.zoom_fly_radius])];
    set(handles.axes_previews(i),'XLim',xlim,'YLim',ylim);
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
    else
      menu_go_previous_frame_Callback(hObject, eventdata, handles);
    end
     
  case 'rightarrow',
    if strcmpi(eventdata.Modifier,'control'),
      menu_go_next_bout_start_Callback(hObject,eventdata,handles);
    else
      menu_go_next_frame_Callback(hObject, eventdata, handles);
    end
  
  case 'uparrow',
    menu_go_back_X_frames_Callback(hObject, eventdata, handles);
    
  case 'downarrow',
    menu_go_forward_X_frames_Callback(hObject, eventdata, handles);
    
  case handles.label_shortcuts,
    behaviori = find(strcmp(eventdata.Key,handles.label_shortcuts),1);
    behaviori = behaviori - 1;
    if behaviori == 0,
      set(handles.togglebutton_label_unknown,'Value',get(handles.togglebutton_label_unknown,'Value')==0);
      togglebutton_label_unknown_Callback(handles.togglebutton_label_unknown, eventdata, handles);
    else
      set(handles.togglebutton_label_behaviors(behaviori),'Value',...
        get(handles.togglebutton_label_behaviors(behaviori),'Value')==0);
      togglebutton_label_behavior1_Callback(handles.togglebutton_label_behaviors(behaviori), eventdata, handles);
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
t = min(max(1,handles.ts(axesi)+1),handles.nframes);
% set current frame
SetCurrentFrame(handles,axesi,t,hObject);

% --------------------------------------------------------------------
function menu_go_previous_frame_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_previous_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;
t = min(max(1,handles.ts(axesi)-1),handles.nframes);
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
t = min(max(1,handles.ts(axesi)+handles.nframes_jump_go),handles.nframes);
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
t = min(max(1,handles.ts(axesi)-handles.nframes_jump_go),handles.nframes);
% set current frame
SetCurrentFrame(handles,axesi,t,hObject);


% --------------------------------------------------------------------
function menu_go_next_bout_start_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_next_bout_start (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;
if handles.ts(axesi) >= handles.t1_curr,
  return;
end
labelidx = handles.data.GetLabelIdx(handles.expi,handles.flies,handles.ts(axesi),handles.t1_curr);
j = find(labelidx ~= labelidx(1),1);
if isempty(j),
  return;
end
k = find(ismember(labelidx(j:end),handles.seek_behaviors_go),1);
if isempty(k),
  return;
end
t = handles.ts(axesi) + j - 1 + k - 1;
SetCurrentFrame(handles,axesi,t,hObject);

% --------------------------------------------------------------------
function menu_go_previous_bout_end_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_previous_bout_end (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes
axesi = 1;
if handles.t0_curr >= handles.ts(axesi),
  return;
end
labelidx = handles.data.GetLabelIdx(handles.expi,handles.flies,handles.t0_curr,handles.ts(axesi));
j = find(labelidx ~= labelidx(end),1,'last');
if isempty(j),
  return;
end
k = find(ismember(labelidx(1:j),handles.seek_behaviors_go),1,'last');
if isempty(k),
  return;
end
t = handles.t0_curr + k - 1;
SetCurrentFrame(handles,axesi,t,hObject);


% --------------------------------------------------------------------
function menu_go_navigation_preferences_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_navigation_preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'figure_NavigationPreferences') && ishandle(handles.figure_NavigationPreferences),
  figure(handles.figure_NavigationPreferences);
else
  handles.figure_NavigationPreferences = NavigationPreferences(handles.figure_JLabel);
  guidata(hObject,handles);
end


% --------------------------------------------------------------------
function menu_edit_label_shortcuts_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_label_shortcuts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


prompts = [{'Unknown'},handles.data.labelnames];
sh = inputdlg(prompts,'Label Shortcuts',1,handles.label_shortcuts);
handles.label_shortcuts = sh;
guidata(hObject,handles);