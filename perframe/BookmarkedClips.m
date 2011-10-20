function varargout = BookmarkedClips(varargin)
% BOOKMARKEDCLIPS MATLAB code for BookmarkedClips.fig
%      BOOKMARKEDCLIPS, by itself, creates a new BOOKMARKEDCLIPS or raises the existing
%      singleton*.
%
%      H = BOOKMARKEDCLIPS returns the handle to a new BOOKMARKEDCLIPS or the handle to
%      the existing singleton*.
%
%      BOOKMARKEDCLIPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BOOKMARKEDCLIPS.M with the given input arguments.
%
%      BOOKMARKEDCLIPS('Property','Value',...) creates a new BOOKMARKEDCLIPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BookmarkedClips_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BookmarkedClips_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BookmarkedClips

% Last Modified by GUIDE v2.5 14-Oct-2011 05:31:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BookmarkedClips_OpeningFcn, ...
                   'gui_OutputFcn',  @BookmarkedClips_OutputFcn, ...
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


% --- Executes just before BookmarkedClips is made visible.
function BookmarkedClips_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to BookmarkedClips (see VARARGIN)

% first input is the JLabel GUI
handles.parent = varargin{1};
parent_handles = guidata(handles.parent);

% second input is the JLabelData
handles.data = varargin{2};

% get an index for this window number
handles.windowi = numel(parent_handles.bookmark_windows)+1;
parent_handles.bookmark_windows(end+1) = hObject;
guidata(handles.parent,parent_handles);

% optional arguments
[handles.clips,...
  handles.nclips_per_page,...
  handles.page] = myparse(varargin(3:end),'clips',[],...
  'nclips_per_page',12,...
  'page',1);

% Choose default command line output for BookmarkedClips
handles.output = hObject;

% colors for each fly
handles.fly_colors = parent_handles.fly_colors;

handles = InitializeState(handles);

handles = InitializePlots(handles);

% which clips are shown
clip0 = (handles.page-1)*handles.nclips_per_page+1;
clip1 = min(handles.page*handles.nclips_per_page,handles.nclips);
for clipi = clip0:clip1,
  handles = SetClipToPanel(handles,clipi-clip0+1,clipi);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BookmarkedClips wait for user response (see UIRESUME)
% uiwait(handles.figure_BookmarkedClips);

function handles = InitializeState(handles)

% number of clips
handles.nclips = numel(handles.clips);

% number of flies plotted
handles.nflies_curr = max(cellfun(@numel,{handles.clips.flies}));

% how many pages of clips
handles.npages = ceil(handles.nclips/handles.nclips_per_page);

% which clip is stored in each panel
handles.panel2clip = nan(1,handles.nclips_per_page);

% which frame we are currently showing
handles.ts = [handles.clips.t0];

% whether we are playing each panel
handles.isplaying = false(1,handles.nclips_per_page);

% readframe functions
handles.readframes = {};
handles.panel2readframei = nan(1,handles.nclips_per_page);
handles.readframe_expis = [];
handles.readframe_fids = [];
handles.movie_heights = [];
handles.movie_widths = [];

% get relative locations of stuffs
handles = GetGUIPositions(handles);

function handles = InitializePlots(handles)

handles.play_color = [.2,.4,0];
handles.stop_color = [.5,0,0];

set(handles.figure_BookmarkedClips,'Colormap',gray(256));

handles.panel_clips = handles.panel_clip1;
v = get(handles.panel_clip1);
fnscopy = {'Units','Position','BackgroundColor','ForegroundColor','FontUnits','ResizeFcn'};
fnscopy_chil = {'Units','Position','BackgroundColor','ForegroundColor',...
  'FontUnits','FontSize','FontWeight','HorizontalAlignment',...
  'Min','Max','String','Value','Callback','SliderStep','TooltipString','UIContextMenu',...
  'Visible'};
fnscopy_axes = {'Units','Position','Color','FontName','FontSize','FontUnits','FontWeight','XColor','YColor'};
  
hchil = v.Children;
for j = 1:numel(hchil),
  tag = get(hchil(j),'Tag');
  fn = [tag(1:end-1),'s'];
  handles.(fn) = handles.(tag);
end
for i = 2:handles.nclips_per_page,
  handles.panel_clips(i) = ...
    uipanel('Parent',handles.figure_BookmarkedClips,...
    'Title','',...
    'Visible','off',...
    'Tag',sprintf('panel_clip%d',i));
  for j = 1:numel(fnscopy),
    set(handles.panel_clips(i),fnscopy{j},v.(fnscopy{j}));
  end
  for j = 1:numel(hchil),
    vchil = get(hchil(j));
    tag = vchil.Tag(1:end-1);
    fn = [tag,'s'];
    if strcmpi(vchil.Type,'uicontrol'),
      handles.(fn)(i) = uicontrol('Parent',handles.panel_clips(i),...
        'Style',vchil.Style,...
        'Tag',[tag,num2str(i)]);
      for k = 1:numel(fnscopy_chil),
        set(handles.(fn)(i),fnscopy_chil{k},vchil.(fnscopy_chil{k}));
      end
    elseif strcmpi(vchil.Type,'axes'),
      handles.(fn)(i) = axes('Parent',handles.panel_clips(i),...
        'Tag',[tag,num2str(i)]);
      for k = 1:numel(fnscopy_axes),
        set(handles.(fn)(i),fnscopy_axes{k},vchil.(fnscopy_axes{k}));
      end
      
    end
  end  
end
for i = 1:handles.nclips_per_page,
  hold(handles.axes_clips(i),'on');
  handles.himage_clips(i) = image(uint8(0),'Parent',handles.axes_clips(i));
  set(handles.axes_clips(i),'CLim',[0,255]);
  axis(handles.axes_clips(i),'image','ij');
  for j = 1:handles.nflies_curr,
    handles.hflies(j,i) = plot(handles.axes_clips(i),nan,nan,'-','color','r','linewidth',2);
  end
end

function handles = SetClipToPanel(handles,paneli,clipi)

if handles.panel2clip(paneli) == clipi,
  return;
end

% store the index
handles.panel2clip(paneli) = clipi;

% set the slider
nframes = handles.clips(clipi).t1-handles.clips(clipi).t0+1;
if handles.clips(clipi).t0 == handles.clips(clipi).t1,
  set(handles.slider_framenumbber_clips(paneli),'Visible','off');
else
  set(handles.slider_framenumber_clips(paneli),'Min',handles.clips(clipi).t0,...
    'Max',handles.clips(clipi).t1,...
    'Value',handles.ts(clipi),...
    'SliderStep',[1/(nframes-1),ceil(nframes/10)/(nframes-1)],...
    'Visible','on');
end

% set the edit frame number
set(handles.edit_framenumber_clips(paneli),'String',num2str(handles.ts(clipi)));

% set to not playing
handles.isplaying(paneli) = false;
set(handles.pushbutton_playstop_clips(paneli),'String','Play',...
  'BackgroundColor',handles.play_color);

% set info string
s = cell(1,3);
s{1} = handles.data.expnames{handles.clips(clipi).expi};
if numel(handles.clips(clipi).flies) > 1,
  s{2} = sprintf('Flies %s, ',mat2str(handles.clips(clipi).flies));
else
  s{2} = sprintf('Fly %d, ',handles.clips(clipi).flies);  
end
s{2} = [s{2},sprintf('Frames %d:%d',handles.clips(clipi).t0,handles.clips(clipi).t1)];
s{3} = sprintf('Label: %s',handles.clips(clipi).label);
set(handles.text_info_clips(paneli),'String',s);

% set notes
% TODO

% get readframe fcn
readframei = find(handles.clips(clipi).expi==handles.readframe_expis,1);
if isempty(readframei),
  moviefilename = handles.data.GetFile('movie',handles.clips(clipi).expi);
  [readframe,~,fid] = get_readframe_fcn(moviefilename);
  readframei = numel(handles.readframes)+1;
  handles.readframes{readframei} = readframe;
  handles.readframe_expis(readframei) = handles.clips(clipi).expi;
  handles.readframe_fids(readframei) = fid;
  im = handles.readframes{readframei}(handles.ts(clipi));
  handles.movie_heights(readframei) = size(im,1);
  handles.movie_widths(readframei) = size(im,2);
end
handles.panel2readframei(paneli) = readframei;

% initialize axes lims
if handles.clips(clipi).zoom_fly,
  ZoomInOnFlies(handles,paneli);
else
  set(handles.axes_clips(paneli),'XLim',handles.clips(clip).xlim,...
    'YLim',handles.clips(clip).ylim);
end

% colors
for i = 1:numel(handles.clips(clipi).flies),
  fly = handles.clips(clipi).flies(i);
  set(handles.hflies(i,paneli),'Visible','on','Color',handles.fly_colors(fly,:));
end
for i = numel(handles.clips(clipi).flies)+1:handles.nflies_curr,
  set(handles.hflies(i,paneli),'Visible','off');
end  

handles = UpdatePanels(handles,'panels',paneli);

function handles = UpdatePanels(handles,varargin)

[panelis,refreshim,refreshflies,refreshtrx,refreshlabels] = ...
  myparse(varargin,'panels',1:handles.nclips_per_page,...
  'refreshim',true,'refreshflies',true,'refreshtrx',true,'refreshlabels',true);

for paneli = panelis,
  
  clipi = handles.panel2clip(paneli);
  if isnan(clipi),
    set(handles.panel_clips(paneli),'Visible','off');
    continue;
  end
  
  t = handles.ts(clipi);
  i = t - handles.clips(clipi).t0 + 1;
  
  if refreshim,
    % read in current frame
    im = handles.readframes{handles.panel2readframei(paneli)}(t);    
    % update frame
    set(handles.himage_clips(paneli),'CData',im);
  end
  
  % update current position
  if refreshflies,
  %     if handles.ts(i) < handles.t0_curr || handles.ts(i) > handles.t1_curr,
  %       labelidx = [];
  %     else
  %       labelidx = handles.data.GetLabelIdx(handles.expi,handles.flies,handles.ts(i),handles.ts(i));
  %     end
  %     inbounds = handles.data.firstframes_per_exp{handles.expi} <= handles.ts(i) & ...
  %       handles.data.endframes_per_exp{handles.expi} >= handles.ts(i);
  %     set(handles.hflies(~inbounds,i),'XData',nan,'YData',nan);
  for j = 1:numel(handles.clips(clipi).flies),
    updatefly(handles.hflies(j,paneli),...
      handles.clips(clipi).trx(j).x(i),...
      handles.clips(clipi).trx(j).y(i),...
      handles.clips(clipi).trx(j).theta(i),...
      handles.clips(clipi).trx(j).a(i),...
      handles.clips(clipi).trx(j).b(i));
  end
  %     for fly = find(inbounds),
  %       t = handles.ts(i);
  %       [xcurr,ycurr,thetacurr,acurr,bcurr] = ...
  %         handles.data.GetTrxPos1(handles.expi,fly,t);
  %       updatefly(handles.hflies(fly,i),...
  %         xcurr,ycurr,thetacurr,acurr,bcurr);
  %       if ismember(fly,handles.flies),
  %         set(handles.hflies(fly,i),'LineWidth',5);
  %         if labelidx <= 0,
  %           set(handles.hflies(fly,i),'Color',handles.labelunknowncolor);
  %         else
  %           set(handles.hflies(fly,i),'Color',handles.labelcolors(labelidx,:));
  %         end
  %       else
  %         set(handles.hflies(fly,i),'LineWidth',2);
  %       end
  end
  
  if handles.clips(clipi).zoom_fly,
    ZoomInOnFlies(handles,paneli);
  end
% 
%   % update trx
%   % TODO: remove hard-coded nprev, npost
%   nprev = 25;
%   npost = 25;
%   if refreshtrx,
%     for j = 1:numel(handles.flies),
%       fly = handles.flies(j);
%       tmp = handles.ts(i);
%       t0 = handles.data.firstframes_per_exp{handles.expi}(fly);
%       t1 = handles.data.endframes_per_exp{handles.expi}(fly);
%       ts = max(t0,tmp-nprev):min(t1,tmp+npost);
%       set(handles.htrx(j,i),'XData',handles.data.GetTrxX1(handles.expi,fly,ts),...
%         'YData',handles.data.GetTrxY1(handles.expi,fly,ts));
%       %j0 = max(1,tmp-nprev);
%       %j1 = min(handles.data.trx(fly).nframes,tmp+npost);
%       %set(handles.htrx(j,i),'XData',handles.data.trx(fly).x(j0:j1),...
%       %  'YData',handles.data.trx(fly).y(j0:j1));
%       %trx = handles.data.GetTrx(handles.expi,fly,handles.ts(i)-nprev:handles.ts(i)+npost);
%       %set(handles.htrx(j,i),'XData',trx.x,'YData',trx.y);
%     end
%   end  
%   
%   % update labels plotted
%   if refreshlabels,
%     for k = 1:numel(handles.flies),
%       fly = handles.flies(k);
%       T0 = handles.data.firstframes_per_exp{handles.expi}(fly);
%       T1 = handles.data.endframes_per_exp{handles.expi}(fly);
% %       T0 = handles.data.GetTrxFirstFrame(handles.expi,fly);
% %       T1 = handles.data.GetTrxEndFrame(handles.expi,fly);
%       t0 = min(T1,max(T0,handles.ts(i)-nprev));
%       t1 = min(T1,max(T0,handles.ts(i)+npost));
%       for j = 1:handles.data.nbehaviors,
%         set(handles.hlabels(j),'XData',handles.labels_plot.x(handles.labels_plot_off+t0:handles.labels_plot_off+t1,j,k),...
%           'YData',handles.labels_plot.y(handles.labels_plot_off+t0:handles.labels_plot_off+t1,j,k));
%       end
%     end
%   end
  
end

function ZoomInOnFlies(handles,panelis)

if nargin < 2,
  panelis = 1:numel(handles.clips_per_page);
end

for paneli = panelis,

  if isnan(handles.panel2clip(paneli)),
    continue;
  end
  clipi = handles.panel2clip(paneli);

  i = handles.ts(clipi) - handles.clips(clipi).t0 + 1;

  xs = nan(1,numel(handles.clips(clipi).flies));
  ys = nan(1,numel(handles.clips(clipi).flies));

  for j = 1:numel(handles.clips(clipi).flies),
    xs(j) = handles.clips(clipi).trx(j).x(i);
    ys(j) = handles.clips(clipi).trx(j).y(i);
  end
  if ~all(isnan(xs)) && ~all(isnan(ys)),
    w = handles.movie_widths(handles.panel2readframei(paneli));
    h = handles.movie_heights(handles.panel2readframei(paneli));
    xlim = [max([.5,xs-handles.clips(clipi).zoom_fly_radius(1)]),min([w+.5,xs+handles.clips(clipi).zoom_fly_radius(1)])];
    ylim = [max([.5,ys-handles.clips(clipi).zoom_fly_radius(2)]),min([h+.5,ys+handles.clips(clipi).zoom_fly_radius(2)])];
    set(handles.axes_clips(paneli),'XLim',xlim,'YLim',ylim);
  end
end

function [nr,nc,n,w,h] = ChooseClipLayout(handles)

fig_pos = get(handles.figure_BookmarkedClips,'Position');
figw = fig_pos(3) - handles.guipos.panel_border;
figh = fig_pos(4) - handles.guipos.page_edit(3) - handles.guipos.page_edit(6) - handles.guipos.panel_border;

rho = handles.guipos.ideal_panel_ratio;
if handles.page == handles.npages,
  n = rem(handles.nclips,handles.nclips_per_page);
  if n == 0,
    n = handles.nclips_per_page;
  end
else
  n = handles.nclips_per_page;
end
for i = 1:n,
  set(handles.panel_clips(i),'Visible','on');
end
for i = n+1:handles.nclips_per_page,
  set(handles.panel_clips(i),'Visible','off');
end

% approx ratio of nc/nr approx figw / (figh * rho)
% nr * nr * figw / (figh * rho) >= n
nr = round(sqrt(n * figh * rho / figw));
nrstry = [nr-1,nr,nr+1];
nrstry = nrstry(nrstry >= 1 & nrstry <= n);
ws = nan(1,numel(nrstry));
hs = nan(1,numel(nrstry));
for i = 1:numel(nrstry),
  nr = nrstry(i);
  nc = ceil(n/nr);
  ws(i) = figw / nc - handles.guipos.panel_border;
  hs(i) = figh / nr - handles.guipos.panel_border;
end
truews = min(ws,hs*rho);
truehs = min(hs,ws/rho);
i = argmax(truews.*truehs);
nr = nrstry(i);
nc = ceil(n/nr);
w = ws(i);
h = hs(i);


function handles = GetGUIPositions(handles)

panel_pos = get(handles.panel_clip1,'Position');
fig_pos = get(handles.figure_BookmarkedClips,'Position');
page_edit_pos = get(handles.popupmenu_show_page,'Position');
page_text_pos = get(handles.text_show_page,'Position');
info_pos = get(handles.text_info_clip1,'Position');
slider_pos = get(handles.slider_framenumber_clip1,'Position');
edit_pos = get(handles.edit_framenumber_clip1,'Position');
play_pos = get(handles.pushbutton_playstop_clip1,'Position');
axes_pos = get(handles.axes_clip1,'Position');
notes_pos = get(handles.edit_notes_clip1,'Position');

% use the left border to set all panel borders
handles.guipos.panel_border = panel_pos(1);

% page edit
handles.guipos.page_edit = [nan,... %l
  fig_pos(3)-(page_edit_pos(1)+page_edit_pos(3)),... %r
  page_edit_pos(2),... %b
  nan,... %t
  page_edit_pos(3),... %w
  page_edit_pos(4)]; %h

% page text
handles.guipos.page_text = [nan,... %l
  page_edit_pos(1)-(page_text_pos(1)+page_text_pos(3)),... %r
  page_text_pos(2),... %b
  nan,... %t
  page_text_pos(3),... %w
  page_text_pos(4)]; %h

% notes (l, r, b, t, w, h)
handles.guipos.notes = [notes_pos(1),... %l
  panel_pos(3)-notes_pos(1)-notes_pos(3),... %r
  notes_pos(2),... %b
  nan,... %t
  nan,... %w
  notes_pos(4)]; %h
% info
handles.guipos.info = [info_pos(1),... %l
  panel_pos(3)-info_pos(1)-info_pos(3),... %r
  info_pos(2)-notes_pos(2)-notes_pos(4),... %b
  nan,... %t
  nan,... %w
  info_pos(4)]; %h
% slider
handles.guipos.slider = [slider_pos(1),... %l
  edit_pos(1)-slider_pos(1)-slider_pos(3),... %r
  slider_pos(2)-info_pos(2)-info_pos(4),... %b
  nan,... %t
  nan,... %w
  slider_pos(4)]; %h
% edit frame number
handles.guipos.edit_frame = [edit_pos(1)-slider_pos(1)-slider_pos(3),... %l
  play_pos(1)-edit_pos(1)-edit_pos(3),... %r
  edit_pos(2)-info_pos(2)-info_pos(4),... %b
  nan,... %t
  edit_pos(3),... %w
  edit_pos(4)]; %h
% play button
handles.guipos.play = [play_pos(1)-edit_pos(1)-edit_pos(3),... %l
  panel_pos(3)-play_pos(1)-play_pos(3),... %r
  play_pos(2)-info_pos(2)-info_pos(4),... %b
  nan,... %t
  play_pos(3),... %w
  play_pos(4)]; %h
handles.guipos.axes = [axes_pos(1),... %l
  panel_pos(3)-axes_pos(1)-axes_pos(3),... %r
  axes_pos(2)-(slider_pos(2)+slider_pos(4)),... %b
  panel_pos(4)-(axes_pos(2)+axes_pos(4)),... %t
  nan,... %w
  nan]; %h

min_slider_width = 100;
min_axes_height = 50;
handles.guipos.min_panel_width = ...
  handles.guipos.slider(1) + min_slider_width + ...
  handles.guipos.edit_frame(1) + handles.guipos.edit_frame(5) + ...
  handles.guipos.play(1) + handles.guipos.play(5);
handles.guipos.min_panel_height = ...
  handles.guipos.notes(3) + handles.guipos.notes(6) + ...
  handles.guipos.info(3) + handles.guipos.info(6) + ...
  handles.guipos.slider(3) + handles.guipos.slider(6) + ...
  handles.guipos.axes(3) + min_axes_height;
handles.guipos.ideal_panel_ratio = panel_pos(3)/panel_pos(4);
handles.guipos.fig_bottom_border = ...
  handles.guipos.page_edit(3) + handles.guipos.page_edit(6) + ...
  handles.guipos.panel_border;

% --- Outputs from this function are returned to the command line.
function varargout = BookmarkedClips_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_framenumber_clip1_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
% hObject    handle to slider_framenumber_clip1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_framenumber_clip1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_framenumber_clip1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit_framenumber_clip1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_framenumber_clip1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_framenumber_clip1 as text
%        str2double(get(hObject,'String')) returns contents of edit_framenumber_clip1 as a double


% --- Executes during object creation, after setting all properties.
function edit_framenumber_clip1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_framenumber_clip1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_playstop_clip1.
function pushbutton_playstop_clip1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_playstop_clip1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_notes_clip1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_notes_clip1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_notes_clip1 as text
%        str2double(get(hObject,'String')) returns contents of edit_notes_clip1 as a double


% --- Executes during object creation, after setting all properties.
function edit_notes_clip1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_notes_clip1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in popupmenu_show_page.
function popupmenu_show_page_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_show_page (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_show_page contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_show_page


% --- Executes during object creation, after setting all properties.
function popupmenu_show_page_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_show_page (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_new_window_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_new_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_view_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_view_clips_per_page_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_clips_per_page (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_view_sort_by_Callback(hObject, eventdata, handles)
% hObject    handle to menu_view_sort_by (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_close_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when figure_BookmarkedClips is resized.
function figure_BookmarkedClips_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure_BookmarkedClips (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fig_pos = get(handles.figure_BookmarkedClips,'Position');

% choose the best layout for the panels
[handles.nr,handles.nc,handles.nclips_curr,handles.panel_width,handles.panel_height] = ...
  ChooseClipLayout(handles);

% compute the minimum figure size
nr_curr = ceil(handles.nclips_curr/handles.nc);
if nr_curr == 1,
  nc_curr = handles.nclips_curr;
else
  nc_curr = handles.nc;
end
min_h = nr_curr * (handles.guipos.min_panel_height + handles.guipos.panel_border) + ...
  handles.guipos.fig_bottom_border;
min_w = nc_curr * (handles.guipos.min_panel_width + handles.guipos.panel_border) + ...
  handles.guipos.panel_border;

if fig_pos(3) < min_w || fig_pos(4) < min_h,
  fig_pos(3) = max(min_w,fig_pos(3));
  fig_pos(4) = max(min_h,fig_pos(4));
  % make sure that the top left corner is on the screen
  screen_pos = get(0,'ScreenSize');
  if fig_pos(2) + fig_pos(4) > screen_pos(4),
    fig_pos(2) = screen_pos(4) - fig_pos(4);
  end
  set(handles.figure_BookmarkedClips,'Position',fig_pos);
  return;
end

% page
new_page_edit_pos = [fig_pos(3)-(handles.guipos.page_edit(2)+handles.guipos.page_edit(5)),...
  handles.guipos.page_edit(3),handles.guipos.page_edit(5),handles.guipos.page_edit(6)];
set(handles.popupmenu_show_page,'Position',new_page_edit_pos);
new_page_text_pos = [new_page_edit_pos(1)-(handles.guipos.page_text(2)+handles.guipos.page_text(5)),...
  handles.guipos.page_text(3),handles.guipos.page_text(5),handles.guipos.page_text(6)];
set(handles.text_show_page,'Position',new_page_text_pos);

% set panel positions

% left border for all panels
off_x0 = handles.guipos.panel_border;
% top border for all panels
off_y0 = handles.guipos.panel_border;

for i = 1:handles.nclips_curr,
  % which panel is this
  [c,r] = ind2sub([handles.nc,handles.nr],i);
  % compute position from left and top
  off_x = off_x0 + (c-1)*(handles.panel_width+handles.guipos.panel_border);
  off_y = off_y0 + (r-1)*(handles.panel_height+handles.guipos.panel_border);
  panel_pos = [off_x,fig_pos(4)-off_y-handles.panel_height,handles.panel_width,handles.panel_height];
  set(handles.panel_clips(i),'Position',panel_pos);
end


% --- Executes when panel_clip1 is resized.
function panel_clip1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to panel_clip1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

paneli = find(handles.panel_clips == hObject,1);
if isempty(paneli),
  warning('Could not find panel corresponding to %f',hObject);
  return;
end

panel_pos = get(hObject,'Position');

% update notes position
new_notes_pos = [handles.guipos.notes(1),...
  handles.guipos.notes(3),...
  panel_pos(3) - (handles.guipos.notes(1)+handles.guipos.notes(2)),...
  handles.guipos.notes(6)];
set(handles.edit_notes_clips(paneli),'Position',new_notes_pos);

% update info position
new_info_pos = [handles.guipos.info(1),...
  new_notes_pos(2)+new_notes_pos(4)+handles.guipos.info(3),...
  panel_pos(3) - (handles.guipos.info(1)+handles.guipos.info(2)),...
  handles.guipos.info(6)];
set(handles.text_info_clips(paneli),'Position',new_info_pos);

% play
new_play_pos = [panel_pos(3)-handles.guipos.play(2)-handles.guipos.play(5),...
  new_info_pos(2)+new_info_pos(4)+handles.guipos.play(3),...
  handles.guipos.play(5),...
  handles.guipos.play(6)];
set(handles.pushbutton_playstop_clips(paneli),'Position',new_play_pos);

% edit frame number
new_edit_frame_pos = [new_play_pos(1)-handles.guipos.edit_frame(2)-handles.guipos.edit_frame(5),...
  new_info_pos(2)+new_info_pos(4)+handles.guipos.edit_frame(3),...
  handles.guipos.edit_frame(5),...
  handles.guipos.edit_frame(6)];
set(handles.edit_framenumber_clips(paneli),'Position',new_edit_frame_pos);

% slider
new_slider_pos = [handles.guipos.slider(1),...
  new_info_pos(2)+new_info_pos(4)+handles.guipos.slider(3),...
  new_edit_frame_pos(1) - (handles.guipos.slider(1)+handles.guipos.slider(2)),...
  handles.guipos.slider(6)];
set(handles.slider_framenumber_clips(paneli),'Position',new_slider_pos);

% axes
b = new_slider_pos(2)+new_slider_pos(4)+handles.guipos.axes(3);
new_axes_pos = [handles.guipos.axes(1),...
  b,...
  panel_pos(3) - (handles.guipos.axes(1) + handles.guipos.axes(2)),...
  panel_pos(4) - b - handles.guipos.axes(4)];
set(handles.axes_clips(paneli),'Position',new_axes_pos);
