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

% Last Modified by GUIDE v2.5 10-Oct-2011 10:42:16

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

% get an index for this window number
handles.windowi = numel(parent_handles.bookmark_windows)+1;
parent_handles.bookmark_windows(end+1) = hObject;
guidata(handles.parent,parent_handles);

% optional arguments
[handles.clips,...
  handles.nclips_per_page,...
  handles.page] = myparse(varargin(2:end),'clips',[],...
  'nclips_per_page',12,...
  'page',1);

% Choose default command line output for BookmarkedClips
handles.output = hObject;

handles = InitializeState(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes BookmarkedClips wait for user response (see UIRESUME)
% uiwait(handles.figure_BookmarkedClips);

function handles = InitializeState(handles)

% number of clips
handles.nclips = numel(handles.clips);

handles.npages = ceil(handles.nclips/handles.nclips_per_page);

% get relative locations of stuffs
handles = GetGUIPositions(handles);


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
  notes_pos(3),... %b
  nan,... %t
  nan,... %w
  notes_pos(3)]; %h
% info
handles.guipos.info = [info_pos(1),... %l
  panel_pos(3)-info_pos(1)-info_pos(3),... %r
  info_pos(2)-notes_pos(2)-notes_pos(4),... %b
  nan,... %t
  nan,... %w
  info_pos(3)]; %h
% slider
handles.guipos.slider = [slider_pos(1),... %l
  nan,... %r
  slider_pos(2)-info_pos(2)-info_pos(4),... %b
  nan,... %t
  nan,... %w
  slider_pos(3)]; %h
% edit frame number
handles.guipos.edit_frame = [edit_pos(1)-slider_pos(1)-slider_pos(3),... %l
  nan,... %r
  nan,... %b
  nan,... %t
  edit_pos(3),... %w
  edit_pos(4)]; %h
% play button
handles.guipos.play = [play_pos(1)-edit_pos(1)-edit_pos(3),... %l
  panel_pos(3)-play_pos(1)-play_pos(3),... %r
  nan,... %b
  nan,... %t
  play_pos(3),... %w
  play_pos(4)]; %h
handles.guipos.axes = [axes_pos(1),... %l
  panel_pos(3)-axes_pos(1)-axes_pos(3),... %r
  axes_pos(2)-(slider_pos(1)+slider_pos(3)),... %b
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


% --- Outputs from this function are returned to the command line.
function varargout = BookmarkedClips_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_framenumber_clip1_Callback(hObject, eventdata, handles)
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

panel_pos = get(handles.panel_clip1,'Position');
fig_pos = get(handles.figure_BookmarkedClips,'Position');
page_pos = get(handles.popupmenu_show_page,'Position');

% choose the best layout for the panels
[handles.nr,handles.nc,handles.nclips_curr,handles.panel_width,handles.panel_height] = ...
  ChooseClipLayout(handles);

% page
new_page_edit_pos = [fig_pos(3)-(handles.guipos.page_edit(1)+handles.guipos.page_edit(5)),...
  handles.guipos.page_edit(2),handles.guipos.page_edit(5),handles.guipos.page_edit(6)];
set(handles.popupmenu_show_page,'Position',new_page_edit_pos);
new_page_text_pos = [new_page_edit_pos(1)-(handles.guipos.page_text(1)+handles.guipos.page_text(5)),...
  handles.guipos.page_text(2),handles.guipos.page_text(5),handles.guipos.page_text(6)];
set(handles.text_show_page,'Position',new_page_text_pos);

% panels
for i = 1:handles.nclips_curr,
end
for i = handles.nclips_curr+1:numel(handles.panel_clips),
  set(handles.panel_clips(i),'Visible','off');
end