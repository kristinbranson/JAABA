function varargout = playfmf_interestpoints(varargin)
% PLAYFMF_INTERESTPOINTS M-file for playfmf_interestpoints.fig
%      PLAYFMF_INTERESTPOINTS, by itself, creates a new PLAYFMF_INTERESTPOINTS or raises the existing
%      singleton*.
%
%      H = PLAYFMF_INTERESTPOINTS returns the handle to a new PLAYFMF_INTERESTPOINTS or the handle to
%      the existing singleton*.
%
%      PLAYFMF_INTERESTPOINTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLAYFMF_INTERESTPOINTS.M with the given input arguments.
%
%      PLAYFMF_INTERESTPOINTS('Property','Value',...) creates a new PLAYFMF_INTERESTPOINTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before playfmf_interestpoints_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to playfmf_interestpoints_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help playfmf_interestpoints

% Last Modified by GUIDE v2.5 03-May-2023 10:38:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @playfmf_interestpoints_OpeningFcn, ...
                   'gui_OutputFcn',  @playfmf_interestpoints_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before playfmf_interestpoints is made visible.
function playfmf_interestpoints_OpeningFcn(hObject, eventdata, handles, varargin)
% Play movie with mouse-annotation.
% playfmf_interestpoints('moviefile',fname,'annDataObj',annDataObj)
%
% annDataObj should be a scalar HandleObj; use this optional argument when
% annotations are desired. The result of annotation will be stored in
% annDataObj.data.


handles.MAXSLIDERRES = 10^6;

[moviename,annDataObj,handles.ipnames,handles.doaskleftright,handles.ipsshow,startframe,figpos] = ...
  myparse(varargin,'moviename',[],'annDataObj',[],...
    'ipnames',{'food','mouth','perch'},'doaskleftright',true,...
    'ipsshow',[],'startframe',1,'figpos',[]);
if ~isempty(annDataObj)
  assert(isa(annDataObj,'HandleObj'),'Handle object for annotations must be supplied.');
end

handles.tfAnnotate = ~isempty(annDataObj);
handles.annDataObj = annDataObj;
handles.tfipsshow = ~isempty(handles.ipsshow);

assert(~(handles.tfAnnotate && handles.tfipsshow),...
  'Cannot simultaneously show existing annotations as well as create new annotations.');

handles.figpos = get(handles.figure1,'Position');

if ~handles.tfAnnotate,
  buttonnames = {'pbAnnotate','pbDone','pbCancel'};
  for i = 1:numel(buttonnames),
    buttonname = buttonnames{i};
    if isfield(handles,buttonname) && ishandle(handles.(buttonname)),
      set(handles.(buttonname),'Visible','off');
    end
  end
  axpos = get(handles.axes_Video,'Position');
  w = handles.figpos(3)-2*axpos(1);
  axpos(3) = w;
  set(handles.axes_Video,'Position',axpos);
end

% set up path
if isempty(which('myparse')),
  if exist('../misc','file'),
    addpath('../misc');
  end
  while isempty(which('myparse')),
    miscdir = uigetdir('.','Choose "misc" folder to add to path');
    try
      addpath(miscdir);
    catch %#ok<CTCH>
    end
  end
end
    
% load previous values
handles.rcfilename = '.playfmfrc.mat';
handles.previous_values = struct;
if exist(handles.rcfilename,'file')
  handles.previous_values = load(handles.rcfilename);
end

tfAutoOpen = false;
tfIsOpen = false;
if ~isempty(moviename)
  if iscell(moviename) && numel(moviename) == 2 && ...
      isa(moviename{1},'function_handle') && isnumeric(moviename{2}),
    handles.filename = '';
    handles.filedir = '';
    handles.filenamebase = '';
    handles.fileext = '';
    handles.readframe = moviename{1};
    handles.nframes = moviename{2};
    handles.headerinfo = struct;
    handles.fid = -1;
    tfAutoOpen = false;
    tfIsOpen = true;
  else
    handles.filename = moviename;
    [handles.filedir,handles.filenamebase,handles.fileext] = ...
      fileparts(handles.filename);
    if exist(handles.filename,'file')
      % user has explicitly supplied an existing movie file
      tfAutoOpen = true;
    else
      warning('playfmf:fileNotFound','Cannot find movie ''%s''.',handles.filename);
      handles.filename = '';
    end
  end
elseif isfield(handles.previous_values,'filename'),
  handles.filename = handles.previous_values.filename;
  [handles.filedir,handles.filenamebase,handles.fileext] = ...
    fileparts(handles.previous_values.filename);
  if ~exist(handles.filename,'file'),
    handles.filename = '';
  end
else
  handles.filename = '';
end

if isfield(handles.previous_values,'CompressionSettings'),
  handles.CompressionSettings = handles.previous_values.CompressionSettings;
else
  handles.CompressionSettings = struct;
end
if ~isfield(handles.CompressionSettings,'OutputFPS'),
  handles.CompressionSettings.OutputFPS = 30;
end
if ~isfield(handles.CompressionSettings,'Compression'),
  handles.CompressionSettings.Compression = 'None';
end
if ~isfield(handles.CompressionSettings,'Quality'),
  handles.CompressionSettings.Quality = 100;
end
handles.CompressionSettings.StartFrame = 1;
handles.CompressionSettings.EndFrame = inf;

if isfield(handles.previous_values,'MaxFPS'),
  handles.MaxFPS = handles.previous_values.MaxFPS;
else
  handles.MaxFPS = 0;
  handles.MinSPF = 0;
end

% set callback for slider motion
fcn = get(handles.slider_Frame,'Callback');

if verLessThan('matlab','8.4.0')
  handles.hslider_listener = handle.listener(handles.slider_Frame,...
    'ActionEvent',fcn);
  set(handles.slider_Frame,'Callback','');
else
  handles.hslider_listener = addlistener(handles.slider_Frame,...
    'ContinuousValueChange',fcn);
  set(handles.slider_Frame,'Callback','');
end

% open video
handles = open_fmf(handles,~tfAutoOpen,tfIsOpen,startframe);

% init annotation state
if handles.tfAnnotate
  handles.impts = [];
else
  set(handles.pbAnnotate,'Enable','off');
end

set(handles.pbDone,'Enable','off');

% set figure size
if ~isempty(figpos),
  set(handles.figure1,'Position',figpos);
  handles = figure1_ResizeFcn(handles.figure1,[],handles);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes playfmf_interestpoints wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function play(hObject)

global ISPLAYING;

handles = guidata(hObject);
set(hObject,'String','Stop','BackgroundColor',[.5,0,0]);
ISPLAYING = true;

fnext = handles.f;
while true,
  tic;
  handles = guidata(hObject);
  if ~ISPLAYING || fnext > handles.nframes,
    break;
  end
  f = fnext;
  handles.f = f;
  handles = update_frame(handles,hObject);
  drawnow;
  guidata(hObject,handles);
  if handles.MaxFPS > 0,
    elapsedtime = toc;
    if elapsedtime < handles.MinSPF,
      dt = handles.MinSPF-elapsedtime;
      pause(dt);
      fnext = f + 1;
    else,
      df = round(elapsedtime*handles.MaxFPS);
      fnext = f + df;
      if fnext > handles.nframes,
        break;
      end
    end
  else
    fnext = f+1;
    drawnow;
  end
  tic;
end

stop(hObject);

function stop(hObject)

global ISPLAYING;

handles = guidata(hObject);
ISPLAYING = false;
guidata(hObject,handles);
set(hObject,'String','Play','BackgroundColor',[0,.5,0]);

% --- Outputs from this function are returned to the command line.
function varargout = playfmf_interestpoints_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.figure1;

function handles = update_frame(handles,hObject)

if nargin == 1,
  hObject = nan;
end
persistent update_frame_f;

if ~isempty(update_frame_f) && update_frame_f == handles.f,
  return;
end
update_frame_f = handles.f;
if handles.f ~= handles.lastfread,
  try
    [handles.im,handles.timestamp] = handles.readframe(handles.f);
    handles.lastfread = handles.f;
    guidata(handles.figure1,handles);
  catch ME
    warning('Failed to read frame %d: %s',handles.f,getReport(ME));
  end  
  set(handles.himage,'CData',handles.im);
end
if hObject ~= handles.edit_Frame,
  set(handles.edit_Frame,'String',num2str(handles.f));
end
if hObject ~= handles.slider_Frame,
  %set(handles.slider_Frame,'Value',(handles.f-1)/(handles.nframes-1));
  set(handles.slider_Frame,'Value',handles.f);
end
update_frame_f = [];

% --- Executes on slider movement.
function slider_Frame_Callback(hObject, eventdata, handles)
% hObject    handle to slider_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

v = get(hObject,'Value');
roundv = round(v);
if roundv == handles.f,
  if v > handles.f,
    handles.f = min(handles.nframes,handles.f+1);
  elseif v < handles.f,
    handles.f = max(1,handles.f-1);
  end
else
  handles.f = roundv;
end

handles = update_frame(handles);
guidata(hObject,handles);

% v = get(hObject,'Value');
% handles.f = round(1 + v * (handles.nframes - 1));
% handles = update_frame(handles,hObject);
% guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function slider_Frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_PlayStop.
function pushbutton_PlayStop_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_PlayStop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global ISPLAYING;

if ~ISPLAYING,
  play(hObject);
else
  stop(hObject);
end

function edit_Frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_Frame as text
%        str2double(get(hObject,'String')) returns contents of edit_Frame as a double
f = str2double(get(hObject,'String'));
if isnan(f),
  set(hObject,'String',num2str(handles.f));
  return;
end
handles.f = max(1,min(f,handles.nframes));
handles = update_frame(handles,hObject);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_Frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_Frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function menu_File_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_File_Open_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File_Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = resetAnnotationInfoIfNec(hObject,handles);
handles = open_fmf(handles);
guidata(hObject,handles);

function handles = open_fmf(handles,tfSelectMovie,tfIsOpen,startframe)

global ISPLAYING;
if ~exist('tfSelectMovie','var')
  tfSelectMovie = true;
end
if ~exist('tfIsOpen','var'),
  tfIsOpen = false;
end
if ~exist('startframe','var'),
  startframe = 1;
end

if ~isempty(ISPLAYING) && ISPLAYING,
  stop(handles.figure1);
  handles = guidata(handles.figure1);
end

if ~tfIsOpen && tfSelectMovie
  handles.filterspec = {  '*.ufmf','MicroFlyMovieFormat (*.ufmf)'; ...
    '*.fmf','FlyMovieFormat (*.fmf)'; ...
    '*.sbfmf','StaticBackgroundFMF (*.sbfmf)'; ...
    '*.avi','AVI (*.avi)'
    '*.mp4','MP4 (*.mp4)'
    '*.mov','MOV (*.mov)'
    '*.mmf','MMF (*.mmf)'
    '*.*','*.*'};
  
  if isfield(handles,'fileext'),
    % default ext is last chosen
    i = find(strcmpi(['*',handles.fileext],handles.filterspec(:,1)),1);
    n = size(handles.filterspec,1);
    if ~isempty(i),
      handles.filterspec = handles.filterspec([i,1:i-1,i+1:n],:);
    end
  end
  
  [filename, pathname] = uigetfile(handles.filterspec, 'Choose FMF video to play',handles.filename);
  if ~ischar(filename),
    return;
  end
  handles.filename = fullfile(pathname,filename);
  [handles.filedir,handles.filenamebase,handles.fileext] = fileparts(handles.filename);
end

if isfield(handles,'fid') && ~isempty(fopen(handles.fid)) && handles.fid > 1,
  fclose(handles.fid);
end
if isfield(handles,'himage') && ishandle(handles.himage),
  delete(handles.himage);
end

if ~tfIsOpen,
  try
    [handles.readframe,handles.nframes,handles.fid,handles.headerinfo] = ...
      get_readframe_fcn(handles.filename);
  catch ME
    s = sprintf('Could not read video %s.',handles.filename);
    uiwait(errordlg(s,'Error opening video'));
    rethrow(ME);
  end
end

% set slider steps
% this seems to be the limit to slider step resolution
handles.stepsize = ceil(handles.nframes/handles.MAXSLIDERRES);
step1 = handles.stepsize/(handles.nframes-1);
sliderstep = [step1,min(1,100*step1)];
set(handles.slider_Frame,'Value',0,'SliderStep',sliderstep,'Min',1,'Max',handles.nframes);

% % set slider steps
% sliderstep = [1/(handles.nframes-1),min(1,100/(handles.nframes-1))];
% set(handles.slider_Frame,'Value',0,'SliderStep',sliderstep);

% show first image
handles.f = min(handles.nframes,max(1,startframe));
[handles.im,handles.timestamp] = handles.readframe(handles.f);
handles.lastfread = handles.f;
if size(handles.im,3) == 1,
  handles.himage = imagesc(handles.im,'Parent',handles.axes_Video,[0,255]);
else
  handles.himage = image(uint8(handles.im),'Parent',handles.axes_Video);
end
colormap(handles.axes_Video,'gray');
axis(handles.axes_Video,'image','off');
handles = update_frame(handles);

set(handles.txMovieName,'String',handles.filename);

if handles.tfipsshow
  flds = fieldnames(handles.ipsshow);
  tfhold = ishold(handles.axes_Video);
  hold(handles.axes_Video,'on');

  for i = 1:numel(flds)
    val = handles.ipsshow.(flds{i});
    if isnumeric(val) && numel(val)==2
      plot(handles.axes_Video,val(1),val(2),'ro','MarkerSize',8,'MarkerFaceColor',[1 0 0]);
    end
  end
  
  if ~tfhold
    hold(handles.axes_Video,'off');
  end
end

ISPLAYING = false;

% --------------------------------------------------------------------
function menu_File_Quit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File_Quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figure1_CloseRequestFcn(handles.figure1, eventdata, handles);

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

global ISPLAYING;
ISPLAYING = false;

if isfield(handles,'fid') && ~isempty(fopen(handles.fid)) && handles.fid > 1,
  fclose(handles.fid);
end

savefns = {'filename','CompressionSettings'};
save(handles.rcfilename,'-struct','handles',savefns{:});
delete(hObject);


% --- Executes when figure1 is resized.
function handles = figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'figpos'),
  return;
end

sliderpos = get(handles.slider_Frame,'Position');
editpos = get(handles.edit_Frame,'Position');
pushbuttonpos = get(handles.pushbutton_PlayStop,'Position');
axespos = get(handles.axes_Video,'Position');
newfigpos = get(handles.figure1,'Position');

woff = newfigpos(3) - handles.figpos(3);
hoff = newfigpos(4) - handles.figpos(4);

sliderpos(3) = sliderpos(3) + woff;
editpos(1) = editpos(1) + woff;
pushbuttonpos(1) = pushbuttonpos(1) + woff;

axespos(3) = axespos(3) + woff;
axespos(4) = axespos(4) + hoff;

set(handles.slider_Frame,'Position',sliderpos);
set(handles.edit_Frame,'Position',editpos);
set(handles.pushbutton_PlayStop,'Position',pushbuttonpos);
set(handles.axes_Video,'Position',axespos);

if handles.tfAnnotate,
  hsmove = {'pbAnnotate','pbDone','pbCancel'};
  for i = 1:numel(hsmove),
    hmove = hsmove{i};
    pos = get(handles.(hmove),'Position');
    pos(1) = pos(1) + woff;
    pos(2) = pos(2) + hoff;
    set(handles.(hmove),'Position',pos);
  end
end
    
handles.figpos = newfigpos;
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_Edit_Preferences_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Edit_Preferences (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompts = {};
defAns = {};
prompts{end+1} = 'Max FPS: ';
defAns{end+1} = num2str(handles.MaxFPS);
while true,
  answer = inputdlg(prompts,'playfmf Preferences',1,defAns,'on');
  if isempty(answer),
    return;
  end
  v = str2double(answer{1});
  iserror = false;
  if isnan(v),
    iserror = true;
  else
    defAns{1} = num2str(v);
    handles.MaxFPS = v;
    handles.MinSPF = 1 / handles.MaxFPS;
  end
  if iserror,
    uiwait(warndlg('Illegal values entered. Max FPS must be a number','Bad Preferences'));
  else
    break;
  end
end
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_Help_About_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Help_About (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = {};
s{end+1} = 'PlayFMF';
s{end+1} = '';
s{end+1} = 'Kristin Branson';
s{end+1} = 'bransonk@janelia.hhmi.org';
s{end+1} = '';
s{end+1} = 'This is a GUI for playing FMF, SBFMF, UFMF, and AVI videos. The maximum frame rate can be set through the "Preferences..." menu. Set to <= 0 for no maximum frame rate. Control the frame shown with the slider or editable text box. The Play/Stop button does what you think it does.';
msgbox(s,'About PlayFMF','help','Replace');

% --------------------------------------------------------------------
function menu_Edit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_Edit_CompressionSettings_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Edit_CompressionSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CompressionSettings = SaveFMFSettings(handles.CompressionSettings);
fns = fieldnames(CompressionSettings);
for i = 1:length(fns),
  fn = fns{i};
  handles.CompressionSettings.(fn) = CompressionSettings.(fn);
end
guidata(hObject,handles);

% --------------------------------------------------------------------
function menu_File_Compress_Callback(hObject, eventdata, handles)
% hObject    handle to menu_File_Compress (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.CompressionSettings.StartFrame > handles.nframes,
  uiwait(errordlg(sprintf('StartFrame set to %d > NFrames = %d. Please fix Compression Settings',...
    handles.CompressionSettings.StartFrame,handles.nframes)));
  return;
end

[pathstr,filestr,ext] = fileparts(handles.filename);
handles.aviname = fullfile(pathstr,[filestr,'.avi']);

[filename,pathname] = uiputfile('*.avi','Save to AVI file',handles.aviname);
if ~ischar(filename),
  return;
end
handles.aviname = fullfile(pathname,filename);

ncolors = size(handles.im,3);
isindexed = ismember(handles.CompressionSettings.Compression,{'MSVC','RLE'}) || ...
  (strcmp(handles.CompressionSettings.Compression,'None') && ...
  ncolors == 1);
if isindexed,
  if ncolors == 1,
    cmap = repmat(linspace(0,1,256)',[1,3]);
  else
    [tmp,cmap] = rgb2ind(handles.im,256,'nodither');
  end
end

params = {'compression',handles.CompressionSettings.Compression,...
  'fps',handles.CompressionSettings.OutputFPS};

if isindexed,
  params(end+1:end+2) = {'colormap',cmap};
end
if ~strcmp(handles.CompressionSettings.Compression,'None'),
  params(end+1:end+2) = {'quality',handles.CompressionSettings.Quality};
end
handles.aviobj = avifile(handles.aviname,params{:});

endframe = min(handles.nframes,handles.CompressionSettings.EndFrame);
nframescompress = endframe - handles.CompressionSettings.StartFrame + 1;

i = 0;
s = sprintf('Compressing frame %d / %d',i,nframescompress);
hwaitbar = waitbar(0,s,'CreateCancelBtn',...
  'setappdata(gcbf,''canceling'',1)');
setappdata(hwaitbar,'canceling',0);
for t = handles.CompressionSettings.StartFrame:endframe,
  if getappdata(hwaitbar,'canceling')
    break
  end
  im = uint8(handles.readframe(t));
  if isindexed,
    if ncolors == 1,
      handles.aviobj = addframe(handles.aviobj,im);
    else
      handles.aviobj = addframe(handles.aviobj,rgb2ind(im,cmap,'nodither'));
    end
  else
    if ncolors == 1,
      handles.aviobj = addframe(handles.aviobj,repmat(im,[1,1,3]));
    else
      handles.aviobj = addframe(handles.aviobj,im);
    end
  end
  i = i + 1;
  if mod(i,50) == 0,
    s = sprintf('Compressing frame %d / %d',i,nframescompress);
    waitbar(i/nframescompress,hwaitbar,s);
  end
end
handles.aviobj = close(handles.aviobj);
if ishandle(hwaitbar),
  delete(hwaitbar);
end
msgbox(sprintf('Successfully compressed %d / %d of frames in the interval [%d,%d].',i,nframescompress,handles.CompressionSettings.StartFrame,endframe),'Compression Complete','modal');

% --------------------------------------------------------------------
function menu_Help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_Help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function handles = resetAnnotationInfoIfNec(hObject,handles)
if handles.tfAnnotate && ~isempty(handles.impts)
  structfun(@delete,handles.impts);
  handles.impts = [];
  guidata(hObject,handles);
  set(handles.pbDone,'Enable','off');
  handles.annotated_frame = [];
end

function pbAnnotate_Callback(hObject, eventdata, handles)

global ISPLAYING;
ISPLAYING = false;

handles = resetAnnotationInfoIfNec(hObject,handles);
  

for i = 1:numel(handles.ipnames),
  title(handles.axes_Video,sprintf('Click on the %s location',handles.ipnames{i}),'fontweight','bold');
  handles.impts.(handles.ipnames{i}) = impoint(handles.axes_Video);
  handles.impts.(handles.ipnames{i}).setColor([1 0 0]);
end

handles.annotated_frame = handles.f;

% title(handles.axes_Video,'Click on the mouth location','fontweight','bold');
% handles.impts.mouth = impoint(handles.axes_Video);
% handles.impts.mouth.setColor([1 0 0]);
% 
% title(handles.axes_Video,'Click on the perch location','fontweight','bold');
% handles.impts.perch = impoint(handles.axes_Video);
% handles.impts.perch.setColor([1 0 0]);

if handles.doaskleftright,
  res = questdlg('Is the mouse facing right?','Choose direction','Yes','No','Yes');
  if strcmpi(res,'yes')
    handles.lr = 'right';
  elseif strcmpi(res,'no'),
    handles.lr = 'left';
  else
    structfun(@delete,handles.impts);
    handles.impts = [];
    error('playfmf:choosedir','A direction must be chosen.');
  end
  tstr = sprintf('Face %s: points may be adjusted.',handles.lr);
else
  handles.lr = 'unknown';
  tstr = 'Drag points to adjust their locations';

end

title(handles.axes_Video,tstr,'fontweight','bold');

guidata(hObject,handles);
set(handles.pbDone,'Enable','on');

function pbDone_Callback(hObject, eventdata, handles)
global ISPLAYING;
ISPLAYING = false;

% successful annotation
s = struct();
for i = 1:numel(handles.ipnames),
  ipname = handles.ipnames{i};
  try
    s.(ipname) = handles.impts.(ipname).getPosition();
  catch ME,
    warning(getReport(ME));
  end
end
s.lr = handles.lr;
s.annotated_frame = handles.annotated_frame;
s.figpos = handles.figpos;
handles.annDataObj.data = s;
figure1_CloseRequestFcn(handles.figure1, [], handles);


% --- Executes on button press in pbCancel.
function pbCancel_Callback(hObject, eventdata, handles)
figure1_CloseRequestFcn(handles.figure1, eventdata, handles);
