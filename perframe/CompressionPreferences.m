function varargout = CompressionPreferences(varargin)
% COMPRESSIONPREFERENCES MATLAB code for CompressionPreferences.fig
%      COMPRESSIONPREFERENCES, by itself, creates a new COMPRESSIONPREFERENCES or raises the existing
%      singleton*.
%
%      H = COMPRESSIONPREFERENCES returns the handle to a new COMPRESSIONPREFERENCES or the handle to
%      the existing singleton*.
%
%      COMPRESSIONPREFERENCES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPRESSIONPREFERENCES.M with the given input arguments.
%
%      COMPRESSIONPREFERENCES('Property','Value',...) creates a new COMPRESSIONPREFERENCES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CompressionPreferences_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CompressionPreferences_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CompressionPreferences

% Last Modified by GUIDE v2.5 24-Nov-2011 13:22:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CompressionPreferences_OpeningFcn, ...
                   'gui_OutputFcn',  @CompressionPreferences_OutputFcn, ...
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


% --- Executes just before CompressionPreferences is made visible.
function CompressionPreferences_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CompressionPreferences (see VARARGIN)

% Choose default command line output for CompressionPreferences
handles.parent_figure = varargin{1};

handles = Initialize(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CompressionPreferences wait for user response (see UIRESUME)
uiwait(handles.figure_CompressionPreferences);

function handles = Initialize(handles)

% set last values
% useVideoWriter
% outavi_compression -> compression_format: profile (VideoWriter) or
% compression (avifile) 
% outavi_fps -> fps
% outavi_quality -> quality: Quality (VideoWriter + Motion JPEG AVI) or
% Compression Ratio and LosslessCompression (VideoWriter + !Motion JPEG AVI
% or quality (avifile)
parent_handles = guidata(handles.parent_figure);

if isfield(parent_handles,'useVideoWriter'),
  if exist('VideoWriter','file'),
    handles.useVideoWriter = parent_handles.useVideoWriter;
  else
    handles.useVideoWriter = true;
  end
else
  handles.useVideoWriter = exist('VideoWriter','file');
end

if isfield(parent_handles,'outavi_compression'),
  handles.compression_format = parent_handles.outavi_compression;
else
  handles.compression_format = 'None';
end

if isfield(parent_handles,'outavi_fps'),
  handles.fps = parent_handles.outavi_fps;
else
  handles.fps = 15;
end

if isfield(parent_handles,'outavi_quality'),
  handles.quality = parent_handles.outavi_quality;
else
  handles.quality = 100;
end

handles = SetCompressionFormats(handles);

handles = UpdateGUIPrompts(handles);

set(handles.radiobutton_useVideoWriter,'Value',handles.useVideoWriter);
set(handles.radiobutton_useAviFile,'Value',~handles.useVideoWriter);
%i = find(strcmpi(handles.compression_format,handles.compression_formats),1);
%set(handles.popupmenu_compression_profile,'Value',i);
set(handles.edit_fps,'String',num2str(handles.fps));
set(handles.edit_quality,'String',num2str(handles.quality));

function handles = UpdateGUIPrompts(handles)

% useVideoWriter
% outavi_compression -> compression_format: profile (VideoWriter) or
% compression (avifile) 
% outavi_fps -> fps
% outavi_quality -> quality: Quality (VideoWriter + Motion JPEG AVI) or
% Compression Ratio and LosslessCompression (VideoWriter + ~Motion JPEG AVI
% or quality (avifile)

if ~exist('VideoWriter','file'),
  if handles.useVideoWriter,
    warning('This should never happen: useVideoWriter == true and VideoWriter does not exist');
    handles.useVideoWriter = false;
  end
  set(handles.uipanel_advanced,'Visible','off');
else
  set(handles.uipanel_advanced,'Visible','on');
end

if handles.useVideoWriter && strcmpi(handles.compression_format,'Motion JPEG AVI'),
  set(handles.text_quality,'String','Compression Quality');
else
  set(handles.text_quality,'String','Compression Ratio');
end

if (handles.useVideoWriter && ...
    ismember(lower(handles.compression_format),{'archival','none'})) || ...
    (~handles.useVideoWriter && strcmpi(handles.compression_format,'None')),
  set(handles.text_quality,'Enable','off');
  set(handles.edit_quality,'Enable','off');
else
  set(handles.text_quality,'Enable','on');
  set(handles.edit_quality,'Enable','on');
end

function handles = SetCompressionFormats(handles)

if handles.useVideoWriter,
  profiles = VideoWriter.getProfiles();
  handles.compression_formats = {profiles.Name};
  isvalid = profiles.isvalid;
  if ~all(isvalid),
    handles.compresion_formats(~isvalid) = [];
  end
  i = find(strcmpi(handles.compression_formats,'Uncompressed AVI'),1);
  if ~isempty(i),
    handles.compression_formats{i} = 'None';
  end
else
  if ispc,
    compressors = mmcompinfo;
    handles.compression_formats = [{'None'},{compressors.videoCompressors.Name}];
  else
    handles.compression_formats = {'None'};
  end
end

i = find(strcmpi(handles.compression_formats,handles.compression_format),1);
if isempty(i),
  i = 1;
  handles.compression_format = handles.compression_formats{1};
end

set(handles.popupmenu_compression_profile,'String',handles.compression_formats,'Value',i);

  
% --- Outputs from this function are returned to the command line.
function varargout = CompressionPreferences_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout = {};
delete(handles.figure_CompressionPreferences);

% --- Executes on button press in radiobutton_useVideoWriter.
function radiobutton_useVideoWriter_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to radiobutton_useVideoWriter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_useVideoWriter
handles.useVideoWriter = get(handles.radiobutton_useVideoWriter,'Value');
set(handles.radiobutton_useAviFile,'Value',~handles.useVideoWriter);
% translatable compression format
handles.compression_format = 'none';
handles.quality = 100;
handles = SetCompressionFormats(handles);
handles = UpdateGUIPrompts(handles);
guidata(hObject,handles);

% --- Executes on button press in radiobutton_useAviFile.
function radiobutton_useAviFile_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_useAviFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_useAviFile
handles.useVideoWriter = get(handles.radiobutton_useAviFile,'Value') == 0;
set(handles.radiobutton_useVideoWriter,'Value',handles.useVideoWriter);
% translatable compression format
handles.compression_format = 'none';
handles.quality = 100;
handles = SetCompressionFormats(handles);
handles = UpdateGUIPrompts(handles);
guidata(hObject,handles);

% --- Executes on selection change in popupmenu_compression_profile.
function popupmenu_compression_profile_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_compression_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_compression_profile contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_compression_profile
old_compression_format = handles.compression_format;
handles.compression_format = handles.compression_formats{get(handles.popupmenu_compression_profile,'Value')};
% translate quality
if handles.useVideoWriter && ...
    ismember(lower(handles.compression_format),{'archival','none'}),
  handles.quality = 1;
  set(handles.edit_quality,'String',num2str(handles.quality));  
elseif ~handles.useVideoWriter && strcmpi(handles.compression_format,'None'),
  handles.quality = 100;
  set(handles.edit_quality,'String',num2str(handles.quality));
elseif strcmpi(old_compression_format,'Motion JPEG AVI') && ...
    ~strcmpi(handles.compression_format,'Motion JPEG AVI'),
  % >99 -> 1, <=99 -> 10
  if handles.quality >= 99,
    handles.quality = 1;
  else
    handles.quality = 10;
  end
  set(handles.edit_quality,'String',num2str(handles.quality));
elseif ~strcmpi(old_compression_format,'Motion JPEG AVI') && ...
    strcmpi(handles.compression_format,'Motion JPEG AVI'),
  % <= 1.5 -> 100, >1.5 -> 75
  if handles.quality <= 1.5,
    handles.quality = 100;
  else
    handles.quality = 75;
  end
  set(handles.edit_quality,'String',num2str(handles.quality));
end

handles = UpdateGUIPrompts(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_compression_profile_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to popupmenu_compression_profile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_quality_Callback(hObject, eventdata, handles)
% hObject    handle to edit_quality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_quality as text
%        str2double(get(hObject,'String')) returns contents of edit_quality as a double

v = str2double(get(hObject,'String'));
if isnan(v),
  uiwait(warndlg(sprintf('%s must be a number',get(handles.text_quality,'String'))));
  set(hObject,'String',num2str(handles.quality));
  return;
end
if handles.useVideoWriter && ~strcmpi(handles.compression_format,'Motion JPEG AVI'),
  if v < 1,
    uiwait(warndlg('Compression ratio should be >= 1'));
    set(hObject,'String',num2str(handles.quality));
    return;
  end
else
  if v <= 0 || v > 100,
    uiwait(warndlg('Compression quality should be > 0 and <= 100'));
    set(hObject,'String',num2str(handles.quality));
    return;
  end
end
handles.quality = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_quality_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_quality (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% store current configuration
parent_handles = guidata(handles.parent_figure);
parent_handles.outavi_fps = handles.fps;
parent_handles.outavi_compression = handles.compression_format;
parent_handles.useVideoWriter = handles.useVideoWriter;
parent_handles.outavi_quality = handles.quality;
guidata(handles.parent_figure,parent_handles);
delete(handles.figure_CompressionPreferences);

% --- Executes during object creation, after setting all properties.
function pushbutton_done_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure_CompressionPreferences);



function edit_fps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fps as text
%        str2double(get(hObject,'String')) returns contents of edit_fps as a double

v = str2double(get(hObject,'String'));
if isnan(v) || v <= 0 || round(v) ~= v,
  uiwait(warndlg('Frames per second must be a whole number greater than zero'));
  set(hObject,'String',num2str(handles.fps));
  return;
end
handles.fps = v;
guidata(hObject,handles);