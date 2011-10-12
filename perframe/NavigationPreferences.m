function varargout = NavigationPreferences(varargin)
% NAVIGATIONPREFERENCES MATLAB code for NavigationPreferences.fig
%      NAVIGATIONPREFERENCES, by itself, creates a new NAVIGATIONPREFERENCES or raises the existing
%      singleton*.
%
%      H = NAVIGATIONPREFERENCES returns the handle to a new NAVIGATIONPREFERENCES or the handle to
%      the existing singleton*.
%
%      NAVIGATIONPREFERENCES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NAVIGATIONPREFERENCES.M with the given input arguments.
%
%      NAVIGATIONPREFERENCES('Property','Value',...) creates a new NAVIGATIONPREFERENCES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NavigationPreferences_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NavigationPreferences_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NavigationPreferences

% Last Modified by GUIDE v2.5 19-Sep-2011 13:55:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NavigationPreferences_OpeningFcn, ...
                   'gui_OutputFcn',  @NavigationPreferences_OutputFcn, ...
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


% --- Executes just before NavigationPreferences is made visible.
function NavigationPreferences_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NavigationPreferences (see VARARGIN)

% first input should be the parent figure from JLabel
handles.figure_JLabel = varargin{1};

% read current values
parent_handles = guidata(handles.figure_JLabel);
handles.nframes_jump_go = parent_handles.nframes_jump_go;
handles.seek_behaviors_go = parent_handles.seek_behaviors_go;
handles.behaviors = [{'Unknown'},parent_handles.data.labelnames];

% set these current values in the GUI
set(handles.edit_nframes_jump,'String',num2str(handles.nframes_jump_go));
set(handles.listbox_seek_behavior,'String',handles.behaviors);
set(handles.listbox_seek_behavior,'Value',handles.seek_behaviors_go+1);

% Choose default command line output for NavigationPreferences
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NavigationPreferences wait for user response (see UIRESUME)
% uiwait(handles.figure_NavigationPreferences);


% --- Outputs from this function are returned to the command line.
function varargout = NavigationPreferences_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_nframes_jump_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to edit_nframes_jump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nframes_jump as text
%        str2double(get(hObject,'String')) returns contents of edit_nframes_jump as a double
v = str2double(get(hObject,'String'));
if isnan(v) || round(v) ~= v || v <= 0,
  warndlg('N. frames jump must be a positive integer','Bad value');
  set(hObject,'String',num2str(handles.nframes_jump_go));
else
  handles.nframes_jump_go = v;
  guidata(hObject,handles);
end

% --- Executes during object creation, after setting all properties.
function edit_nframes_jump_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to edit_nframes_jump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_seek_behavior.
function listbox_seek_behavior_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_seek_behavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_seek_behavior contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_seek_behavior
v = get(hObject,'Value');
if isempty(v),
  warndlg('At least one must behavior must be selected','Bad value');
  set(hObject,'Value',handles.seek_behaviors_go+1);
else
  handles.seek_behaviors_go = v-1;
  guidata(hObject,handles);
end


% --- Executes during object creation, after setting all properties.
function listbox_seek_behavior_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_seek_behavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent_handles = guidata(handles.figure_JLabel);
parent_handles.nframes_jump_go = handles.nframes_jump_go;
parent_handles.seek_behaviors_go = handles.seek_behaviors_go;
guidata(handles.figure_JLabel,parent_handles);
JLabel('SetJumpGoMenuLabels',parent_handles);
delete(handles.figure_NavigationPreferences);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

delete(handles.figure_NavigationPreferences);

% --- Executes on button press in pushbutton_apply.
function pushbutton_apply_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

parent_handles = guidata(handles.figure_JLabel);
parent_handles.nframes_jump_go = handles.nframes_jump_go;
parent_handles.seek_behaviors_go = handles.seek_behaviors_go;
guidata(handles.figure_JLabel,parent_handles);
JLabel('SetJumpGoMenuLabels',parent_handles);
