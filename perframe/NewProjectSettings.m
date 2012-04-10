function varargout = NewProjectSettings(varargin)
% NEWPROJECTSETTINGS MATLAB code for NewProjectSettings.fig
%      NEWPROJECTSETTINGS, by itself, creates a new NEWPROJECTSETTINGS or raises the existing
%      singleton*.
%
%      H = NEWPROJECTSETTINGS returns the handle to a new NEWPROJECTSETTINGS or the handle to
%      the existing singleton*.
%
%      NEWPROJECTSETTINGS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEWPROJECTSETTINGS.M with the given input arguments.
%
%      NEWPROJECTSETTINGS('Property','Value',...) creates a new NEWPROJECTSETTINGS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before NewProjectSettings_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to NewProjectSettings_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help NewProjectSettings

% Last Modified by GUIDE v2.5 05-Apr-2012 13:16:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @NewProjectSettings_OpeningFcn, ...
                   'gui_OutputFcn',  @NewProjectSettings_OutputFcn, ...
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


% --- Executes just before NewProjectSettings is made visible.
function NewProjectSettings_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to NewProjectSettings (see VARARGIN)

% Choose default command line output for NewProjectSettings
handles.output = hObject;
handles.prevProjects = varargin{1};

handles.projname = '';
handles.configfile = '';
handles.defaultConfigfile = '';
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NewProjectSettings wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = NewProjectSettings_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.projname;
varargout{3} = handles.configfile;
varargout{4} = handles.defaultConfigfile;
delete(handles.figure1);

function editproj_Callback(hObject, eventdata, handles)
% hObject    handle to editproj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editproj as text
%        str2double(get(hObject,'String')) returns contents of editproj as a double
t = get(hObject,'String');
if iscell(t)
  projname = t{1};
else
  projname = t;
end
if any(strcmp(projname,handles.prevProjects))
  uiwait(warndlg('A project already exists with that name'));
  return
end

handles.projname = projname;
handles.configfile = fullfile('params',[projname '_params.xml']);
set(handles.editconfig,'String',handles.configfile);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function editproj_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editproj (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editconfig_Callback(hObject, eventdata, handles)
% hObject    handle to editconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editconfig as text
%        str2double(get(hObject,'String')) returns contents of editconfig as a double
configfile = get(hObject,'String');
handles.configfile = configfile;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editconfig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editdefaultconfig_Callback(hObject, eventdata, handles)
% hObject    handle to editdefaultconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editdefaultconfig as text
%        str2double(get(hObject,'String')) returns contents of editdefaultconfig as a double
defaultConfigfile = get(hObject,'String');
handles.defaultConfigfile = defaultConfigfile;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function editdefaultconfig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editdefaultconfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonConfig.
function pushbuttonConfig_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonConfig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
configfile = handles.configfile;
[fname,pname] = uigetfile(configfile,'Select configuration file..');
if ~fname, return; end
handles.configfile = fullfile(pname,fname);
set(handles.editconfig,'String',handles.configfile);
guidata(hObject,handles);

% --- Executes on button press in pushbuttonDefault.
function pushbuttonDefault_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname] = uigetfile('params/*.xml','Copy configuration settings from..');
if ~fname, return; end
handles.defaultConfigfile = fullfile(pname,fname);
set(handles.editdefaultconfig,'String',handles.defaultConfigfile);
guidata(hObject,handles);


% --- Executes on button press in pushbuttonDone.
function pushbuttonDone_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonDone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);

% --- Executes on button press in pushbuttonCancel.
function pushbuttonCancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.projname = '';
handles.configfile = '';
handles.defaultConfigfile = '';
uiresume(handles.figure1);
