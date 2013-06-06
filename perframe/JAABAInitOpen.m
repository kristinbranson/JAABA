function varargout = JAABAInitOpen(varargin)
% JAABAINITOPEN MATLAB code for JAABAInitOpen.fig
%      JAABAINITOPEN, by itself, creates a new JAABAINITOPEN or raises the existing
%      singleton*.
%
%      H = JAABAINITOPEN returns the handle to a new JAABAINITOPEN or the handle to
%      the existing singleton*.
%
%      JAABAINITOPEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JAABAINITOPEN.M with the given input arguments.
%
%      JAABAINITOPEN('Property','Value',...) creates a new JAABAINITOPEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JAABAInitOpen_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JAABAInitOpen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JAABAInitOpen

% Last Modified by GUIDE v2.5 04-Jun-2013 16:23:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JAABAInitOpen_OpeningFcn, ...
                   'gui_OutputFcn',  @JAABAInitOpen_OutputFcn, ...
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


% --- Executes just before JAABAInitOpen is made visible.
function JAABAInitOpen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JAABAInitOpen (see VARARGIN)

% Choose default command line output for JAABAInitOpen
handles.output = hObject;
handles.sel = struct('val',[],'edit',false);

figureJLabel = varargin{1};
% Change a few things so they still work well on Mac
adjustColorsIfMac(hObject);

centerOnParentFigure(hObject,figureJLabel);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JAABAInitOpen wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = JAABAInitOpen_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.sel;
delete(handles.figure1);

% --- Executes on button press in pushbutton_new.
function pushbutton_new_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_new (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sel.val = 'New';
guidata(hObject,handles);
figure1_CloseRequestFcn(handles.figure1,[],handles);

% --- Executes on button press in pushbutton_open.
function pushbutton_open_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sel.val = 'Open';
guidata(hObject,handles);
figure1_CloseRequestFcn(handles.figure1,[],handles);


% --- Executes on button press in pushbutton_opengt.
function pushbutton_opengt_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_opengt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.sel.val = 'OpenGT';
guidata(hObject,handles);
figure1_CloseRequestFcn(handles.figure1,[],handles);


% --- Executes on button press in checkbox_edit.
function checkbox_edit_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_edit

handles.sel.edit = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
% delete(hObject);
uiresume(hObject);
