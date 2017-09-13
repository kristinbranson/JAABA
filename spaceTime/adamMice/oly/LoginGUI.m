function varargout = LoginGUI(varargin)
% LOGINGUI MATLAB code for LoginGUI.fig
%      LOGINGUI, by itself, creates a new LOGINGUI or raises the existing
%      singleton*.
%
%      H = LOGINGUI returns the handle to a new LOGINGUI or the handle to
%      the existing singleton*.
%
%      LOGINGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LOGINGUI.M with the given input arguments.
%
%      LOGINGUI('Property','Value',...) creates a new LOGINGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LoginGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LoginGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LoginGUI

% Last Modified by GUIDE v2.5 30-May-2011 10:59:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LoginGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @LoginGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [], ...
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

% loginGUI(fcnHandle)
% fcnHandle has signature fcnHandle(username) and is executed when the
% Login button is pressed. If the window is killed, fcnHandle is not
% called.
function LoginGUI_OpeningFcn(hObject, eventdata, handles, varargin)
assert(numel(varargin)==1,'Expected exactly one input argument.');
handles.output = hObject;
handles.loginFcn = varargin{1};
guidata(hObject, handles);
centerfig(hObject);
set(handles.figure1,'WindowStyle','modal')

function varargout = LoginGUI_OutputFcn(hObject, eventdata, handles) %#ok<*INUSL>

% putting these two lines here rather than in OpeningFcn gets the GUI to
% open with the focus in etUsername
uicontrol(handles.etUsername); 
uiwait(handles.figure1);

varargout{1} = [];

function etUsername_Callback(hObject, eventdata, handles) %#ok<*INUSD,*DEFNU>
% If user enters a username and hits return, move focus to pbLogin
c = get(handles.figure1,'CurrentCharacter');
if strcmp(c,sprintf('\r'))
    uicontrol(handles.pbLogin);
end

function edit2_Callback(hObject, eventdata, handles)

function pbLogin_Callback(hObject, eventdata, handles)
username = get(handles.etUsername,'String');
username = strtrim(username);
if isempty(username)
    errordlg('Please enter a username.','Missing Username','modal');
    return;
end
delete(handles.figure1);
handles.loginFcn(username);

function figure1_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'escape'
        delete(handles.figure1);        
end

function pbLogin_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'return'
        pbLogin_Callback(hObject,[],handles);
    case 'escape'
        % this is not working
        set(handles.figure1,'CurrentObject',handles.figure1);
end
