function varargout = ClassifierChange(varargin)
% CLASSIFIERCHANGE MATLAB code for ClassifierChange.fig
%      CLASSIFIERCHANGE, by itself, creates a new CLASSIFIERCHANGE or raises the existing
%      singleton*.
%
%      H = CLASSIFIERCHANGE returns the handle to a new CLASSIFIERCHANGE or the handle to
%      the existing singleton*.
%
%      CLASSIFIERCHANGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLASSIFIERCHANGE.M with the given input arguments.
%
%      CLASSIFIERCHANGE('Property','Value',...) creates a new CLASSIFIERCHANGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ClassifierChange_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ClassifierChange_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ClassifierChange

% Last Modified by GUIDE v2.5 17-Jan-2012 10:53:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ClassifierChange_OpeningFcn, ...
                   'gui_OutputFcn',  @ClassifierChange_OutputFcn, ...
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


% --- Executes just before ClassifierChange is made visible.
function ClassifierChange_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ClassifierChange (see VARARGIN)

% Choose default command line output for ClassifierChange
handles.output = hObject;

handles.JLDObj = varargin{1};
handles.JLabelhandles = varargin{2};
% set(handles.axes1,'axes',off);
handles.curExp = handles.JLabelhandles.expi;
handles.curFly = handles.JLabelhandles.flies;
set(handles.listboxExp,'String',handles.JLDObj.expnames,'Value',handles.curExp,'max',0,'min',0);
guidata(hObject, handles);
UpdateListBoxExp(handles);
% Update handles structure

% UIWAIT makes ClassifierChange wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ClassifierChange_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listboxExp.
function listboxExp_Callback(hObject, eventdata, handles)
% hObject    handle to listboxExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxExp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxExp
selExp = get(handles.listboxExp,'Value');
if selExp == handles.curExp; return; end
handles.curExp = selExp;
handles.curFly = 1;
% Update handles structure
guidata(hObject, handles);
UpdateListBoxExp(handles);


function UpdateListBoxExp(handles)
selExp = handles.curExp;
s = cell(1,handles.JLDObj.nflies_per_exp(selExp));
for flyNum = 1:handles.JLDObj.nflies_per_exp(selExp)
  flyStats = handles.JLDObj.GetFlyStats(selExp,flyNum);
    s{flyNum} = sprintf('Fly %3d, Trajectory length %5d, First frame %5d, N bouts labeled %2d',...
      flyNum,flyStats.trajLength,flyStats.firstframe,flyStats.nbouts);
    if flyStats.hassex,
      if ~isempty(flyStats.sexfrac),
        s{flyNum} = [s{flyNum},sprintf(', Sex: %3d%%M, %3d%%F',...
          round(flyStats.sexfrac.M*100),round(flyStats.sexfrac.F*100))];
      else
        s{flyNum} = [s{flyNum},sprintf(', Sex: %s',flyStats.sex)];
      end
    end
    if ~isempty(flyStats.nscoreframes)
      s{flyNum} = [s{flyNum},sprintf(', Frames Predicted as %s:%d, Total Frames Predicted:%d',...
        handles.JLDObj.labelnames{1},flyStats.nscorepos,flyStats.nscoreframes)];
    end
end
set(handles.listboxFly,'String',s,'Value',handles.curFly);



% --- Executes during object creation, after setting all properties.
function listboxExp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxExp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listboxFly.
function listboxFly_Callback(hObject, eventdata, handles)
% hObject    handle to listboxFly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listboxFly contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listboxFly
handles.curFly = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function listboxFly_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listboxFly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushSwitchfly.
function pushSwitchfly_Callback(hObject, eventdata, handles)
% hObject    handle to pushSwitchfly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
JLabel('SetCurrentMovie',handles.JLabelhandles,handles.curExp);
JLabel('SetCurrentFlies',handles.JLabelhandles,handles.curFly);


% --- Executes on button press in pushClose.
function pushClose_Callback(hObject, eventdata, handles)
% hObject    handle to pushClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);