function varargout = ConfidenceThresholds(varargin)
% CONFIDENCETHRESHOLDS MATLAB code for ConfidenceThresholds.fig
%      CONFIDENCETHRESHOLDS, by itself, creates a new CONFIDENCETHRESHOLDS or raises the existing
%      singleton*.
%
%      H = CONFIDENCETHRESHOLDS returns the handle to a new CONFIDENCETHRESHOLDS or the handle to
%      the existing singleton*.
%
%      CONFIDENCETHRESHOLDS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONFIDENCETHRESHOLDS.M with the given input arguments.
%
%      CONFIDENCETHRESHOLDS('Property','Value',...) creates a new CONFIDENCETHRESHOLDS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ConfidenceThresholds_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ConfidenceThresholds_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ConfidenceThresholds

% Last Modified by GUIDE v2.5 10-Jan-2012 17:47:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ConfidenceThresholds_OpeningFcn, ...
                   'gui_OutputFcn',  @ConfidenceThresholds_OutputFcn, ...
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


% --- Executes just before ConfidenceThresholds is made visible.
function ConfidenceThresholds_OpeningFcn(hObject, eventdata, handles)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ConfidenceThresholds (see VARARGIN)

% Choose default command line output for ConfidenceThresholds
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ConfidenceThresholds wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ConfidenceThresholds_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function SetJLabelData(hObject,obj,JLabelHandle)
handles = guidata(hObject);
handles.JLDobj = obj;
handles.JLabelHandle = JLabelHandle;
guidata(hObject,handles);

function SetSliderColor(hObject,sliderNum,color)
handles = guidata(hObject);
curSlider = sprintf('slider%d',sliderNum);
set(handles.(curSlider),'BackgroundColor',color);
guidata(hObject,handles);

function SetConfidenceThreshold(hObject,sliderNum,value)
handles = guidata(hObject);
curSlider = sprintf('slider%d',sliderNum);
set(handles.(curSlider),'Value',value);
curEdit = sprintf('edit%d',sliderNum);
set(handles.(curEdit),'String',num2str(value));
guidata(hObject,handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
set(handles.edit1,'String',num2str(val));
handles.JLDobj.SetConfidenceThreshold(val,1);
JLabel('predict',handles.JLabelHandle);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
val = get(hObject,'Value');
set(handles.edit2,'String',num2str(val));
handles.JLDobj.SetConfidenceThreshold(val,2);
JLabel('predict',handles.JLabelHandle);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val1 = get(handles.slider1,'value');
val2 = get(handles.slider2,'value');
handles.JLDobj.SetConfidenceThreshold(val1,1);
handles.JLDobj.SetConfidenceThreshold(val2,2);
JLabel('predict',handles.JLabelHandle);
close(handles.figure1);



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
val = str2double(get(hObject,'String'));
if isnan(val)
  errordlg('Input value between 0 and 1');
  return;
end
set(handles.slider1,'Value',val);
handles.JLDobj.SetConfidenceThreshold(val,1);
JLabel('predict',handles.JLabelHandle);



% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
val = str2double(get(hObject,'String'));
if isnan(val)
  errordlg('Input value between 0 and 1');
  return;
end
set(handles.slider2,'Value',val);
handles.JLDobj.SetConfidenceThreshold(val,2);
JLabel('predict',handles.JLabelHandle);


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
