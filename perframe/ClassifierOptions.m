function varargout = ClassifierOptions(varargin)
% CLASSIFIEROPTIONS MATLAB code for ClassifierOptions.fig
%      CLASSIFIEROPTIONS, by itself, creates a new CLASSIFIEROPTIONS or raises the existing
%      singleton*.
%
%      H = CLASSIFIEROPTIONS returns the handle to a new CLASSIFIEROPTIONS or the handle to
%      the existing singleton*.
%
%      CLASSIFIEROPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLASSIFIEROPTIONS.M with the given input arguments.
%
%      CLASSIFIEROPTIONS('Property','Value',...) creates a new CLASSIFIEROPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ClassifierOptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ClassifierOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ClassifierOptions

% Last Modified by GUIDE v2.5 12-Jan-2012 14:25:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ClassifierOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @ClassifierOptions_OutputFcn, ...
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


% --- Executes just before ClassifierOptions is made visible.
function ClassifierOptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ClassifierOptions (see VARARGIN)



handles.data = varargin{1};
handles.classifier_params = handles.data.classifier_params;
set(handles.edit_iter,'String',num2str(handles.data.classifier_params.iter));
set(handles.iter_updates,'String',num2str(handles.data.classifier_params.iter_updates));
set(handles.edit_sample,'String',num2str(handles.data.classifier_params.numSample));
set(handles.edit_bins,'String',num2str(handles.data.classifier_params.numBins));
set(handles.folds,'String',num2str(handles.data.classifier_params.CVfolds));
set(handles.popup_baseclassifier,'String',handles.data.classifier_params.baseClassifierTypes);

% Choose default command line output for ClassifierOptions
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ClassifierOptions wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ClassifierOptions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
delete(handles.figure1);



function edit_iter_Callback(hObject, eventdata, handles)
% hObject    handle to edit_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_iter as text
%        str2double(get(hObject,'String')) returns contents of edit_iter as a double
val = str2double(get(hObject,'String'));
if isempty(val) || double(uint64(val))~=val 
  warndlg('Enter positive value for number of iterations');
  set(hObject,'String',sprintf('%d',handles.classifier_params.iter));
  return;
end;
handles.classifier_params.iter = val;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_iter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_iter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_sample_Callback(hObject, eventdata, handles)
% hObject    handle to edit_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_sample as text
%        str2double(get(hObject,'String')) returns contents of edit_sample as a double
val = str2double(get(hObject,'String'));
if isempty(val) || double(uint64(val))~=val 
  warndlg('Enter positive value for number of sample points');
  set(hObject,'String',sprintf('%d',handles.classifier_params.numSample));
  return;
end;
handles.classifier_params.numSample = val;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_sample_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_sample (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_bins_Callback(hObject, eventdata, handles)
% hObject    handle to edit_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_bins as text
%        str2double(get(hObject,'String')) returns contents of edit_bins as a double
val = str2double(get(hObject,'String'));
if isempty(val) || double(uint64(val))~=val 
  warndlg('Enter positive integer value for number of iterations');
  set(hObject,'String',sprintf('%d',handles.classifier_params.numBins));
  return;
end;

if val>255 || val<10
  warndlg('Maximum number of bins is 255 and minimum number is 10');
  set(hObject,'String',sprintf('%d',handles.classifier_params.numBins));
  return;
end;
handles.classifier_params.numBins = val;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_bins_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_bins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_cancel.
function button_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to button_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);


% --- Executes on button press in pushbutton_apply.
function pushbutton_apply_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_apply (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.data.classifier_params.numBins ~= handles.classifier_params.numBins,
  handles.data.windowdata.binVals = [];
  handles.data.windowdata.bins = [];  
end
handles.data.classifier_params = handles.classifier_params;
uiresume(handles.figure1);


% --- Executes on selection change in popup_baseclassifier.
function popup_baseclassifier_Callback(hObject, eventdata, handles)
% hObject    handle to popup_baseclassifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_baseclassifier contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_baseclassifier
handles.classifier_params.baseClassifierSelected = get(hObject,'Value');

% --- Executes during object creation, after setting all properties.
function popup_baseclassifier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_baseclassifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function folds_Callback(hObject, eventdata, handles)
% hObject    handle to folds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of folds as text
%        str2double(get(hObject,'String')) returns contents of folds as a double
val = str2double(get(hObject,'String'));
if isempty(val) || double(uint64(val))~=val 
  warndlg('Enter positive value for number of iterations');
  set(hObject,'String',sprintf('%d',handles.classifier_params.iter));
  return;
end;
handles.classifier_params.CVfolds = val;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function folds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to folds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function iter_updates_Callback(hObject, eventdata, handles)
% hObject    handle to iter_updates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of iter_updates as text
%        str2double(get(hObject,'String')) returns contents of iter_updates as a double
val = str2double(get(hObject,'String'));
if isempty(val) || double(uint64(val))~=val 
  warndlg('Enter positive value for number of iterations');
  set(hObject,'String',sprintf('%d',handles.classifier_params.iter_updates));
  return;
end;
handles.classifier_params.iter_updates = val;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function iter_updates_CreateFcn(hObject, eventdata, handles)
% hObject    handle to iter_updates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
