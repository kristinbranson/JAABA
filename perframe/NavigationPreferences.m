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

% Last Modified by GUIDE v2.5 24-Jul-2012 14:56:46

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

% if ismac, % On mac change the foreground color to black.
%   allpopups = findall(hObject,'Style','popup');
%   set(allpopups,'ForegroundColor',[0 0 0]);
% end

handles.figure_JLabel = varargin{1};
handles.NJObj = varargin{2};

% read current values
parent_handles = guidata(handles.figure_JLabel);
handles.nframes_jump_go = parent_handles.guidata.nframes_jump_go;
handles.seek_behaviors_go = handles.NJObj.GetSeekBehaviorsGo();
labelNamesCapitalized= ...
  cellfun(@upperFirstLowerRest, ...
          parent_handles.data.labelnames, ...
          'UniformOutput',false);
handles.behaviors = [{'Unknown'},labelNamesCapitalized];
% set these current values in the GUI
set(handles.edit_nframes_jump,'String',num2str(handles.nframes_jump_go));
set(handles.listbox_seek_behavior, ...
    'String',labelNamesCapitalized, ...
    'Value',handles.seek_behaviors_go, ...
    'ListboxTop',1);
set(handles.jumpToPopUp,'String',handles.NJObj.GetAllTypes());
jumptondx = find(strcmp(handles.NJObj.GetCurrentType(),handles.NJObj.GetAllTypes()));
if numel(jumptondx) ==1
  set(handles.jumpToPopUp,'Value',jumptondx);
else
  set(handles.jumpToPopUp,'Value',1);
end

set(handles.thresholdPopup1,'String',[{'Select'}; handles.NJObj.GetPerframefns]);
updateThresholdButtons(handles);
updateConfidenceButtons(handles);
% Choose default command line output for NavigationPreferences
handles.output = hObject;

% Adjust all the widget colors to look OK on Mac
adjustColorsIfMac(hObject);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes NavigationPreferences wait for user response (see UIRESUME)
% uiwait(handles.figure_NavigationPreferences);
return


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
  set(hObject,'Value',handles.seek_behaviors_go);
else
  handles.seek_behaviors_go = v;
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

pushbutton_apply_Callback(hObject,eventdata,handles);
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
parent_handles.guidata.nframes_jump_go = handles.nframes_jump_go;
% parent_handles.guidata.seek_behaviors_go = handles.seek_behaviors_go;
handles.NJObj.SetSeekBehaviorsGo(handles.seek_behaviors_go);
allJTypes = handles.NJObj.GetAllTypes;
curJType = allJTypes{get(handles.jumpToPopUp,'Value')};
handles.NJObj.SetCurrentType(curJType);
JLabel('SetJumpGoMenuLabels',parent_handles);

if strcmp(handles.NJObj.GetCurrentType,'Thresholds')
  selFeatures = get(handles.thresholdPopup1,'Value');
  if selFeatures==1
    handles.NJObj.perframeSelFeatures = [];
    handles.NJObj.perframeSelThresholds = [];
    handles.NJObj.perframeComparisonType = [];
  else
    handles.NJObj.perframeSelFeatures = selFeatures-1;
    handles.NJObj.perframeSelThresholds = ...
      str2double(get(handles.thresholdValue1,'String'));
    handles.NJObj.perframeComparisonType = get(handles.thresholdType1,'Value');
  end
end


% --- Executes on selection change in jumpToPopUp.
function jumpToPopUp_Callback(hObject, eventdata, handles)
% hObject    handle to jumpToPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns jumpToPopUp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from jumpToPopUp
contents = cellstr(get(hObject,'String'));
handles.NJObj.SetCurrentType(contents{get(hObject,'Value')});
updateThresholdButtons(handles);
updateConfidenceButtons(handles);

% --- Executes during object creation, after setting all properties.
function jumpToPopUp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to jumpToPopUp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in thresholdPopup1.
function thresholdPopup1_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdPopup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns thresholdPopup1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from thresholdPopup1
handles = guidata(hObject);
selFeatures = get(hObject,'Value');
if selFeatures==1
  handles.NJObj.perframeSelFeatures = [];
  handles.NJObj.perframeSelThresholds = [];
  handles.NJObj.perframeComparisonType = [];
else
  handles.NJObj.perframeSelFeatures = selFeatures-1;
  handles.NJObj.perframeSelThresholds = ...
    str2double(get(handles.thresholdValue1,'Value'));
  handles.NJObj.perframeComparisonType = get(handles.thresholdType1,'Value');
end
   
   
 

% --- Executes during object creation, after setting all properties.
function thresholdPopup1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdPopup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function thresholdValue1_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdValue1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of thresholdValue1 as text
%        str2double(get(hObject,'String')) returns contents of thresholdValue1 as a double
handles = guidata(hObject);
handles.NJObj.perframeSelThresholds = ...
  str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function thresholdValue1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdValue1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function updateThresholdButtons(handles)

if ~strcmp(handles.NJObj.GetCurrentType,'Thresholds on perframe values')
  for ndx = 1
    set(handles.(sprintf('thresholdPopup%d',ndx)),'Enable','off');
    set(handles.(sprintf('thresholdType%d',ndx)),'Enable','off');
    set(handles.(sprintf('thresholdValue%d',ndx)),'Enable','off');
  end
else
  if ~isempty(handles.NJObj.perframeSelFeatures)
    for ndx = 1
      set(handles.(sprintf('thresholdPopup%d',ndx)),'Enable','on',...
        'Value',handles.NJObj.perframeSelFeatures(ndx)+1);
      set(handles.(sprintf('thresholdType%d',ndx)),'Enable','on',...
        'Value',handles.NJObj.perframeComparisonType(ndx));
      set(handles.(sprintf('thresholdValue%d',ndx)),'Enable','on',...
        'String',sprintf('%d',handles.NJObj.perframeSelThresholds(ndx)));
    end
  else
    for ndx = 1
      set(handles.(sprintf('thresholdPopup%d',ndx)),'Enable','on',...
        'Value',1);
      set(handles.(sprintf('thresholdType%d',ndx)),'Enable','on',...
        'Value',1);
      set(handles.(sprintf('thresholdValue%d',ndx)),'Enable','on',...
        'Value',0);
    end
  end    
end

function updateConfidenceButtons(handles)

if ~strcmp(handles.NJObj.GetCurrentType,'Low Confidence') && ~strcmp(handles.NJObj.GetCurrentType,'High Confidence Errors')
  set(handles.text_conf1,'Enable','off');
  set(handles.text_conf2,'Enable','off');
  set(handles.edit_thresh_high,'Enable','off');
  set(handles.edit_thresh_low,'Enable','off');
else
  set(handles.text_conf1,'Enable','on');
  set(handles.text_conf2,'Enable','on');
  hthresh = handles.NJObj.hthresh;
  lthresh = handles.NJObj.lthresh;
  set(handles.edit_thresh_high,'Enable','on','String',sprintf('%.2f',hthresh));
  set(handles.edit_thresh_low,'Enable','on','String',sprintf('%.2f',lthresh));
end


% --- Executes on selection change in thresholdType1.
function thresholdType1_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdType1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns thresholdType1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from thresholdType1
handles = guidata(hObject);
handles.NJObj.perframeComparisonType = get(hObject,'Value');


% --- Executes during object creation, after setting all properties.
function thresholdType1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdType1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_thresh_high_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thresh_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thresh_high as text
%        str2double(get(hObject,'String')) returns contents of edit_thresh_high as a double
val = str2double(get(hObject,'String'));
if isempty(val) || isnan(val) || val<-1 || val>1,
  uiwait(warndlg('Enter a number between -1 and 1'));
  return;
end
if val < handles.NJObj.lthresh,
  uiwait(warndlg('High threshold cannot be smaller than low threshold'));
  return;
end

handles.NJObj.SetHighThresh(val);

% --- Executes during object creation, after setting all properties.
function edit_thresh_high_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thresh_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_thresh_low_Callback(hObject, eventdata, handles)
% hObject    handle to edit_thresh_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_thresh_high as text
%        str2double(get(hObject,'String')) returns contents of edit_thresh_high as a double

val = str2double(get(hObject,'String'));
if isempty(val) || isnan(val) || val<-1 || val>1,
  uiwait(warndlg('Enter a number between -1 and 1'));
  return;
end
if val > handles.NJObj.lthresh,
  uiwait(warndlg('Low threshold cannot be larger than high threshold'));
  return;
end

handles.NJObj.SetLowThresh(val);


% --- Executes during object creation, after setting all properties.
function edit_thresh_low_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_thresh_high (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
