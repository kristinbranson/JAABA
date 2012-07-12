function varargout = JAABA_plot(varargin)
%JAABA_PLOT M-file for JAABA_plot.fig
%      JAABA_PLOT, by itself, creates a new JAABA_PLOT or raises the existing
%      singleton*.
%
%      H = JAABA_PLOT returns the handle to a new JAABA_PLOT or the handle to
%      the existing singleton*.
%
%      JAABA_PLOT('Property','Value',...) creates a new JAABA_PLOT using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to JAABA_plot_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      JAABA_PLOT('CALLBACK') and JAABA_PLOT('CALLBACK',hObject,...) call the
%      local function named CALLBACK in JAABA_PLOT.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JAABA_plot

% Last Modified by GUIDE v2.5 12-Jul-2012 15:57:38

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JAABA_plot_OpeningFcn, ...
                   'gui_OutputFcn',  @JAABA_plot_OutputFcn, ...
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


% --- Executes just before JAABA_plot is made visible.
function JAABA_plot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

if(exist('matlabpool')==2 && matlabpool('size')==0)
  matlabpool open
end

%keyboard
if(exist('JAABA_plot.mat')==2)
  handles_saved=load('JAABA_plot.mat');
  handles_saved=handles_saved.handles;
  handles.experimentlist=handles_saved.experimentlist;
  handles.experimentvalue=handles_saved.experimentvalue;
  handles.experimentlist2=handles_saved.experimentlist2;
  handles.experimentvalue2=handles_saved.experimentvalue2;
  handles.behaviors=handles_saved.behaviors;
  handles.behaviorlist=handles_saved.behaviorlist;
  handles.behaviorvalue=handles_saved.behaviorvalue;
  handles.behaviorlogic=handles_saved.behaviorlogic;
  handles.behaviorvalue2=handles_saved.behaviorvalue2;
  handles.features=handles_saved.features;
  handles.featurelist=handles_saved.featurelist;
  handles.featurevalue=handles_saved.featurevalue;
  handles.timing=handles_saved.timing;
  handles.statistic=handles_saved.statistic;
  handles.windowradius=handles_saved.windowradius;
  handles.bar=handles_saved.bar;
  handles.zoom=handles_saved.zoom;
  handles.sex=handles_saved.sex;
  set(handles.ExperimentList,'String',handles.experimentlist);
  set(handles.ExperimentList,'Value',handles.experimentvalue);
  set(handles.ExperimentList2,'String',handles.experimentlist2);
  set(handles.ExperimentList2,'Value',handles.experimentvalue2);
  set(handles.BehaviorList,'String',handles.behaviorlist);
  set(handles.BehaviorList,'Value',handles.behaviorvalue);
  set(handles.BehaviorLogic,'Value',handles.behaviorlogic);
  set(handles.BehaviorList2,'String',handles.behaviorlist);
  set(handles.BehaviorList2,'Value',handles.behaviorvalue2);
  set(handles.FeatureList,'String',handles.featurelist);
  set(handles.FeatureList,'Value',handles.featurevalue);
  if(length(handles.experimentlist)==0)
    set(handles.ExperimentList,'enable','off');
  end
  if(length(handles.experimentlist2)==0)
    set(handles.ExperimentList2,'enable','off');
  end
  if((length(handles.experimentlist)==0) && (length(handles.experimentlist2)==0))
    set(handles.BehaviorList,'enable','off');
    set(handles.BehaviorLogic,'enable','off');
    set(handles.BehaviorList2,'enable','off');
    set(handles.FeatureList,'enable','off');
    set(handles.SexList,'enable','off');
  end
  if(handles.behaviorlogic==1)
    set(handles.BehaviorList2,'enable','off');
  else
    set(handles.BehaviorList2,'enable','on');
  end
  switch(handles.timing)
    case 1
      MenuTimeSeriesPrefsTimingOnset_Callback(hObject, eventdata, handles);
    case 2
      MenuTimeSeriesPrefsTimingOffset_Callback(hObject, eventdata, handles);
  end
  switch(handles.statistic)
    case 1
      MenuTimeSeriesPrefsStatisticMean_Callback(hObject, eventdata, handles);
    case 2
      MenuTimeSeriesPrefsStatisticMedian_Callback(hObject, eventdata, handles);
  end
  switch(handles.bar)
    case 1
      MenuTimeSeriesPrefsBarsNone_Callback(hObject, eventdata, handles);
    case 2
      MenuTimeSeriesPrefsBarsDeviation_Callback(hObject, eventdata, handles);
    case 3
      MenuTimeSeriesPrefsBarsError_Callback(hObject, eventdata, handles);
  end
  switch(handles.zoom)
    case 1
      MenuTimeSeriesPrefsZoomBlack_Callback(hObject, eventdata, handles);
    case 2
      MenuTimeSeriesPrefsZoomRed_Callback(hObject, eventdata, handles);
  end
  set(handles.SexList,'Value',handles.sex);
else
  handles.experimentlist={};
  handles.experimentvalue=[];
  handles.experimentlist2={};
  handles.experimentvalue2=[];
  handles.behaviors={};
  handles.behaviorlist={''};
  handles.behaviorvalue=1;
  handles.behaviorlogic=1;
  handles.behaviorvalue2=1;
  handles.features={};
  handles.featurelist={''};
  handles.featurevalue=1;
  handles.timing=1;
  handles.statistic=1;
  handles.windowradius=10;
  handles.bar=1;
  handles.zoom=1;
  handles.sex=1;
  set(handles.ExperimentList,'enable','off');
  set(handles.ExperimentList2,'enable','off');
  set(handles.BehaviorList,'enable','off');
  set(handles.BehaviorLogic,'enable','off');
  set(handles.BehaviorList2,'enable','off');
  set(handles.FeatureList,'enable','off');
  set(handles.SexList,'enable','off');
end

% Choose default command line output for JAABA_plot
handles.output = hObject;

set(hObject,'CloseRequestFcn',@figure_CloseRequestFcn);
set(handles.Table,'CellSelectionCallback',@Interesting_CellSelectionCallback);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JAABA_plot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Executes on button press in Reset.
function Reset_Callback(hObject, eventdata, handles)
% hObject    handle to Reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

handles.experimentlist={};
handles.experimentvalue=[];
handles.experimentlist2={};
handles.experimentvalue2=[];
handles.behaviors={};
handles.behaviorlist={''};
handles.behaviorvalue=1;
handles.behaviorlogic=1;
handles.behaviorvalue2=1;
handles.features={};
handles.featurelist={''};
handles.featurevalue=1;
handles.timing=1;
handles.statistic=1;
handles.windowradius=10;
handles.bar=1;
handles.zoom=1;
handles.sex=1;

set(handles.ExperimentList,'String',handles.experimentlist,'Value',handles.experimentvalue,'enable','off');
set(handles.ExperimentList2,'String',handles.experimentlist2,'Value',handles.experimentvalue2,'enable','off');
set(handles.BehaviorList,'String',handles.behaviorlist,'enable','off');
set(handles.BehaviorLogic,'Value',handles.behaviorlogic,'enable','off');
set(handles.BehaviorList2,'String',handles.behaviorlist,'enable','off');
set(handles.FeatureList,'String',handles.featurelist,'enable','off');
MenuTimeSeriesPrefsTimingOnset_Callback(hObject, eventdata, handles);
MenuTimeSeriesPrefsStatisticMean_Callback(hObject, eventdata, handles);
MenuTimeSeriesPrefsBarsNone_Callback(hObject, eventdata, handles);
MenuTimeSeriesPrefsZoomBlack_Callback(hObject, eventdata, handles);
set(handles.SexList,'Value',handles.sex);

save('JAABA_plot.mat','handles');

guidata(hObject, handles);


% ---
function figure_CloseRequestFcn(hObject, eventdata)

handles=guidata(hObject);
save('JAABA_plot.mat','handles');
delete(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = JAABA_plot_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in UpdatePlot.
function UpdatePlot_Callback(hObject, eventdata, handles)
% hObject    handle to UpdatePlot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)


% --- Executes on selection change in FeatureList.
function FeatureList_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns FeatureList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from FeatureList

handles.featurevalue=get(handles.FeatureList,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function FeatureList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to FeatureList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BehaviorList.
function BehaviorList_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorList

handles.behaviorvalue=get(handles.BehaviorList,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function BehaviorList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BehaviorLogic.
function BehaviorLogic_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorLogic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorLogic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorLogic

handles.behaviorlogic=get(handles.BehaviorLogic,'Value');
if(handles.behaviorlogic==1)
  set(handles.BehaviorList2,'enable','off');
else
  set(handles.BehaviorList2,'enable','on');
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function BehaviorLogic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorLogic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in BehaviorList2.
function BehaviorList2_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns BehaviorList2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from BehaviorList2


% --- Executes during object creation, after setting all properties.
function BehaviorList2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BehaviorList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ExperimentList2.
function ExperimentList2_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ExperimentList2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ExperimentList2

handles.experimentvalue2=get(handles.ExperimentList2,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ExperimentList2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExperimentList2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ExperimentList.
function ExperimentList_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ExperimentList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ExperimentList

handles.experimentvalue=get(handles.ExperimentList,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ExperimentList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExperimentList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ---
function ret_val=check_for_diff_and_return_intersection(arg)

if(length(arg)==1)  ret_val=arg{1};  return;  end

flatten=[arg{:}];
anotb=setdiff(flatten,arg{end});
bnota=setdiff(arg{end},flatten);
if(numel(anotb)>0)
  tmp=char(anotb);
  tmp=[tmp repmat(', ',size(tmp,1),1)];
  tmp=reshape(tmp',1,numel(tmp));
  tmp=tmp(1:end-2);
  warning([tmp ' are/is in prior experiments but not new one.  proceeding with intersection.']);
end
if(numel(bnota)>0)
  tmp=char(bnota);
  tmp=[tmp repmat(', ',size(tmp,1),1)];
  tmp=reshape(tmp',1,numel(tmp));
  tmp=tmp(1:end-2);
  warning([tmp ' are/is in new experiment but not prior ones.  proceeding with intersection.']);
end
ret_val=intersect(arg{end},flatten);


% ---
function experiment_add(hObject,handles,arg)

set(handles.Status,'string','Thinking...');

tmp=uigetdir([],'Select Experiment Directory');
if(tmp==0)  return;  end
if(arg==1)
  handles.experimentlist{end+1}=tmp;
  handles.experimentvalue=length(handles.experimentlist);
  set(handles.ExperimentList,'String',handles.experimentlist);
  set(handles.ExperimentList,'Value',handles.experimentvalue);
  set(handles.ExperimentList,'enable','on');
  ptr=handles.experimentlist;
else
  handles.experimentlist2{end+1}=tmp;
  handles.experimentvalue2=length(handles.experimentlist2);
  set(handles.ExperimentList2,'String',handles.experimentlist2);
  set(handles.ExperimentList2,'Value',handles.experimentvalue2);
  set(handles.ExperimentList2,'enable','on');
  ptr=handles.experimentlist2;
end

set(handles.BehaviorList,'enable','on');
set(handles.BehaviorLogic,'enable','on');
if(handles.behaviorlogic>1)
  set(handles.BehaviorList2,'enable','on');
end
set(handles.FeatureList,'enable','on');
set(handles.SexList,'enable','on');

tmp=dir(fullfile(ptr{end},'scores*.mat'));
[handles.behaviors{length(ptr)+(arg-1)*length(handles.experimentlist)}{1:length(tmp)}]=deal(tmp.name);
handles.behaviorlist=check_for_diff_and_return_intersection(handles.behaviors);
set(handles.BehaviorList,'String',handles.behaviorlist);
set(handles.BehaviorList2,'String',handles.behaviorlist);

tmp=dir(fullfile(ptr{end},'perframe','*.mat'));
[handles.features{length(ptr)+(arg-1)*length(handles.experimentlist)}{1:length(tmp)}]=deal(tmp.name);
handles.featurelist=check_for_diff_and_return_intersection(handles.features);
set(handles.FeatureList,'String',handles.featurelist);

guidata(hObject,handles);

set(handles.Status,'string','Ready.');


% --- Executes on button press in ExperimentAdd2.
function ExperimentAdd2_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentAdd2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

experiment_add(hObject,handles,2);


% --- Executes on button press in ExperimentAdd.
function ExperimentAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

experiment_add(hObject,handles,1);


% ---
function experiment_delete(hObject,handles,idx)

handles.behaviors(idx)=[];
handles.behaviorlist=unique([handles.behaviors{:}]);
handles.behaviorvalue=max(1,min(handles.behaviorvalue,length(handles.behaviorlist)));
handles.behaviorvalue2=max(1,min(handles.behaviorvalue2,length(handles.behaviorlist)));
if(length(handles.behaviorlist)==0)
  handles.behaviorlist={''};
end
set(handles.BehaviorList,'String',handles.behaviorlist);
set(handles.BehaviorList2,'String',handles.behaviorlist);
set(handles.BehaviorList,'Value',handles.behaviorvalue);
set(handles.BehaviorList2,'Value',handles.behaviorvalue2);

handles.features(idx)=[];
handles.featurelist=unique([handles.features{:}]);
handles.featurevalue=max(1,min(handles.featurevalue,length(handles.featurelist)));
if(length(handles.featurelist)==0)
  handles.featurelist={''};
end
set(handles.FeatureList,'String',handles.featurelist);
set(handles.FeatureList,'Value',handles.featurevalue);

guidata(hObject,handles);

set(handles.Status,'string','Ready.');


% --- Executes on button press in ExperimentDelete2.
function ExperimentDelete2_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentDelete2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(length(handles.experimentlist2)==0)  return;  end

set(handles.Status,'string','Thinking...');

idx=get(handles.ExperimentList2,'Value');
handles.experimentlist2(idx)=[];
handles.experimentvalue2=min(handles.experimentvalue2,length(handles.experimentlist2));
handles.experimentvalue2=max(handles.experimentvalue2,1);
if(length(handles.experimentlist2)==0)
  handles.experimentlist2={};
  handles.experimentvalue2=[];
  set(handles.ExperimentList2,'enable','off');
  if(length(handles.experimentlist)==0)
    set(handles.BehaviorList,'enable','off');
    set(handles.BehaviorLogic,'enable','off');
    set(handles.BehaviorList2,'enable','off');
    set(handles.FeatureList,'enable','off');
    set(handles.SexList,'enable','off');
  end
end
set(handles.ExperimentList2,'String',handles.experimentlist2);
set(handles.ExperimentList2,'Value',handles.experimentvalue2);

experiment_delete(hObject,handles,idx+length(handles.experimentlist));


% --- Executes on button press in ExperimentDelete.
function ExperimentDelete_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

if(length(handles.experimentlist)==0)  return;  end

set(handles.Status,'string','Thinking...');

idx=get(handles.ExperimentList,'Value');
handles.experimentlist(idx)=[];
handles.experimentvalue=min(handles.experimentvalue,length(handles.experimentlist));
handles.experimentvalue=max(handles.experimentvalue,1);
if(length(handles.experimentlist)==0)
  handles.experimentlist={};
  handles.experimentvalue=[];
  set(handles.ExperimentList,'enable','off');
  if(length(handles.experimentlist2)==0)
    set(handles.BehaviorList,'enable','off');
    set(handles.BehaviorLogic,'enable','off');
    set(handles.BehaviorList2,'enable','off');
    set(handles.FeatureList,'enable','off');
    set(handles.SexList,'enable','off');
  end
end
set(handles.ExperimentList,'String',handles.experimentlist);
set(handles.ExperimentList,'Value',handles.experimentvalue);

experiment_delete(hObject,handles,idx);


% --- Executes on button press in ExperimentMove2.
function ExperimentMove2_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentMove2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(length(handles.experimentlist2)==0)  return;  end

set(handles.Status,'string','Thinking...');

idx=get(handles.ExperimentList2,'Value');
idx2=length(handles.experimentlist);
handles.experimentlist={handles.experimentlist{:} handles.experimentlist2{idx}};
tmp=length(handles.experimentlist);
handles.experimentvalue=(tmp-length(idx)+1):tmp;
handles.experimentlist2(idx)=[];
handles.experimentvalue2=1;
if(length(handles.experimentlist2)==0)
  handles.experimentlist2={};
  handles.experimentvalue2=[];
  set(handles.ExperimentList2,'enable','off');
end
set(handles.ExperimentList,'String',handles.experimentlist);
set(handles.ExperimentList,'Value',handles.experimentvalue,'enable','on');
set(handles.ExperimentList2,'String',handles.experimentlist2);
set(handles.ExperimentList2,'Value',handles.experimentvalue2);
handles.behaviors=handles.behaviors([1:idx2 idx2+idx setdiff((idx2+1):length(handles.behaviors),idx2+idx)]);
handles.features=handles.features([1:idx2 idx2+idx setdiff((idx2+1):length(handles.behaviors),idx2+idx)]);

guidata(hObject,handles);

set(handles.Status,'string','Ready.');


% --- Executes on button press in ExperimentMove.
function ExperimentMove_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentMove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(length(handles.experimentlist)==0)  return;  end

set(handles.Status,'string','Thinking...');

idx=get(handles.ExperimentList,'Value');
handles.experimentlist2={handles.experimentlist2{:} handles.experimentlist{idx}};
tmp=length(handles.experimentlist2);
handles.experimentvalue2=(tmp-length(idx)+1):tmp;
handles.experimentlist(idx)=[];
handles.experimentvalue=1;
if(length(handles.experimentlist)==0)
  handles.experimentlist={};
  handles.experimentvalue=[];
  set(handles.ExperimentList,'enable','off');
end
set(handles.ExperimentList2,'String',handles.experimentlist2);
set(handles.ExperimentList2,'Value',handles.experimentvalue2,'enable','on');
set(handles.ExperimentList,'String',handles.experimentlist);
set(handles.ExperimentList,'Value',handles.experimentvalue);
handles.behaviors=handles.behaviors([setdiff(1:length(handles.behaviors),idx) idx]);
handles.features=handles.features([setdiff(1:length(handles.features),idx) idx]);

guidata(hObject,handles);

set(handles.Status,'string','Ready.');


% --- Executes on selection change in SexList.
function SexList_Callback(hObject, eventdata, handles)
% hObject    handle to SexList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns SexList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SexList


% --- Executes during object creation, after setting all properties.
function SexList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SexList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ---
function ret_val=get_label(feature_name,feature_units)

ret_val=[strrep(char(feature_name),'_','-') ' ('];
if(~isempty(feature_units.num))
  if(length(feature_units.num)>1)
    ret_val=[ret_val '['];
  end
  tmp=char(feature_units.num);
  tmp=[tmp repmat('-',size(tmp,1),1)];
  tmp=reshape(tmp',1,numel(tmp));
  tmp=tmp(1:end-1);
  ret_val=[ret_val tmp];
  if(length(feature_units.num)>1)
    ret_val=[ret_val ']'];
  end
  if(~isempty(feature_units.den))
    ret_val=[ret_val ' / ' ];
  end
end
if(~isempty(feature_units.den))
  if(length(feature_units.den)>1)
    ret_val=[ret_val '['];
  end
  tmp=char(feature_units.den);
  tmp=[tmp repmat('-',size(tmp,1),1)];
  tmp=reshape(tmp',1,numel(tmp));
  tmp=tmp(1:end-1);
  ret_val=[ret_val tmp];
  if(length(feature_units.den)>1)
    ret_val=[ret_val ']'];
  end
  if(isempty(feature_units.num))
    ret_val=[ret_val ' ^ -1 ' ];
  end
end
ret_val=[ret_val ')'];


% --- 
function [during not_during]=calculate_histogram(behavior_data,behavior_logic,behavior_data2,feature_data,sex_data)

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)~=2)  error('uhoh');  end
  for i=1:length(feature_data.data)
    feature_data.data{i}=strcmp(feature_data.data{i},vals{1});
  end
end

during={};  not_during={};
%partition_idx=cell(1,length(feature_data.data));
partition_idx=cell(1,length(behavior_data.allScores.t0s));
for i=1:length(behavior_data.allScores.t0s)
  tmp1=false(1,length(feature_data.data{i}));
  for j=1:length(behavior_data.allScores.t0s{i})
    tmp1((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
        -behavior_data.allScores.tStart(i)+1)=true;
  end
  if(behavior_logic>1)
    tmp2=false(1,length(feature_data.data{i}));
    for j=1:length(behavior_data2.allScores.t0s{i})
      tmp2((behavior_data2.allScores.t0s{i}(j):(behavior_data2.allScores.t1s{i}(j)-1))...
          -behavior_data2.allScores.tStart(i)+1)=true;
    end
  end
  switch(behavior_logic)
    case(1)
      partition_idx{i}=tmp1;
    case(2)
      partition_idx{i}=tmp1 & tmp2;
    case(3)
      partition_idx{i}=tmp1 & ~tmp2;
    case(4)
      partition_idx{i}=tmp1 | tmp2;
  end
  during{i}=feature_data.data{i}(partition_idx{i} & sex_data.data{i}(1:length(partition_idx{i})));
  not_during{i}=feature_data.data{i}((~partition_idx{i}) & sex_data.data{i}(1:length(partition_idx{i})));
end

during=[during{:}];
not_during=[not_during{:}];


% ---
function feature_units=plot_histogram(experiment_value,experiment_list,...
    behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
    feature_value,feature_list,sex,color)

during={};  not_during={};  feature_units={};
parfor j=1:length(experiment_value)
  behavior_data=load(fullfile(char(experiment_list(experiment_value(j))),...
        char(behavior_list(behavior_value))));
  if(behavior_logic>1)
    behavior_data2=load(fullfile(char(experiment_list(experiment_value(j))),...
        char(behavior_list2(behavior_value2))));
  else
    behavior_data2=[];
  end
  feature_data=load(fullfile(char(experiment_list(experiment_value(j))),'perframe',...
        char(feature_list(feature_value))));
  feature_units{j}=feature_data.units;

  sex_data=load(fullfile(char(experiment_list(experiment_value(j))),'perframe',...
      char(feature_list(find(strcmp(feature_list,'sex.mat'))))));
  for i=1:length(sex_data.data)
    switch(sex)
      case(1)
        sex_data.data{i}=ones(1,length(sex_data.data{i}));
      case(2)
        sex_data.data{i}=strcmp(sex_data.data{i},'M');
      case(3)
        sex_data.data{i}=strcmp(sex_data.data{i},'F');
    end
  end

  [during{j} not_during{j}]=calculate_histogram(behavior_data,behavior_logic,behavior_data2,...
      feature_data,sex_data);
end
during=[during{:}];
not_during=[not_during{:}];
feature_units=feature_units{1};

tmp=linspace(min([during not_during]),max([during not_during]));
hist_during=hist(during,tmp);
hist_not_during=hist(not_during,tmp);
plot(tmp,hist_not_during./max(hist_not_during),[color '-']);
plot(tmp,hist_during./max(hist_during),[color '-'],'linewidth',3);


% --- Executes on button press in PlotHistogram.
function PlotHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to PlotHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...');  drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_value=get(handles.BehaviorList,'Value');
behavior_list=get(handles.BehaviorList,'String');
behavior_logic=get(handles.BehaviorLogic,'Value');
behavior_value2=get(handles.BehaviorList2,'Value');
behavior_list2=get(handles.BehaviorList2,'String');
feature_value=get(handles.FeatureList,'Value');
feature_list=get(handles.FeatureList,'String');
sex=get(handles.SexList,'Value');

axes(handles.Axes);
cla;  hold on;

if(length(experiment_value)>0)
  feature_units=plot_histogram(experiment_value,experiment_list,...
      behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      feature_value,feature_list,sex,'r');
end
if(length(experiment_value2)>0)
  feature_units=plot_histogram(experiment_value2,experiment_list2,...
      behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      feature_value,feature_list,sex,'b');
end

%legend('during','not during');
%text(tmp(5),0.95,['during = ' num2str(100*length(during)/(length(during)+length(not_during)),3) '%'],'color','b');
xlabel(get_label(feature_list(feature_value),feature_units));
ylabel('normalized');
axis tight

set(handles.Status,'string','Ready.');


% ---
function table_data=calculate_interesting_histograms(experiment_value,experiment_list,...
    behavior_list,feature_list)

table_data=zeros(length(behavior_list),length(feature_list),6);
parfor b=1:length(behavior_list)
  k=1;
  parfor_tmp=zeros(length(feature_list),6);
  for f=1:length(feature_list)
    during={};  not_during={};
    for e=1:length(experiment_value)
      behavior_data=load(fullfile(char(experiment_list(experiment_value(e))),char(behavior_list(b))));
      feature_data=load(fullfile(char(experiment_list(experiment_value(e))),'perframe',char(feature_list(f))));
      sex_data=[];
      for i=1:length(feature_data.data)
        sex_data.data{i}=ones(1,length(feature_data.data{i}));
      end
      [during{e} not_during{e}]=calculate_histogram(behavior_data,1,[],feature_data,sex_data);
    end
    during=[during{:}];
    not_during=[not_during{:}];
    %parfor_tmp(k,:)=[b f (mean(during)-mean(not_during))/sqrt((std(during)^2+std(not_during)^2)/2)];
    parfor_tmp(k,:)=[b f mean(during) mean(not_during) std(during) std(not_during)];
    k=k+1;
  end
  table_data(b,:,:)=parfor_tmp;
  disp([num2str(b) ' of ' num2str(length(behavior_list))]);
end
table_data=reshape(table_data,prod(size(table_data))/6,6);


% --- Executes on button press in InterestingHistograms.
function InterestingHistograms_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingHistograms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...');  drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_list=get(handles.BehaviorList,'String');
feature_list=get(handles.FeatureList,'String');

if(length(experiment_value)>0)
  table_data=calculate_interesting_histograms(experiment_value,experiment_list,...
      behavior_list,feature_list);
end
if(length(experiment_value2)>0)
  table_data2=calculate_interesting_histograms(experiment_value2,experiment_list2,...
      behavior_list,feature_list);
end

if((length(experiment_value)>0) && (length(experiment_value2)>0))
  tmp2=[table_data(:,1:2) (table_data2(:,3)-table_data(:,3))./sqrt(table_data2(:,5).^2+table_data(:,5).^2)];
elseif(length(experiment_value)>0)
  tmp2=[table_data(:,1:2) (table_data(:,3)-table_data(:,4))./sqrt(table_data(:,5).^2+table_data(:,6).^2)];
elseif(length(experiment_value2)>0)
  tmp2=[table_data2(:,1:2) (table_data2(:,3)-table_data2(:,4))./sqrt(table_data2(:,5).^2+table_data2(:,6).^2)];
end

tmp2=sortrows(tmp2,-3);

tmp=cell(size(tmp2));
tmp(:,1)=behavior_list(tmp2(:,1));
tmp(:,2)=feature_list(tmp2(:,2));
tmp(:,3)=num2cell(tmp2(:,3));
set(handles.Table,'Data',tmp);
set(handles.Table,'ColumnName',{'Behavior' 'Feature' 'd'''});
set(handles.Table,'ColumnWidth',{150 150 50});
set(handles.Table,'RowStriping','on','BackgroundColor',[1 1 1; 0.95 0.95 0.95]);

handles.table_data=tmp2;
handles.table='histogram';
guidata(hObject,handles);

set(handles.Status,'string','Ready.');


% --- 
function triggered_data=calculate_timeseries(behavior_data,behavior_logic,behavior_data2,feature_data,timing,statistic,windowradius)

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)~=2)  error('uhoh');  end
  for i=1:length(feature_data.data)
    feature_data.data{i}=strcmp(feature_data.data{i},vals{1});
  end
end

k=1;
triggered_data=zeros(length([behavior_data.allScores.t0s{:}]),1+2*windowradius);
for i=1:length(behavior_data.allScores.t0s)
  feature_data_padded=[nan(1,windowradius) feature_data.data{i} nan(1,windowradius)];
  if(behavior_logic>1)
    idx2=[behavior_data2.allScores.t0s{i}'-behavior_data2.allScores.tStart(i) ...
       behavior_data2.allScores.t1s{i}'-behavior_data2.allScores.tStart(i)];
  end
  for j=1:length(behavior_data.allScores.t0s{i})
    switch(timing)
      case(1)
        idx=behavior_data.allScores.t0s{i}(j)-behavior_data.allScores.tStart(i);
      case(2)
        idx=behavior_data.allScores.t1s{i}(j)-behavior_data.allScores.tStart(i);
    end
    if((behavior_logic==2) && (diff(sum(idx2<idx))==0))
      continue;
    end
    if((behavior_logic==3) && (diff(sum(idx2<idx))==-1))
      continue;
    end
    triggered_data(k,:)=feature_data_padded(idx+(0:(2*windowradius))+(2-timing));
    k=k+1;
  end
end


%---
function [feature_units range h]=plot_timeseries(experiment_value,experiment_list,...
    behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
    feature_value,feature_list,sex,timing,statistic,windowradius,bar,color)

triggered_data=[];  feature_units={};
parfor j=1:length(experiment_value)
  behavior_data=load(fullfile(char(experiment_list(experiment_value(j))),...
      char(behavior_list(behavior_value))));
  if(behavior_logic>1)
    behavior_data2=load(fullfile(char(experiment_list(experiment_value(j))),...
        char(behavior_list2(behavior_value2))));
  else
    behavior_data2=[];
  end
  feature_data=load(fullfile(char(experiment_list(experiment_value(j))),'perframe',...
      char(feature_list(feature_value))));
  feature_units{j}=feature_data.units;

  tmp=calculate_timeseries(behavior_data,behavior_logic,behavior_data2,...
      feature_data,timing,statistic,windowradius);
  triggered_data=[triggered_data; tmp];
end

feature_units=feature_units{1};

k=size(triggered_data,1);
switch(statistic)
  case(1)
    triggered_data(k+3,:)=nanmean(triggered_data(1:k,:),1);
    tmp=nanstd(triggered_data(1:k,:),1)./sqrt(size(triggered_data,1));
    triggered_data(k+2,:)=triggered_data(k+3,:)+tmp;
    triggered_data(k+4,:)=triggered_data(k+3,:)-tmp;
    tmp=nanstd(triggered_data(1:k,:),1);
    triggered_data(k+1,:)=triggered_data(k+3,:)+tmp;
    triggered_data(k+5,:)=triggered_data(k+3,:)-tmp;
  case(2)
    triggered_data(k+1,:)=prctile(triggered_data(1:k,:),95,1);
    triggered_data(k+2,:)=prctile(triggered_data(1:k,:),75,1);
    triggered_data(k+3,:)=nanmedian(triggered_data(1:k,:),1);
    triggered_data(k+4,:)=prctile(triggered_data(1:k,:),25,1);
    triggered_data(k+5,:)=prctile(triggered_data(1:k,:), 5,1);
end

plot(-windowradius:windowradius,triggered_data(1:end-5,:)','k-');
h(1)=plot(-windowradius:windowradius,triggered_data(end-2,:)',[color '-'],'linewidth',3);
range=[min(triggered_data(end-2,:)) max(triggered_data(end-2,:))];
switch(bar)
  case(2)
    h(2)=plot(-windowradius:windowradius,triggered_data(end-4,:)',[color '-'],'linewidth',1);
    h(3)=plot(-windowradius:windowradius,triggered_data(end-0,:)',[color '-'],'linewidth',1);
    range=[range; min(triggered_data(end-4,:)) max(triggered_data(end-4,:))];
    range=[range; min(triggered_data(end-0,:)) max(triggered_data(end-0,:))];
  case(3)
    h(2)=plot(-windowradius:windowradius,triggered_data(end-3,:)',[color '-'],'linewidth',1);
    h(3)=plot(-windowradius:windowradius,triggered_data(end-1,:)',[color '-'],'linewidth',1);
    range=[range; min(triggered_data(end-3,:)) max(triggered_data(end-3,:))];
    range=[range; min(triggered_data(end-1,:)) max(triggered_data(end-1,:))];
end


% --- Executes on button press in PlotTimeSeries.
function PlotTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...');  drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_value=get(handles.BehaviorList,'Value');
behavior_list=get(handles.BehaviorList,'String');
behavior_logic=get(handles.BehaviorLogic,'Value');
behavior_value2=get(handles.BehaviorList2,'Value');
behavior_list2=get(handles.BehaviorList2,'String');
feature_value=get(handles.FeatureList,'Value');
feature_list=get(handles.FeatureList,'String');
sex=get(handles.SexList,'Value');
timing=handles.timing;
statistic=handles.statistic;
windowradius=handles.windowradius;

axes(handles.Axes);
cla;  hold on;

range=[];
if(length(experiment_value)>0)
  [feature_units tmp h]=plot_timeseries(experiment_value,experiment_list,...
      behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      feature_value,feature_list,sex,timing,statistic,windowradius,handles.bar,'r');
  range=[range; tmp];
end
if(length(experiment_value2)>0)
  [feature_units tmp h2]=plot_timeseries(experiment_value2,experiment_list2,...
      behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      feature_value,feature_list,sex,timing,statistic,windowradius,handles.bar,'b');
  range=[range; tmp];
  if(length(experiment_value)>0)
    uistack(h,'top');
    uistack(h2,'top');
  end
end

xlabel('time (frames)');
ylabel(get_label(feature_list(feature_value),feature_units));
axis tight
if(handles.zoom==2)
  v=axis;  axis([v(1) v(2) min(range(:,1)) max(range(:,2))]);
end

set(handles.Status,'string','Ready.');


% ---
function table_data=calculate_interesting_timeseries(experiment_value,experiment_list,...
    behavior_list,feature_list,statistic,windowradius)

table_data=cell(length(behavior_list),length(feature_list),2);
parfor b=1:length(behavior_list)
  parfor_tmp=cell(length(feature_list),2);
  for f=1:length(feature_list)
    for t=1:2  % timing
      parfor_tmp{f,t}=[];
      for e=1:length(experiment_value)
        behavior_data=load(fullfile(char(experiment_list(experiment_value(e))),char(behavior_list(b))));
        feature_data=load(fullfile(char(experiment_list(experiment_value(e))),'perframe',char(feature_list(f))));

        tmp=calculate_timeseries(behavior_data,1,[],feature_data,t,statistic,windowradius);
        parfor_tmp{f,t}=[parfor_tmp{f,t}; tmp];
      end
      foo=nanmean(parfor_tmp{f,t}(:,1:windowradius),2);
      foo=parfor_tmp{f,t}-repmat(foo,1,size(parfor_tmp{f,t},2));
      parfor_tmp{f,t}=[...
          sqrt(nanmean(parfor_tmp{f,t}(:,1:windowradius).^2,2))...
          sqrt(nanmean(parfor_tmp{f,t}(:,(windowradius+1):end).^2,2))];
    end
  end
  table_data(b,:,:)=parfor_tmp;
  disp([num2str(b) ' of ' num2str(length(behavior_list))]);
end
%table_data=reshape(table_data,prod(size(table_data))/5,5);


% --- Executes on button press in InterestingTimeSeries.
function InterestingTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...');  drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_list=get(handles.BehaviorList,'String');
feature_list=get(handles.FeatureList,'String');
statistic=handles.statistic;
windowradius=handles.windowradius;

if(length(experiment_value)>0)
  table_data=calculate_interesting_timeseries(experiment_value,experiment_list,...
      behavior_list,feature_list,statistic,windowradius);
end
if(length(experiment_value2)>0)
  table_data2=calculate_interesting_timeseries(experiment_value2,experiment_list2,...
      behavior_list,feature_list,statistic,windowradius);
end

if((length(experiment_value)>0) && (length(experiment_value2)>0))
  tmp=(cellfun(@(x) nanmean(x(:,2)),table_data2) - cellfun(@(x) nanmean(x(:,2)),table_data)) ./ ...
      sqrt(cellfun(@(x) nanstd(x(:,2)),table_data2).^2 + cellfun(@(x) nanstd(x(:,2)),table_data).^2);
elseif(length(experiment_value)>0)
  tmp=cellfun(@(x) (nanmean(x(:,2))-nanmean(x(:,1)))/sqrt(nanstd(x(:,2)).^2+nanstd(x(:,1)).^2),table_data);
elseif(length(experiment_value2)>0)
  tmp=cellfun(@(x) (nanmean(x(:,2))-nanmean(x(:,1)))/sqrt(nanstd(x(:,2)).^2+nanstd(x(:,1)).^2),table_data2);
end
[tmp2(:,1) tmp2(:,2) tmp2(:,3)]=ind2sub(size(tmp),1:prod(size(tmp)));
tmp2(:,4)=tmp(1:prod(size(tmp)));
tmp2=sortrows(tmp2,-4);

tmp=cell(size(tmp2));
tmp(:,1)=behavior_list(tmp2(:,1));
tmp(:,2)=feature_list(tmp2(:,2));
tmp(:,3)=repmat({'Onset'},size(tmp2,1),1);
  find(tmp2(:,3)==2);
  tmp(ans,3)=repmat({'Offset'},size(ans),1);
tmp(:,4)=num2cell(tmp2(:,4));

set(handles.Table,'Data',tmp);
set(handles.Table,'ColumnName',{'Behavior' 'Feature' 'Timing' 'dRMS'});
set(handles.Table,'ColumnWidth',{150 150 50 50 50});
set(handles.Table,'RowStriping','on','BackgroundColor',[1 1 1; 0.95 0.95 0.95]);

handles.table_data=tmp2;
handles.table='timeseries';
guidata(hObject,handles);

set(handles.Status,'string','Ready.');


% --- 
function [male_count female_count male_total female_total individual_count]=...
    calculate_behaviorstats2(behavior_data,behavior_logic,behavior_data2,sex_data)

male_count=0;  female_count=0;  individual_count=[];
male_total=0;  female_total=0;
%partition_idx=cell(1,length(behavior_data.allScores.t0s));
for i=1:length(behavior_data.allScores.t0s)
  tmp1=false(1,length(sex_data.data{i}));
  for j=1:length(behavior_data.allScores.t0s{i})
    tmp1((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
        -behavior_data.allScores.tStart(i)+1)=true;
  end
  if(behavior_logic>1)
    tmp2=false(1,length(sex_data.data{i}));
    for j=1:length(behavior_data2.allScores.t0s{i})
      tmp2((behavior_data2.allScores.t0s{i}(j):(behavior_data2.allScores.t1s{i}(j)-1))...
          -behavior_data2.allScores.tStart(i)+1)=true;
    end
  end
  switch(behavior_logic)
    case(1)
      partition_idx=tmp1;
    case(2)
      partition_idx=tmp1 & tmp2;
    case(3)
      partition_idx=tmp1 & ~tmp2;
    case(4)
      partition_idx=tmp1 | tmp2;
  end
  male_count=male_count+sum(partition_idx & sex_data.data{i}(1:length(partition_idx)));
  female_count=female_count+sum(partition_idx & (~sex_data.data{i}(1:length(partition_idx))));
  male_total=male_total+sum(sex_data.data{i}(1:length(partition_idx)));
  female_total=female_total+sum(~sex_data.data{i}(1:length(partition_idx)));
  individual_count(i)=100*sum(partition_idx)./length(partition_idx);
end


% ---
function ret_val=calculate_behaviorstats(experiment_value,experiment_list,...
    behavior_list,behavior_logic,behavior_value2,behavior_list2,feature_list)

table_data=cell(1,length(experiment_value));
parfor e=1:length(experiment_value)
  sex_data=load(fullfile(char(experiment_list(experiment_value(e))),'perframe',...
      char(feature_list(find(strcmp(feature_list,'sex.mat'))))));
  for i=1:length(sex_data.data)
    sex_data.data{i}=strcmp(sex_data.data{i},'M');
  end
  behavior_data=load(fullfile(char(experiment_list(experiment_value(e))),char(behavior_list(1))));
  table_data{e}=nan(length(behavior_list),4+length(behavior_data.allScores.t0s));
  if(behavior_logic>1)
    behavior_data2=load(fullfile(char(experiment_list(experiment_value(e))),...
        char(behavior_list2(behavior_value2))));
  else
    behavior_data2=[];
  end
  for b=1:length(behavior_list)
    behavior_data=load(fullfile(char(experiment_list(experiment_value(e))),char(behavior_list(b))));

    [male_count female_count male_total female_total individual_count]=...
        calculate_behaviorstats2(behavior_data,behavior_logic,behavior_data2,sex_data);
    if((4+length(individual_count))>size(table_data{e},2))
      table_data{e}(:,(end+1):(4+length(individual_count)))=...
          nan(size(table_data{e},1):((4+length(individual_count))-size(table_data{e},2)));
    end
    table_data{e}(b,:)=[male_count female_count male_total female_total individual_count];
  end
end

table_data2=nan(length(behavior_list),4+sum(cellfun(@(x) size(x,2),table_data))-4*length(experiment_value));
for b=1:length(behavior_list)
  table_data2(b,1)=b;
  tmp1=cellfun(@(x) x(b,1:2),table_data,'uniformoutput',false);
  tmp2=cellfun(@(x) x(b,3:4),table_data,'uniformoutput',false);
  table_data2(b,2)=sum(sum([tmp1{:}])) / sum(sum([tmp2{:}])) .* 100;
  table_data2(b,3:4)=sum(cat(1,tmp1{:}),1) ./ sum(cat(1,tmp2{:}),1) .* 100;
  cellfun(@(x) x(b,5:end),table_data,'uniformoutput',false);
  table_data2(b,5:end)=cat(2,ans{:});
end

ret_val=cell(size(table_data2));
ret_val(:,1)=behavior_list(table_data2(:,1));
ret_val(:,2:end)=num2cell(table_data2(:,2:end));


% --- Executes on button press in BehaviorStats.
function BehaviorStats_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...');  drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_list=get(handles.BehaviorList,'String');
behavior_logic=get(handles.BehaviorLogic,'Value');
behavior_value2=get(handles.BehaviorList2,'Value');
behavior_list2=get(handles.BehaviorList2,'String');
feature_list=get(handles.FeatureList,'String');

if(length(experiment_value)>0)
  table_data=calculate_behaviorstats(experiment_value,experiment_list,...
      behavior_list,behavior_logic,behavior_value2,behavior_list2,feature_list);
end
if(length(experiment_value2)>0)
  table_data2=calculate_behaviorstats(experiment_value2,experiment_list2,...
      behavior_list,behavior_logic,behavior_value2,behavior_list2,feature_list);
end

if((length(experiment_value)>0) && (length(experiment_value2)>0))
  set(handles.Table,'RowStriping','on','BackgroundColor',[1 0.8 0.8; 0.8 0.8 1]);
  tmp(1:2:2*size(table_data2,1),1:size(table_data, 2))=table_data;
  tmp(2:2:2*size(table_data2,1),1:size(table_data2,2))=table_data2;
  ii=max(size(table_data,2),size(table_data2,2));
elseif(length(experiment_value)>0)
  set(handles.Table,'RowStriping','off','BackgroundColor',[1 0.8 0.8]);
  tmp=table_data;
  ii=size(table_data,2);
else
  set(handles.Table,'RowStriping','off','BackgroundColor',[0.8 0.8 1]);
  tmp=table_data2;
  ii=size(table_data2,2);
end

set(handles.Table,'Data',tmp);
tmp={'Behavior' 'Total %' 'Male %' 'Female %'};
for i=1:(ii-4)
  tmp{i+4}=['Indi #' num2str(i) ' %'];
end
set(handles.Table,'ColumnName',tmp);
set(handles.Table,'ColumnWidth',{150 75 75});

%handles.table_data2=table_data2;
handles.table='behavior_stats';
guidata(hObject,handles);

set(handles.Status,'string','Ready.');


% --- 
function Interesting_CellSelectionCallback(hObject,eventdata)

if(size(eventdata.Indices,1)==0)  return;  end

handles=guidata(hObject);

if(strcmp(handles.table,'behavior_stats'))
  if(eventdata.Indices(end,2)==1)  return;  end
  set(handles.Status,'string','Thinking...');  drawnow;

  experiment_value=get(handles.ExperimentList,'Value');
  experiment_list=get(handles.ExperimentList,'String');
  experiment_value2=get(handles.ExperimentList2,'Value');
  experiment_list2=get(handles.ExperimentList2,'String');
  behavior_list=get(handles.BehaviorList,'String');
  feature_list=get(handles.FeatureList,'String');

  if((length(experiment_value)>0) && (length(experiment_value2)>0))
    b=ceil(eventdata.Indices(end,1)/2);
    if(mod(eventdata.Indices(end,1),2))
      experiment_value0=experiment_value;
      experiment_list0=experiment_list;
      color='r';
    else
      experiment_value0=experiment_value2;
      experiment_list0=experiment_list2;
      color='b';
    end
  elseif(length(experiment_value)>0)
    b=eventdata.Indices(end,1);
    experiment_value0=experiment_value;
    experiment_list0=experiment_list;
      color='r';
  else
    b=eventdata.Indices(end,1);
    experiment_value0=experiment_value2;
    experiment_list0=experiment_list2;
    color='b';
  end

  sex_data=cell(1,length(experiment_value0));
  for e=1:length(experiment_value0)
    sex_data{e}=load(fullfile(char(experiment_list0(experiment_value0(e))),'perframe',...
        char(feature_list(find(strcmp(feature_list,'sex.mat'))))));
  end

  cellfun(@(x) size(x.data{1},2),sex_data,'uniformoutput',false);
  behavior_cumulative=zeros(1,max([ans{:}]));
  for e=1:length(experiment_value0)
    behavior_data=load(fullfile(char(experiment_list0(experiment_value0(e))),...
        char(behavior_list(b))));
    if(eventdata.Indices(end,2)>4)
      ii=eventdata.Indices(end,2);
    else
      ii=1:length(behavior_data.allScores.t0s);
      if(eventdata.Indices(end,2)>2)
        for i=1:length(sex_data{e}.data)
          sex_data{e}.data{i}=strcmp(sex_data{e}.data{i},'M');
        end
      end
    end
    for i=ii   % individual
      partition_idx=false(1,length(sex_data{e}.data{i}));
      for j=1:length(behavior_data.allScores.t0s{i})  % bout
        partition_idx((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
            -behavior_data.allScores.tStart(i)+1)=true;
      end
      switch(eventdata.Indices(end,2))
        case(3)
          partition_idx = partition_idx & sex_data{e}.data{i}(1:length(partition_idx));
        case(4)
          partition_idx = partition_idx & (~sex_data{e}.data{i}(1:length(partition_idx)));
      end
      idx=find(partition_idx);
      behavior_cumulative(idx)=behavior_cumulative(idx)+1;
    end
  end
  behavior_cumulative=conv(behavior_cumulative,ones(1,100)./100,'same');
  axes(handles.Axes);
  cla;  hold on;
  plot(behavior_cumulative,[color '-']);
  xlabel('time (frames)');
  ylabel(strrep(char(behavior_list(b)),'_','-'));
  axis tight;

  set(handles.Status,'string','Ready.');

else
  handles.behaviorvalue=handles.table_data(eventdata.Indices(end,1),1);
  set(handles.BehaviorList,'Value',handles.behaviorvalue);
  set(handles.BehaviorLogic,'Value',1);
  set(handles.BehaviorList2,'enable','off');
  handles.featurevalue=handles.table_data(eventdata.Indices(end,1),2);
  set(handles.FeatureList,'Value',handles.featurevalue);
  handles.sex=1;
  set(handles.SexList,'Value',handles.sex);

  if(strcmp(handles.table,'timeseries'))
    handles.timing=handles.table_data(eventdata.Indices(end,1),3);
    %set(handles.Timing,'Value',handles.timing);
    %handles.statistic=handles.table_data(eventdata.Indices(end,1),4);
    %set(handles.Statistic,'Value',handles.statistic);
    PlotTimeSeries_Callback(hObject, eventdata, handles);

  elseif(strcmp(handles.table,'histogram'))
    PlotHistogram_Callback(hObject, eventdata, handles);
  end
end

guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsTiming_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsTiming (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsStatistic_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsStatistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsRadius_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.windowradius=str2num(char(inputdlg({'Window radius:'},'',1,{num2str(10)})));
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsZoom_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsZoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsBars_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsBars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsBarsNone_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsBarsNone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.bar=1;
set(handles.MenuTimeSeriesPrefsBarsNone,'Checked','on');
set(handles.MenuTimeSeriesPrefsBarsDeviation,'Checked','off');
set(handles.MenuTimeSeriesPrefsBarsError,'Checked','off');
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsBarsDeviation_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsBarsDeviation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.bar=2;
set(handles.MenuTimeSeriesPrefsBarsNone,'Checked','off');
set(handles.MenuTimeSeriesPrefsBarsDeviation,'Checked','on');
set(handles.MenuTimeSeriesPrefsBarsError,'Checked','off');
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsBarsError_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsBarsError (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.bar=3;
set(handles.MenuTimeSeriesPrefsBarsNone,'Checked','off');
set(handles.MenuTimeSeriesPrefsBarsDeviation,'Checked','off');
set(handles.MenuTimeSeriesPrefsBarsError,'Checked','on');
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsZoomBlack_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsZoomBlack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.zoom=1;
set(handles.MenuTimeSeriesPrefsZoomBlack,'Checked','on');
set(handles.MenuTimeSeriesPrefsZoomRed,'Checked','off');
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsZoomRed_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsZoomRed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.zoom=2;
set(handles.MenuTimeSeriesPrefsZoomBlack,'Checked','off');
set(handles.MenuTimeSeriesPrefsZoomRed,'Checked','on');
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsStatisticMean_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsStatisticMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.statistic=1;
set(handles.MenuTimeSeriesPrefsStatisticMean,'Checked','on');
set(handles.MenuTimeSeriesPrefsStatisticMedian,'Checked','off');
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsStatisticMedian_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsStatisticMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.statistic=2;
set(handles.MenuTimeSeriesPrefsStatisticMean,'Checked','off');
set(handles.MenuTimeSeriesPrefsStatisticMedian,'Checked','on');
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsTimingOnset_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsTimingOnset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.timing=1;
set(handles.MenuTimeSeriesPrefsTimingOnset,'Checked','on');
set(handles.MenuTimeSeriesPrefsTimingOffset,'Checked','off');
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefsTimingOffset_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefsTimingOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.timing=2;
set(handles.MenuTimeSeriesPrefsTimingOnset,'Checked','off');
set(handles.MenuTimeSeriesPrefsTimingOffset,'Checked','on');
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuTimeSeriesPrefs_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesPrefs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Zoom.
function Zoom_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=zoom(gcf);
if(strcmp(get(h,'enable'),'off'))
  zoom on;
  pan off;
  set(handles.Zoom,'backgroundcolor',0.1*[1 1 1]);
  set(handles.Pan,'backgroundcolor',get(gcf,'color'));
else
  zoom off;
  pan off;
  set(handles.Zoom,'backgroundcolor',get(gcf,'color'));
end


% --- Executes on button press in Pan.
function Pan_Callback(hObject, eventdata, handles)
% hObject    handle to Pan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

h=pan(gcf);
if(strcmp(get(h,'enable'),'off'))
  pan on;
  zoom off;
  set(handles.Pan,'backgroundcolor',0.1*[1 1 1]);
  set(handles.Zoom,'backgroundcolor',get(gcf,'color'));
else
  pan off;
  zoom off;
  set(handles.Pan,'backgroundcolor',get(gcf,'color'));
end
