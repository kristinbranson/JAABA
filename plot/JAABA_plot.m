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

% Last Modified by GUIDE v2.5 02-Jul-2012 17:27:02

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
  handles.behaviorlist=handles_saved.behaviorlist;
  handles.behaviorvalue=handles_saved.behaviorvalue;
  handles.behaviorlogic=handles_saved.behaviorlogic;
  handles.behaviorlist2=handles_saved.behaviorlist2;
  handles.behaviorvalue2=handles_saved.behaviorvalue2;
  handles.featurelist=handles_saved.featurelist;
  handles.featurevalue=handles_saved.featurevalue;
  handles.timing=handles_saved.timing;
  handles.statistic=handles_saved.statistic;
  handles.windowradius=handles_saved.windowradius;
  handles.bar=handles_saved.bar;
  handles.zoom=handles_saved.zoom;
  handles.individuals=handles_saved.individuals;
  if(length(handles.experimentlist)>0)
    set(handles.ExperimentList,'String',handles.experimentlist);
  else
    set(handles.ExperimentList,'String',{'Listbox'});
  end
  set(handles.ExperimentList,'Value',handles.experimentvalue);
  set(handles.BehaviorList,'String',handles.behaviorlist);
  set(handles.BehaviorList,'Value',handles.behaviorvalue);
  set(handles.BehaviorLogic,'Value',handles.behaviorlogic);
  if(handles.behaviorlogic==1)
    set(handles.BehaviorList2,'enable','off');
  else
    set(handles.BehaviorList2,'enable','on');
  end
  set(handles.BehaviorList2,'String',handles.behaviorlist2);
  set(handles.BehaviorList2,'Value',handles.behaviorvalue2);
  set(handles.FeatureList,'String',handles.featurelist);
  set(handles.FeatureList,'Value',handles.featurevalue);
  set(handles.Timing,'Value',handles.timing);
  set(handles.Statistic,'Value',handles.statistic);
  set(handles.WindowRadius,'String',handles.windowradius);
  set(handles.Bar,'Value',handles.bar);
  set(handles.Zoom,'Value',handles.zoom);
  set(handles.Individuals,'Value',handles.individuals);
else
  handles.experimentlist={};
  handles.experimentvalue=1;
  handles.behaviorlist={''};
  handles.behaviorvalue=1;
  handles.behaviorlogic=1;
  handles.behaviorlist2={''};
  handles.behaviorvalue2=1;
  handles.featurelist={''};
  handles.featurevalue=1;
  handles.timing=1;
  handles.statistic=1;
  handles.windowradius=10;
  handles.bar=1;
  handles.zoom=1;
  handles.individuals=1;
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
handles.experimentvalue=1;
handles.behaviorlist={''};
handles.behaviorvalue=1;
handles.behaviorlogic=1;
handles.behaviorlist2={''};
handles.behaviorvalue2=1;
handles.featurelist={''};
handles.featurevalue=1;
handles.timing=1;
handles.statistic=1;
handles.windowradius=10;
handles.bar=1;
handles.zoom=1;
handles.individuals=1;

set(handles.ExperimentList,'String',{'Listbox'});
set(handles.BehaviorList,'String',handles.behaviorlist);
set(handles.BehaviorLogic,'Value',handles.behaviorlogic);
  set(handles.BehaviorList2,'enable','off');
set(handles.BehaviorList2,'String',handles.behaviorlist2);
set(handles.FeatureList,'String',handles.featurelist);
set(handles.Timing,'Value',handles.timing);
set(handles.Statistic,'Value',handles.statistic);
set(handles.WindowRadius,'String',handles.windowradius);
set(handles.Bar,'Value',handles.bar);
set(handles.Zoom,'Value',handles.zoom);
set(handles.Individuals,'Value',handles.individuals);

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


% --- Executes on button press in ExperimentDelete.
function ExperimentDelete_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)


% --- Executes on button press in ExperimentAdd.
function ExperimentAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...');

handles.experimentlist{end+1}=uigetdir([],'Select Experiment Directory');
set(handles.ExperimentList,'String',handles.experimentlist);

tmp=dir(fullfile(handles.experimentlist{end},'scores*.mat'));
[handles.behaviorlist{1:length(tmp)}]=deal(tmp.name);
set(handles.BehaviorList,'String',handles.behaviorlist);
[handles.behaviorlist2{1:length(tmp)}]=deal(tmp.name);
set(handles.BehaviorList2,'String',handles.behaviorlist2);

tmp=dir(fullfile(handles.experimentlist{end},'perframe','*.mat'));
[handles.featurelist{1:length(tmp)}]=deal(tmp.name);
set(handles.FeatureList,'String',handles.featurelist);

guidata(hObject,handles);

set(handles.Status,'string','Ready.');


function WindowRadius_Callback(hObject, eventdata, handles)
% hObject    handle to WindowRadiusText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowRadiusText as text
%        str2double(get(hObject,'String')) returns contents of WindowRadiusText as a double

handles.windowradius=get(handles.WindowRadius,'String');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function WindowRadius_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Timing.
function Timing_Callback(hObject, eventdata, handles)
% hObject    handle to Timing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Timing contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Timing

handles.timing=get(handles.Timing,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Timing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Timing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Statistic.
function Statistic_Callback(hObject, eventdata, handles)
% hObject    handle to Statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Statistic contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Statistic

handles.statistic=get(handles.Statistic,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Statistic_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Statistic (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Bar.
function Bar_Callback(hObject, eventdata, handles)
% hObject    handle to Bar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Bar contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Bar

handles.bar=get(handles.Bar,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Bar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Zoom.
function Zoom_Callback(hObject, eventdata, handles)
% hObject    handle to Zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Zoom contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Zoom

handles.zoom=get(handles.Zoom,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function Zoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Zoom (see GCBO)
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

during=[];  not_during=[];
partition_idx=cell(1,length(feature_data.data));
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
  during=[during feature_data.data{i}(partition_idx{i} & sex_data.data{i}(1:length(partition_idx{i})))];
  not_during=[not_during feature_data.data{i}((~partition_idx{i}) & sex_data.data{i}(1:length(partition_idx{i})))];
end


% --- Executes on button press in PlotHistogram.
function PlotHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to PlotHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...');  drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
behavior_value=get(handles.BehaviorList,'Value');
behavior_list=get(handles.BehaviorList,'String');
behavior_logic=get(handles.BehaviorLogic,'Value');
behavior_value2=get(handles.BehaviorList2,'Value');
behavior_list2=get(handles.BehaviorList2,'String');
feature_value=get(handles.FeatureList,'Value');
feature_list=get(handles.FeatureList,'String');
individuals=get(handles.Individuals,'Value');

behavior_data=load(fullfile(char(experiment_list(experiment_value)),...
      char(behavior_list(behavior_value))));
if(behavior_logic>1)
  behavior_data2=load(fullfile(char(experiment_list(experiment_value)),...
      char(behavior_list2(behavior_value2))));
else
  behavior_data2=[];
end
feature_data=load(fullfile(char(experiment_list(experiment_value)),'perframe',...
      char(feature_list(feature_value))));

sex_data=load(fullfile(char(experiment_list(experiment_value)),'perframe',...
    char(feature_list(find(strcmp(feature_list,'sex.mat'))))));
for i=1:length(sex_data.data)
  switch(individuals)
    case(1)
      sex_data.data{i}=ones(1,length(sex_data.data{i}));
    case(2)
      sex_data.data{i}=strcmp(sex_data.data{i},'M');
    case(3)
      sex_data.data{i}=strcmp(sex_data.data{i},'F');
  end
end

[during not_during]=calculate_histogram(behavior_data,behavior_logic,behavior_data2,feature_data,sex_data);

axes(handles.Axes);
cla;  hold on;
tmp=linspace(min([during not_during]),max([during not_during]));
hist_during=hist(during,tmp);
hist_not_during=hist(not_during,tmp);
plot(tmp,hist_not_during./max(hist_not_during),'r.-');
plot(tmp,hist_during./max(hist_during),'b.-');
text(tmp(5),0.95,['during = ' num2str(100*length(during)/(length(during)+length(not_during)),3) '%'],'color','b');
xlabel(get_label(feature_list(feature_value),feature_data.units));
ylabel('normalized');
axis tight

set(handles.Status,'string','Ready.');


% --- Executes on button press in InterestingHistograms.
function InterestingHistograms_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingHistograms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...');  drawnow;

ncols=3;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
behavior_list=get(handles.BehaviorList,'String');
feature_list=get(handles.FeatureList,'String');

feature_data=load(fullfile(char(experiment_list(experiment_value)),'perframe',char(feature_list(1))));
for i=1:length(feature_data.data)
  sex_data.data{i}=ones(1,length(feature_data.data{i}));
end

table_data=zeros(length(behavior_list),length(feature_list),ncols);
parfor b=1:length(behavior_list)
  k=1;
  tmp=zeros(length(feature_list),ncols);
  behavior_data=load(fullfile(char(experiment_list(experiment_value)),char(behavior_list(b))));
  for f=1:length(feature_list)
    feature_data=load(fullfile(char(experiment_list(experiment_value)),'perframe',char(feature_list(f))));
    [during not_during]=calculate_histogram(behavior_data,1,[],feature_data,sex_data);
    tmp(k,:)=[b f (mean(during)-mean(not_during))/sqrt((std(during)^2+std(not_during)^2)/2)];
    k=k+1;
  end
  table_data(b,:,:)=tmp;
  disp([num2str(b) ' of ' num2str(length(behavior_list))]);
end
table_data=reshape(table_data,prod(size(table_data))/ncols,ncols);
table_data=sortrows(table_data,-ncols);

tmp=cell(size(table_data));
tmp(:,1)=behavior_list(table_data(:,1));
tmp(:,2)=feature_list(table_data(:,2));
tmp(:,3)=num2cell(table_data(:,3));
set(handles.Table,'Data',tmp);
set(handles.Table,'ColumnName',{'Behavior' 'Feature' 'd'''});
set(handles.Table,'ColumnWidth',{150 150 50});

handles.table_data=table_data;
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
triggered_data=zeros(5+length([behavior_data.allScores.t0s{:}]),1+2*windowradius);
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
switch(statistic)
  case(1)
    triggered_data(end-2,:)=nanmean(triggered_data(1:k-1,:),1);
    tmp=nanstd(triggered_data(1:k-1,:),1)./sqrt(size(triggered_data,1));
    triggered_data(end-3,:)=triggered_data(end-2,:)+tmp;
    triggered_data(end-1,:)=triggered_data(end-2,:)-tmp;
    tmp=nanstd(triggered_data(1:k-1,:),1);
    triggered_data(end-4,:)=triggered_data(end-2,:)+tmp;
    triggered_data(end-0,:)=triggered_data(end-2,:)-tmp;
  case(2)
    triggered_data(end-4,:)=prctile(triggered_data(1:k-1,:),95,1);
    triggered_data(end-3,:)=prctile(triggered_data(1:k-1,:),75,1);
    triggered_data(end-2,:)=nanmedian(triggered_data(1:k-1,:),1);
    triggered_data(end-1,:)=prctile(triggered_data(1:k-1,:),25,1);
    triggered_data(end-0,:)=prctile(triggered_data(1:k-1,:), 5,1);
end


% --- Executes on button press in PlotTimeSeries.
function PlotTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...');  drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
behavior_value=get(handles.BehaviorList,'Value');
behavior_list=get(handles.BehaviorList,'String');
behavior_logic=get(handles.BehaviorLogic,'Value');
behavior_value2=get(handles.BehaviorList2,'Value');
behavior_list2=get(handles.BehaviorList2,'String');
feature_value=get(handles.FeatureList,'Value');
feature_list=get(handles.FeatureList,'String');
timing=get(handles.Timing,'Value');
statistic=get(handles.Statistic,'Value');
windowradius=str2num(get(handles.WindowRadius,'string'));

behavior_data=load(fullfile(char(experiment_list(experiment_value)),...
    char(behavior_list(behavior_value))));
if(behavior_logic>1)
  behavior_data2=load(fullfile(char(experiment_list(experiment_value)),...
      char(behavior_list2(behavior_value2))));
else
  behavior_data2=[];
end
feature_data=load(fullfile(char(experiment_list(experiment_value)),'perframe',...
    char(feature_list(feature_value))));

%figure;  hold on;
%plot(behavior_data.allScores.scores{1});
%axis tight;  v=axis;
%for i=1:length(behavior_data.allScores.t0s{1})
%  plot(behavior_data.allScores.t0s{1}([i i]),[v(3) v(4)],'r-');
%  plot(behavior_data.allScores.t1s{1}([i i]),[v(3) v(4)],'r:');
%end

triggered_data=calculate_timeseries(behavior_data,behavior_logic,behavior_data2,feature_data,timing,statistic,windowradius);

axes(handles.Axes);
cla;  hold on;
plot(-windowradius:windowradius,triggered_data(1:end-5,:)','k-');
plot(-windowradius:windowradius,triggered_data(end-2,:)','r-','linewidth',3);
range_red=[min(triggered_data(end-2,:)) max(triggered_data(end-2,:))];
switch(handles.bar)
  case(2)
    plot(-windowradius:windowradius,triggered_data(end-4,:)','r-','linewidth',1);
    plot(-windowradius:windowradius,triggered_data(end-0,:)','r-','linewidth',1);
    range_red=[range_red; min(triggered_data(end-4,:)) max(triggered_data(end-4,:))];
    range_red=[range_red; min(triggered_data(end-0,:)) max(triggered_data(end-0,:))];
  case(3)
    plot(-windowradius:windowradius,triggered_data(end-3,:)','r-','linewidth',1);
    plot(-windowradius:windowradius,triggered_data(end-1,:)','r-','linewidth',1);
    range_red=[range_red; min(triggered_data(end-3,:)) max(triggered_data(end-3,:))];
    range_red=[range_red; min(triggered_data(end-1,:)) max(triggered_data(end-1,:))];
end
xlabel('time (frames)');
ylabel(get_label(feature_list(feature_value),feature_data.units));
axis tight
if(handles.zoom==2)
  v=axis;  axis([v(1) v(2) min(range_red(:,1)) max(range_red(:,2))]);
end

set(handles.Status,'string','Ready.');


% --- Executes on button press in InterestingTimeSeries.
function InterestingTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...');  drawnow;

ncols=4;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
behavior_list=get(handles.BehaviorList,'String');
feature_list=get(handles.FeatureList,'String');
windowradius=str2num(get(handles.WindowRadius,'string'));
statistic=get(handles.Statistic,'Value');

table_data=zeros(length(behavior_list),length(feature_list)*2,ncols);
parfor b=1:length(behavior_list)
  k=1;
  parfor_tmp=zeros(length(feature_list)*2,ncols);
  behavior_data=load(fullfile(char(experiment_list(experiment_value)),char(behavior_list(b))));
  for f=1:length(feature_list)
    feature_data=load(fullfile(char(experiment_list(experiment_value)),'perframe',char(feature_list(f))));
    for t=1:2  % timing
      triggered_data=calculate_timeseries(behavior_data,1,[],feature_data,t,statistic,windowradius);
      foo=nanmean(triggered_data(1:end-5,1:windowradius),2);
      foo=triggered_data(1:end-5,:)-repmat(foo,1,size(triggered_data,2));
      foo=sqrt(nanmean(triggered_data(1:end-5,(windowradius+1):end).^2,2)) ./ ...
          sqrt(nanmean(triggered_data(1:end-5,1:windowradius).^2,2));
      parfor_tmp(k,:)=[b f t nanmedian(foo)];
%      parfor_tmp(k,:)=[b f t ...
%          nanmean(triggered_data(end-2,(windowradius+1):end))/nanmean(triggered_data(end-2,1:windowradius)) ...
%          nanstd(triggered_data(end-2,(windowradius+1):end))/ nanstd(triggered_data(end-2,1:windowradius))];
      k=k+1;
    end
  end
  table_data(b,:,:)=parfor_tmp;
  disp([num2str(b) ' of ' num2str(length(behavior_list))]);
end
table_data=reshape(table_data,prod(size(table_data))/ncols,ncols);
table_data=sortrows(table_data,-ncols);
%table_data=sortrows(table_data,-(ncols-1));

tmp=cell(size(table_data));
tmp(:,1)=behavior_list(table_data(:,1));
tmp(:,2)=feature_list(table_data(:,2));
tmp(:,3)=repmat({'Onset'},size(table_data,1),1);
  find(table_data(:,3)==2);
  tmp(ans,3)=repmat({'Offset'},size(ans),1);
%tmp(:,4)=repmat({'Mean'},size(table_data,1),1);
%  find(table_data(:,4)==2);
%  tmp(ans,4)=repmat({'Median'},size(ans),1);
tmp(:,4)=num2cell(table_data(:,4));
%tmp(:,6)=num2cell(table_data(:,6));

set(handles.Table,'Data',tmp);
set(handles.Table,'ColumnName',{'Behavior' 'Feature' 'Timing' 'dRMS'});
%set(handles.Table,'ColumnName',{'Behavior' 'Feature' 'Timing' 'Statistic' 'Mean' 'Std'});
set(handles.Table,'ColumnWidth',{150 150 50 50 50});

handles.table_data=table_data;
handles.table='timeseries';
guidata(hObject,handles);

set(handles.Status,'string','Ready.');


% --- Executes on selection change in Individuals.
function Individuals_Callback(hObject, eventdata, handles)
% hObject    handle to Individuals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Individuals contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Individuals


% --- Executes during object creation, after setting all properties.
function Individuals_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Individuals (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- 
function [male_pct female_pct individual_pct]=calculate_stats(behavior_data,sex_data)

male_pct=0;  female_pct=0;  individual_pct=[];
male_total=0;  female_total=0;
for i=1:length(behavior_data.allScores.t0s)
  partition_idx=false(1,length(sex_data.data{i}));
  for j=1:length(behavior_data.allScores.t0s{i})
    partition_idx((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
        -behavior_data.allScores.tStart(i)+1)=true;
  end
  male_pct=male_pct+sum(partition_idx & sex_data.data{i}(1:length(partition_idx)));
  female_pct=female_pct+sum(partition_idx & (~sex_data.data{i}(1:length(partition_idx))));
  male_total=male_total+sum(sex_data.data{i}(1:length(partition_idx)));
  female_total=female_total+sum(~sex_data.data{i}(1:length(partition_idx)));
  individual_pct(i)=100*sum(partition_idx)./length(partition_idx);
end
male_pct=100*male_pct/male_total;
female_pct=100*female_pct/female_total;


% --- Executes on button press in BehaviorStats.
function BehaviorStats_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...');  drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
behavior_list=get(handles.BehaviorList,'String');
feature_list=get(handles.FeatureList,'String');

sex_data=load(fullfile(char(experiment_list(experiment_value)),'perframe',...
    char(feature_list(find(strcmp(feature_list,'sex.mat'))))));
for i=1:length(sex_data.data)
  sex_data.data{i}=strcmp(sex_data.data{i},'M');
end

behavior_data=load(fullfile(char(experiment_list(experiment_value)),char(behavior_list(1))));
table_data=nan(length(behavior_list),3+length(behavior_data.allScores.t0s));
for b=1:length(behavior_list)
  behavior_data=load(fullfile(char(experiment_list(experiment_value)),char(behavior_list(b))));

  [male_pct female_pct individual_pct]=calculate_stats(behavior_data,sex_data);
  if((3+length(individual_pct))>size(table_data,2))
    table_data(:,(end+1):(3+length(individual_pct)))=...
        nan(size(table_data,1):((3+length(individual_pct))-size(table_data,2)));
  end
  table_data(b,:)=[b male_pct female_pct individual_pct];
  %disp([num2str(b) ' of ' num2str(length(behavior_list))]);
end

tmp=cell(size(table_data));
tmp(:,1)=behavior_list(table_data(:,1));
tmp(:,2:end)=num2cell(table_data(:,2:end));
set(handles.Table,'Data',tmp);
tmp={'Behavior' 'Male %' 'Female %'};
for i=1:(size(table_data,2)-3)
  tmp{i+3}=['Indi #' num2str(i)];
end
set(handles.Table,'ColumnName',tmp);
set(handles.Table,'ColumnWidth',{150 75 75});

handles.table_data=table_data;
handles.table='behavior_stats';
guidata(hObject,handles);

set(handles.Status,'string','Ready.');


% --- 
function Interesting_CellSelectionCallback(hObject,eventdata)

if(size(eventdata.Indices,1)==0)  return;  end

handles=guidata(hObject);

if(strcmp(handles.table,'behavior_stats'))
  set(handles.Status,'string','Thinking...');  drawnow;

  experiment_value=get(handles.ExperimentList,'Value');
  experiment_list=get(handles.ExperimentList,'String');
  behavior_list=get(handles.BehaviorList,'String');
  feature_list=get(handles.FeatureList,'String');

  behavior_data=load(fullfile(char(experiment_list(experiment_value)),...
      char(behavior_list(eventdata.Indices(end,1)))));
  sex_data=load(fullfile(char(experiment_list(experiment_value)),'perframe',...
      char(feature_list(find(strcmp(feature_list,'sex.mat'))))));
  if(eventdata.Indices(end,2)>3)
    ii=eventdata.Indices(end,2);
  else
    ii=1:length(behavior_data.allScores.t0s);
    if(eventdata.Indices(end,2)>1)
      for i=1:length(sex_data.data)
        sex_data.data{i}=strcmp(sex_data.data{i},'M');
      end
    end
  end
  behavior_cumulative=zeros(1,length(sex_data.data{1}));
  for i=ii
    partition_idx=false(1,length(sex_data.data{i}));
    for j=1:length(behavior_data.allScores.t0s{i})
      partition_idx((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
          -behavior_data.allScores.tStart(i)+1)=true;
    end
    switch(eventdata.Indices(end,2))
      case(2)
        partition_idx = partition_idx & sex_data.data{i}(1:length(partition_idx));
      case(3)
        partition_idx = partition_idx & (~sex_data.data{i}(1:length(partition_idx)));
    end
    behavior_cumulative(partition_idx)=behavior_cumulative(partition_idx)+1;
  end
  behavior_cumulative=conv(behavior_cumulative,ones(1,100)./100,'same');
  axes(handles.Axes);
  cla;  hold on;
  plot(behavior_cumulative);
  xlabel('time (frames)');
  ylabel(strrep(char(behavior_list(eventdata.Indices(end,1))),'_','-'));
  axis tight;

  set(handles.Status,'string','Ready.');

else
  handles.behaviorvalue=handles.table_data(eventdata.Indices(end,1),1);
  set(handles.BehaviorList,'Value',handles.behaviorvalue);
  set(handles.BehaviorLogic,'Value',1);
  set(handles.BehaviorList2,'enable','off');
  handles.featurevalue=handles.table_data(eventdata.Indices(end,1),2);
  set(handles.FeatureList,'Value',handles.featurevalue);
  handles.individuals=1;
  set(handles.Individuals,'Value',handles.individuals);

  if(strcmp(handles.table,'timeseries'))
    handles.timing=handles.table_data(eventdata.Indices(end,1),3);
    set(handles.Timing,'Value',handles.timing);
    %handles.statistic=handles.table_data(eventdata.Indices(end,1),4);
    %set(handles.Statistic,'Value',handles.statistic);
    PlotTimeSeries_Callback(hObject, eventdata, handles);

  elseif(strcmp(handles.table,'histogram'))
    PlotHistogram_Callback(hObject, eventdata, handles);
  end
end

guidata(hObject,handles);
