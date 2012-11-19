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

% Last Modified by GUIDE v2.5 08-Nov-2012 17:28:34

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


% ---
function handles=initialize(handles)

handles.grouplist={};
handles.groupvalue=1;
handles.experimentlist={{}};
handles.experimentvalue={1};
handles.behaviors={};
handles.behaviorlist={};
handles.behaviorvalue=1;
handles.behaviorlogic=1;
handles.behaviorvalue2=1;
handles.features={};
handles.featurelist={};
handles.featurevalue=1;
handles.individuals=[];
handles.individuallist={'All' 'Male' 'Female'};
handles.individualvalue=1;
%handles.individual=1;
handles.sexdata={};
handles.behaviorbarchart_perwhat=1;
handles.behaviortimeseries_style=1;
handles.featurehistogram_perwhat=1;
handles.featurehistogram_style=1;
handles.featurehistogram_logbinsize=0;
handles.featurehistogram_notduring=0;
handles.featurehistogram_nbins=100;
handles.featuretimeseries_timing=1;
handles.featuretimeseries_style=1;
handles.featuretimeseries_subtractmean=0;
handles.timeseries_windowradius=10;
handles.interestingfeaturehistograms_omitnan=1;
handles.interestingfeaturehistograms_omitinf=1;
%handles.timeseries_tight=0;
handles.prefs_centraltendency=1;
handles.prefs_dispersion=1;
handles.prefs_convolutionwidth=1000;
%handles.logy=0;
%handles.stats=0;
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
handles.colors={'r' 'g' 'b' 'c' 'm' 'y' 'k';
                'red' 'green' 'blue' 'cyan' 'magenta' 'yellow' 'black'};


% ---
function handles=load_configuration_file(filename,hObject,eventdata,handles)

handles_saved=load(filename);
handles_saved=handles_saved.handles;
handles.grouplist=handles_saved.grouplist;
handles.groupvalue=handles_saved.groupvalue;
handles.experimentlist=handles_saved.experimentlist;
handles.experimentvalue=handles_saved.experimentvalue;
handles.behaviors=handles_saved.behaviors;
handles.behaviorlist=handles_saved.behaviorlist;
handles.behaviorvalue=handles_saved.behaviorvalue;
handles.behaviorlogic=handles_saved.behaviorlogic;
handles.behaviorvalue2=handles_saved.behaviorvalue2;
handles.features=handles_saved.features;
handles.featurelist=handles_saved.featurelist;
handles.featurevalue=handles_saved.featurevalue;
handles.individuals=handles_saved.individuals;
handles.individuallist=handles_saved.individuallist;
handles.individualvalue=handles_saved.individualvalue;
handles.sexdata=handles_saved.sexdata;
handles.behaviorbarchart_perwhat=handles_saved.behaviorbarchart_perwhat;
handles.behaviortimeseries_style=handles_saved.behaviortimeseries_style;
handles.featurehistogram_perwhat=handles_saved.featurehistogram_perwhat;
handles.featurehistogram_style=handles_saved.featurehistogram_style;
handles.featurehistogram_logbinsize=handles_saved.featurehistogram_logbinsize;
handles.featurehistogram_notduring=handles_saved.featurehistogram_notduring;
handles.featurehistogram_nbins=handles_saved.featurehistogram_nbins;
handles.featuretimeseries_timing=handles_saved.featuretimeseries_timing;
handles.featuretimeseries_style=handles_saved.featuretimeseries_style;
handles.featuretimeseries_subtractmean=handles_saved.featuretimeseries_subtractmean;
handles.timeseries_windowradius=handles_saved.timeseries_windowradius;
handles.interestingfeaturehistograms_omitnan=handles_saved.interestingfeaturehistograms_omitnan;
handles.interestingfeaturehistograms_omitinf=handles_saved.interestingfeaturehistograms_omitinf;
%handles.timeseries_tight=handles_saved.timeseries_tight;
handles.prefs_centraltendency=handles_saved.prefs_centraltendency;
handles.prefs_dispersion=handles_saved.prefs_dispersion;
handles.prefs_convolutionwidth=handles_saved.prefs_convolutionwidth;
%handles.logy=handles_saved.logy;
%handles.stats=handles_saved.stats;
handles.interestingfeaturehistograms_cache=handles_saved.interestingfeaturehistograms_cache;
handles.interestingfeaturetimeseries_cache=handles_saved.interestingfeaturetimeseries_cache;
handles.colors=handles_saved.colors;
handles.table=[];

%set(handles.MenuTimeSeriesTight,'Checked','off');
%if(handles.logy)   set(handles.LogX, 'backgroundcolor',0.4*[1 1 1]);  end
%if(handles.stats)  set(handles.Stats,'backgroundcolor',0.4*[1 1 1]);  end


% ---
function update_figure(handles)

if(isempty(handles.grouplist) || length(handles.experimentlist{handles.groupvalue})==0)
  set(handles.ExperimentList,'enable','off');
else
  set(handles.ExperimentList,'enable','on');
end
if(isempty(handles.grouplist))
  set(handles.GroupList,'enable','off');
  set(handles.ExperimentAdd,'enable','off');
  set(handles.ExperimentDelete,'enable','off');
  set(handles.ExperimentMove,'enable','off');
else
  set(handles.GroupList,'enable','on');
  set(handles.ExperimentAdd,'enable','on');
  set(handles.ExperimentDelete,'enable','on');
  set(handles.ExperimentMove,'enable','on');
end
if(sum(cellfun(@length,handles.experimentlist))==0)
  set(handles.BehaviorList,'enable','off');
  set(handles.BehaviorLogic,'enable','off');
  set(handles.BehaviorList2,'enable','off');
  set(handles.FeatureList,'enable','off');
  set(handles.IndividualList,'enable','off');
else
  set(handles.BehaviorList,'enable','on');
  set(handles.BehaviorLogic,'enable','on');
  set(handles.BehaviorList2,'enable','on');
  set(handles.FeatureList,'enable','on');
  set(handles.IndividualList,'enable','on');
end
if(handles.behaviorlogic==1)
  set(handles.BehaviorList2,'enable','off');
else
  set(handles.BehaviorList2,'enable','on');
end

if(isempty(handles.grouplist))
  set(handles.GroupList,'String',{''},'Value',1);
else
  tmp=length(handles.grouplist);
  cellstr(strcat(repmat('<html><font color="',tmp,1),...
      {handles.colors{2,1+mod(0:length(handles.grouplist)-1,length(handles.colors))}}',...
      repmat('">',tmp,1),handles.grouplist',repmat('</font></html>',tmp,1)));
  set(handles.GroupList,'String',ans,'Value',handles.groupvalue);
end
if(isempty(handles.experimentlist))
  set(handles.ExperimentList,'String',{''},'Value',1);
else
  set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue});
  set(handles.ExperimentList,'Value',handles.experimentvalue{handles.groupvalue});
end
if(isempty(handles.behaviorlist))
  set(handles.BehaviorList,'String',{''},'Value',1);
  set(handles.BehaviorList2,'String',{''},'Value',1);
else
  set(handles.BehaviorList,'String',handles.behaviorlist,'Value',handles.behaviorvalue);
  set(handles.BehaviorList2,'String',handles.behaviorlist,'Value',handles.behaviorvalue2);
end
  set(handles.BehaviorLogic,'Value',handles.behaviorlogic);
if(isempty(handles.featurelist))
  set(handles.FeatureList,'String',{''},'Value',1);
else
  set(handles.FeatureList,'String',handles.featurelist,'Value',handles.featurevalue);
end
if(isempty(handles.individuallist))
  set(handles.IndividualList,'String',{''},'Value',1);
else
  set(handles.IndividualList,'String',handles.individuallist,'Value',handles.individualvalue);
end
set(handles.Table,'Data',[]);

menu_behaviorbarchart_perwhat_set(handles.behaviorbarchart_perwhat);
menu_behaviortimeseries_style_set(handles.behaviortimeseries_style);
menu_featurehistogram_perwhat_set(handles.featurehistogram_perwhat);
menu_featurehistogram_style_set(handles.featurehistogram_style);
menu_featurehistogram_logbinsize_set(handles.featurehistogram_logbinsize);
menu_featurehistogram_notduring_set(handles.featurehistogram_notduring);
menu_featuretimeseries_timing_set(handles.featuretimeseries_timing);
menu_featuretimeseries_style_set(handles.featuretimeseries_style);
menu_featuretimeseries_subtractmean_set(handles.featuretimeseries_subtractmean);
menu_interestingfeaturehistograms_omitnan_set(handles.interestingfeaturehistograms_omitnan);
menu_interestingfeaturehistograms_omitinf_set(handles.interestingfeaturehistograms_omitinf);
menu_prefscentraltendency_set(handles.prefs_centraltendency);
menu_prefsdispersion_set(handles.prefs_dispersion);


% --- Executes just before JAABA_plot is made visible.
function JAABA_plot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

jlabelpath = fileparts(mfilename('fullpath'));
baseDir = fileparts(jlabelpath);
addpath(fullfile(baseDir,'misc'));

if(exist('matlabpool')==2 && matlabpool('size')==0)
  matlabpool open
end

if(exist('most_recent_config.mat')==2)
  handles=load_configuration_file('most_recent_config.mat',hObject,eventdata,handles);
else
  handles=initialize(handles);
end
update_figure(handles);


% Choose default command line output for JAABA_plot
handles.output = hObject;

set(hObject,'CloseRequestFcn',@figure_CloseRequestFcn);
set(handles.Table,'CellSelectionCallback',@CellSelectionCallback);
set(handles.ExperimentList,'Callback',@ListboxCallback);
%set(handles.ExperimentList2,'Callback',@ListboxCallback);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JAABA_plot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function ListboxCallback(hObject,eventdata)

handles=guidata(hObject);
handles.experimentvalue{handles.groupvalue}=get(handles.ExperimentList,'Value');
%handles.experimentvalue2=get(handles.ExperimentList2,'Value');
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


% ---
function figure_CloseRequestFcn(hObject, eventdata)

handles=guidata(hObject);
save('most_recent_config.mat','handles');
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


% --- Executes on selection change in IndividualList.
function IndividualList_Callback(hObject, eventdata, handles)
% hObject    handle to IndividualList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns IndividualList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from IndividualList

handles.individualvalue=get(handles.IndividualList,'Value');
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function IndividualList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to IndividualList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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

handles.behaviorvalue2=get(handles.BehaviorList2,'Value');
guidata(hObject,handles);


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
function handles=fillin_individuallist(handles)

cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];

tmp=cell(1,3+sum(handles.individuals));
tmp(1:3)={'All' 'Male' 'Female'};
k=4;
for e=1:length(handles.individuals)
  g=find(cumsum_num_exp_per_group<e,1,'last');
  for i=1:handles.individuals(e)
    %c='R';  n=e;  if(e>length(handles.experimentlist))  c='B';  n=n-length(handles.experimentlist);  end
    tmp{k}=['Grp ' handles.grouplist{g} ', Exp #' num2str(e-cumsum_num_exp_per_group(g)) ', Indi #' num2str(i)];
    k=k+1;
  end
end
handles.individuallist=tmp;
set(handles.IndividualList,'String',handles.individuallist);


% --- Executes on button press in ExperimentAdd.
function ExperimentAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

persistent directory
if(isempty(directory))  directory=pwd;  end

newrecordings=uipickfiles('prompt','Select experiment directory','filterspec',directory);
if(~iscell(newrecordings) || (length(newrecordings)==0))  return;  end

set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

[directory,~,~]=fileparts(newrecordings{1});
handles.experimentlist{handles.groupvalue}={handles.experimentlist{handles.groupvalue}{:} newrecordings{:}};
handles.experimentvalue{handles.groupvalue}=1:length(handles.experimentlist{handles.groupvalue});
set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue});
set(handles.ExperimentList,'Value',handles.experimentvalue{handles.groupvalue});
set(handles.ExperimentList,'enable','on');

set(handles.BehaviorList,'enable','on');
set(handles.BehaviorLogic,'enable','on');
if(handles.behaviorlogic>1)
  set(handles.BehaviorList2,'enable','on');
end
set(handles.FeatureList,'enable','on');
set(handles.IndividualList,'enable','on');

handlesbehaviors=cell(1,length(newrecordings));
handlesfeatures=cell(1,length(newrecordings));
handlesindividuals=[];
handlessexdata=cell(1,length(newrecordings));
parfor n=1:length(newrecordings)
  tmp=dir(fullfile(newrecordings{n},'scores*.mat'));
  [handlesbehaviors{n}{1:length(tmp)}]=deal(tmp.name);
  handlesbehaviors{n}=cellfun(@(x) x(1:(end-4)),handlesbehaviors{n},'uniformoutput',false);

  tmp=dir(fullfile(newrecordings{n},'perframe','*.mat'));
  [handlesfeatures{n}{1:length(tmp)}]=deal(tmp.name);
  handlesfeatures{n}=cellfun(@(x) x(1:(end-4)),handlesfeatures{n},'uniformoutput',false);

  behavior_data=load(fullfile(newrecordings{n},[handlesbehaviors{n}{1} '.mat']));
  handlesindividuals(n)=length(behavior_data.allScores.t0s);

  sexdata=load(fullfile(newrecordings{n},'perframe','sex.mat'));
  sexdata.data=cellfun(@(x) strcmp(x,'M'),sexdata.data,'uniformoutput',false);
  handlessexdata(n)={sexdata.data};
end
handles.behaviors={handles.behaviors{:} handlesbehaviors{:}};
handles.features={handles.features{:} handlesfeatures{:}};
handles.individuals=[handles.individuals handlesindividuals];
handles.sexdata={handles.sexdata{:} handlessexdata{:}};

tmp=length(handles.behaviors);
idx=[1 : (sum(cellfun(@length,handles.experimentlist(1:handles.groupvalue)))-length(newrecordings)) ...
    ((tmp-length(newrecordings)+1) : tmp) ...
    ((tmp-length(newrecordings)-sum(cellfun(@length,handles.experimentlist((handles.groupvalue+1):end)))+1) : ...
        (tmp-length(newrecordings)))];
handles.behaviors=handles.behaviors(idx);
handles.features=handles.features(idx);
handles.individuals=handles.individuals(idx);
handles.sexdata=handles.sexdata(idx);

handles.behaviorlist=check_for_diff_and_return_intersection(handles.behaviors);
set(handles.BehaviorList,'String',handles.behaviorlist);
set(handles.BehaviorList2,'String',handles.behaviorlist);

handles.featurelist=check_for_diff_and_return_intersection(handles.features);
set(handles.FeatureList,'String',handles.featurelist);

handles=fillin_individuallist(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');


% --- Executes on button press in ExperimentDelete.
function ExperimentDelete_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

if(length(handles.experimentlist)==0)  return;  end

set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

idx=handles.experimentvalue{handles.groupvalue};
handles.experimentlist{handles.groupvalue}(idx)=[];
handles.experimentvalue{handles.groupvalue}=1:length(handles.experimentlist{handles.groupvalue});
if(isempty(handles.experimentlist{handles.groupvalue}))
  handles.experimentlist(handles.groupvalue)=[];
  handles.experimentvalue(handles.groupvalue)=[];
  handles.grouplist(handles.groupvalue)=[];
  handles.groupvalue=max(1,min([handles.groupvalue length(handles.grouplist)]));
  if(isempty(handles.grouplist))
    set(handles.GroupList,'enable','off');
    set(handles.ExperimentAdd,'enable','off');
    set(handles.ExperimentDelete,'enable','off');
    set(handles.ExperimentMove,'enable','off');
  end
  if(sum(cellfun(@length,handles.experimentlist(:)))==0)
    set(handles.BehaviorList,'enable','off');
    set(handles.BehaviorLogic,'enable','off');
    set(handles.BehaviorList2,'enable','off');
    set(handles.FeatureList,'enable','off');
    set(handles.IndividualList,'enable','off');
  end
end
if(isempty(handles.experimentlist))
  set(handles.ExperimentList,'String',{''},'Value',1);
  set(handles.GroupList,'String',{''},'Value',1);
else
  set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue},...
      'Value',handles.experimentvalue{handles.groupvalue});
  tmp=length(handles.grouplist);
  cellstr(strcat(repmat('<html><font color="',tmp,1),...
      {handles.colors{2,1+mod(0:length(handles.grouplist)-1,length(handles.colors))}}',...
      repmat('">',tmp,1),handles.grouplist',repmat('</font></html>',tmp,1)));
  set(handles.GroupList,'String',ans,'Value',handles.groupvalue);
end

idx=idx+sum(cellfun(@length,handles.experimentlist(1:(handles.groupvalue-1))));

handles.behaviors(idx)=[];
handles.behaviorlist=unique([handles.behaviors{:}]);
handles.behaviorvalue=max(1,min(handles.behaviorvalue,length(handles.behaviorlist)));
handles.behaviorvalue2=max(1,min(handles.behaviorvalue2,length(handles.behaviorlist)));
if(isempty(handles.behaviorlist))
  handles.behaviorlist={''};
end
set(handles.BehaviorList,'String',handles.behaviorlist,'Value',handles.behaviorvalue);
set(handles.BehaviorList2,'String',handles.behaviorlist,'Value',handles.behaviorvalue2);

handles.features(idx)=[];
handles.featurelist=unique([handles.features{:}]);
handles.featurevalue=max(1,min(handles.featurevalue,length(handles.featurelist)));
if(isempty(handles.featurelist))
  handles.featurelist={''};
end
set(handles.FeatureList,'String',handles.featurelist,'Value',handles.featurevalue);

handles.individuals(idx)=[];
handles=fillin_individuallist(handles);

handles.sexdata(idx)=[];

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');


% --- Executes on button press in ExperimentMove.
function ExperimentMove_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentMove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if(isempty(handles.experimentlist{handles.groupvalue}))  return;  end
[to_group,ok]=listdlg('liststring',...
    {handles.grouplist{setdiff(1:length(handles.grouplist),handles.groupvalue)}},...
    'selectionmode','single');
if(~ok)  return;  end
if(to_group>=handles.groupvalue)  to_group=to_group+1;  end

set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

from_group=handles.groupvalue;
idx=handles.experimentvalue{from_group};
idxF=idx+sum(cellfun(@length,handles.experimentlist(1:(from_group-1))));
idxT=sum(cellfun(@length,handles.experimentlist(1:to_group)));

handles.experimentlist{to_group}=...
    {handles.experimentlist{to_group}{:} handles.experimentlist{from_group}{idx}};
handles.experimentvalue{to_group}=1:length(handles.experimentlist{to_group});
handles.experimentlist{from_group}(idx)=[];
handles.experimentvalue{from_group}=1:length(handles.experimentlist{from_group});
if(isempty(handles.experimentlist{from_group}))
  set(handles.ExperimentList,'String',{''},'Value',1);
else
  set(handles.ExperimentList,'String',handles.experimentlist{from_group});
  set(handles.ExperimentList,'Value',handles.experimentvalue{from_group});
end

if(idxF(1)<idxT)  idxT=idxT-length(idxF);  end
tmp=setdiff(1:length(handles.behaviors),idxF);
tmp=[tmp(1:idxT) idxF tmp((idxT+1):end)];
handles.behaviors=handles.behaviors(tmp);
handles.features=handles.features(tmp);
handles.individuals=handles.individuals(tmp);
handles.sexdata=handles.sexdata(tmp);

handles=fillin_individuallist(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');


% --- Executes on selection change in GroupList.
function GroupList_Callback(hObject, eventdata, handles)
% hObject    handle to GroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns GroupList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from GroupList

handles.groupvalue=get(handles.GroupList,'Value');
set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue},...
    'Value',handles.experimentvalue{handles.groupvalue});
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function GroupList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GroupList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in GroupNew.
function GroupNew_Callback(hObject, eventdata, handles)
% hObject    handle to GroupNew (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%char(inputdlg({'Name:'},'Create new experiment group'));
%handles.grouplist{end+1}=['<html><font color="' ...
%    handles.colors{2,1+mod(length(handles.grouplist),length(handles.colors))} '">' ans '</font></html>'];
;
inputdlg({'Name:'},'Create new experiment group');
if(isempty(ans))  return;  end
handles.grouplist{end+1}=char(ans);
handles.groupvalue=length(handles.grouplist);
handles.experimentlist{handles.groupvalue}={};
handles.experimentvalue{handles.groupvalue}=[];
tmp=length(handles.grouplist);
cellstr(strcat(repmat('<html><font color="',tmp,1),...
    {handles.colors{2,1+mod(0:length(handles.grouplist)-1,length(handles.colors))}}',...
    repmat('">',tmp,1),handles.grouplist',repmat('</font></html>',tmp,1)));
set(handles.GroupList,'String',ans,'Value',handles.groupvalue);
set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue},...
    'Value',handles.experimentvalue{handles.groupvalue});
set(handles.GroupList,'enable','on');
set(handles.ExperimentAdd,'enable','on');
set(handles.ExperimentDelete,'enable','on');
set(handles.ExperimentMove,'enable','on');
guidata(hObject,handles);


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
function print_csv_help(fid,type,tstr,xstr,ystr)

tstr=strrep(tstr,'%','%%');
xstr=strrep(xstr,'%','%%');
ystr=strrep(ystr,'%','%%');

fprintf(fid,['%% type=' type '\n%% title=' tstr '\n%% xlabel=' xstr '\n%% ylabel=' ystr '\n%%\n']);

%for g=1:length(handles.grouplist)
%  fprintf(fid,['%% xdata, group ' handles.grouplist{g} ...
%             '\n%% ydata1, group ' handles.grouplist{g} ...
%             '\n%% ydata2, group ' handles.grouplist{g} ...
%             '\n%%...\n%% ydataN, group ' handles.grouplist{g} '\n']);
%  if(g<length(handles.grouplist));  fprintf(fid,'%%\n');  end
%end
%fprintf(fid,'\n');


% ---
function plot_it(xdata,ydata,style,centraltendency,dispersion,color,linewidth,fid,experimentlist)

fprintf(fid,'%% xdata\n');  fprintf(fid,'%g, ',xdata);  fprintf(fid,'\n');
if(style~=3)
  switch(centraltendency)
    case 1
      data_ct=nanmean(ydata,1);
      str_ct='mean';
    case 2
      data_ct=nanmedian(ydata,1);
      str_ct='median';
    case 3
      data_ct=mode(ydata,1);
      str_ct='mode';
  end
end
switch(style)
  case 1
    plot(xdata,data_ct,color,'linewidth',linewidth);
    fprintf(fid,['%% ydata, ' str_ct '\n']);  fprintf(fid,'%g, ',data_ct);  fprintf(fid,'\n');
  case 2
    switch(dispersion)
      case 1
        tmp=nanstd(ydata,[],1);
        data_dp=data_ct+tmp;
        data_dn=data_ct-tmp;
        str_dp=[str_ct ' + std dev'];
        str_dn=[str_ct ' - std dev'];
      case 2
        tmp=nanstd(ydata,[],1)./sqrt(sum(~isnan(ydata),1));
        data_dp=data_ct+tmp;
        data_dn=data_ct-tmp;
        str_dp=[str_ct ' + std err'];
        str_dn=[str_ct ' - std err'];
      case 3
        data_dp=prctile(ydata,95);
        data_dn=prctile(ydata,5);
        str_dp='95%%';
        str_dn='5%%';
      case 4
        data_dp=prctile(ydata,75);
        data_dn=prctile(ydata,25);
        str_dp='75%%';
        str_dn='25%%';
    end
    %plot(xdata,data_dp,color,'linewidth',linewidth);
    %plot(xdata,data_dn,color,'linewidth',linewidth);
    %plot(xdata,data_ct,color,'linewidth',3*linewidth);
    h=plot(xdata,data_ct,color,'linewidth',3*linewidth);
    idx=isnan(data_dp) | isnan(data_dn);
    xdata=xdata(~idx);  data_dp=data_dp(~idx);  data_dn=data_dn(~idx);
    color2=(get(h,'color')+[4 4 4])/5;
    k=1;  m=0;  step=10000;
    while(k<=length(xdata))
      idx=k:min(k+step,length(xdata));
      patch([xdata(idx) fliplr(xdata(idx))],[data_dp(idx) fliplr(data_dn(idx))],...
            color2,'edgecolor','none');
      k=k+step+1;  m=m+1;
    end
    get(gca,'children');  set(gca,'children',circshift(ans,-m));  % send to back
    %get(gca,'children');  set(gca,'children',ans([m+1 1:m (m+2):end]));  % send to back
    fprintf(fid,['%% ydata, ' str_dp '\n']);  fprintf(fid,'%g, ',data_dp);  fprintf(fid,'\n');
    fprintf(fid,['%% ydata, ' str_dn '\n']);  fprintf(fid,'%g, ',data_dn);  fprintf(fid,'\n');
    fprintf(fid,['%% ydata, ' str_ct '\n']);  fprintf(fid,'%g, ',data_ct);  fprintf(fid,'\n');
  case 3
    plot(xdata,ydata',color,'linewidth',linewidth);
    for e=1:size(ydata,1)
      fprintf(fid,['%% ydata, exp ' experimentlist{e} '\n']);
      fprintf(fid,'%g, ',ydata(e,:));  fprintf(fid,'\n');
    end
end
fprintf(fid,'\n');


% --- 
function [during not_during]=calculate_feature_histogram(behavior_data,behavior_logic,behavior_data2,...
    feature_data,sexdata,individual,perwhat)

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  %for i=1:length(feature_data.data)
  %  feature_data.data{i}=strcmp(feature_data.data{i},vals{1});
  %end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

during={};  not_during={};
for i=1:length(behavior_data.allScores.t0s)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end
%  tmp1=false(1,length(feature_data.data{i}));
%  for j=1:length(behavior_data.allScores.t0s{i})
%    tmp1((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
%        -behavior_data.allScores.tStart(i)+1)=true;
%  end
  tmp1=zeros(1,length(feature_data.data{i}));
  tmp1(behavior_data.allScores.t0s{i}-behavior_data.allScores.tStart(i)+1)=1;
  tmp1(behavior_data.allScores.t1s{i}-behavior_data.allScores.tStart(i)+1)=-1;
  tmp1=logical(cumsum(tmp1(1:length(feature_data.data{i}))));
  
  if(behavior_logic>1)
    %tmp2=false(1,length(feature_data.data{i}));
    %for j=1:length(behavior_data2.allScores.t0s{i})
    %  tmp2((behavior_data2.allScores.t0s{i}(j):(behavior_data2.allScores.t1s{i}(j)-1))...
    %      -behavior_data2.allScores.tStart(i)+1)=true;
    %end
    tmp2=zeros(1,length(feature_data.data{i}));
    tmp2(behavior_data2.allScores.t0s{i}-behavior_data2.allScores.tStart(i)+1)=1;
    tmp2(behavior_data2.allScores.t1s{i}-behavior_data2.allScores.tStart(i)+1)=-1;
    tmp2=logical(cumsum(tmp2(1:length(feature_data.data{i}))));
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
  if length(partition_idx)==(length(feature_data.data{i})+1)  % n-1 features
    partition_idx=partition_idx(1:end-1);
  end
  if(perwhat==1)  % per frame
    during{i}=feature_data.data{i}(partition_idx & sexdata{i}(1:length(partition_idx)));
    not_during{i}=feature_data.data{i}((~partition_idx) & sexdata{i}(1:length(partition_idx)));
  else  % per bout
    partition_idx=[0 partition_idx 0];
    start=1+find(~partition_idx(1:(end-1)) &  partition_idx(2:end))-1;
    stop =  find( partition_idx(1:(end-1)) & ~partition_idx(2:end))-1;
    during{i}=nan;  not_during{i}=nan;
    if(length(start)>0)
      for j=1:length(start)
        if(sum(sexdata{i}(start(j):stop(j))) < ((stop(j)-start(j)+1)/2))  continue;  end
        switch(perwhat)
          case 2
            during{i}(j)=mean(feature_data.data{i}(start(j):stop(j)));
          case 3
            during{i}(j)=median(feature_data.data{i}(start(j):stop(j)));
          case 4
            during{i}(j)=max(feature_data.data{i}(start(j):stop(j)));
          case 5
            during{i}(j)=min(feature_data.data{i}(start(j):stop(j)));
          case 6
            during{i}(j)=std(feature_data.data{i}(start(j):stop(j)));
        end
      end
      if(length(start)>1)
        for j=1:(length(start)-1)
          if(sum(sexdata{i}(stop(j):start(j+1))) < ((start(j+1)-stop(j)+1)/2))  continue;  end
          switch(perwhat)
            case 2
              not_during{i}(j)=mean(feature_data.data{i}(stop(j):start(j+1)));
            case 3
              not_during{i}(j)=median(feature_data.data{i}(stop(j):start(j+1)));
            case 4
              not_during{i}(j)=max(feature_data.data{i}(stop(j):start(j+1)));
            case 5
              not_during{i}(j)=min(feature_data.data{i}(stop(j):start(j+1)));
            case 6
              not_during{i}(j)=std(feature_data.data{i}(stop(j):start(j+1)));
          end
        end
      end
    end
  end
end

during=[during{:}];
not_during=[not_during{:}];


% ---
%function [table_data feature_units]=plot_feature_histogram(experiment_value,experiment_list,...
%function [feature_units]=plot_feature_histogram(experiment_value,experiment_list,...
function plot_feature_histogram(experiment_value,experiment_list,...
    behavior_value,behavior_list,behavior_logic,behavior_value2,...
    feature_value,feature_list,individual,sexdata,perwhat,style,notduring,logbinsize,nbins,...
    centraltendency,dispersion,color,fid)

during=cell(1,length(experiment_value));
not_during=cell(1,length(experiment_value));
bad=zeros(1,length(experiment_value));
parfor e=1:length(experiment_value)
%for e=1:length(experiment_value)
  behavior_data=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list{behavior_value} '.mat']));
  if(behavior_logic>1)
    behavior_data2=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list{behavior_value2} '.mat']));
  else
    behavior_data2=[];
  end
  feature_data=load(fullfile(experiment_list{experiment_value(e)},'perframe',...
        [feature_list{feature_value} '.mat']));

  if((length(behavior_data.allScores.scores)~=length(feature_data.data)) || ...
      ((behavior_logic>1) && (length(behavior_data2.allScores.scores)~=length(feature_data.data))))
    bad(e)=1;
    continue;
  end

  tmp2=sexdata{experiment_value(e)};
  for i=1:length(tmp2)
    switch(individual)
      case(2)
      case(3)
        tmp2{i}=~tmp2{i};
      otherwise
        tmp2{i}=ones(1,length(tmp2{i}));
    end
  end
  tmp=nan;  if(individual>3)  tmp=individual-3;  end

  [during{e} not_during{e}]=calculate_feature_histogram(behavior_data,behavior_logic,behavior_data2,...
      feature_data,tmp2,tmp,perwhat);
end

tmp=find(bad);
if(length(tmp)>0)
  during(tmp)=[];  not_during(tmp)=[];
  msg{1}=['the following experiments had different numbers of individuals ' ...
      'for the selected behavior and feature.  they will be ommitted from the analysis.'];
  msg{2}='';
  for i=1:length(tmp)
    [~,msg{end+1},~]=fileparts(experiment_list{i});
  end
  uiwait(errordlg(msg));
end

max(cellfun(@(x) size(x,2),during));
cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],during,'uniformoutput',false);
during=cat(1,ans{:});
max(cellfun(@(x) size(x,2),not_during));
cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],not_during,'uniformoutput',false);
not_during=cat(1,ans{:});

if(logbinsize)
  if(notduring)
    unique([reshape(abs(during),1,prod(size(during))) reshape(abs(not_during),1,prod(size(not_during)))]);
    nearzero=ans(1);  if(ans(1)==0)  nearzero=ans(2);  end
    low=min(min([during not_during]));
    high=max(max([during not_during]));
  else
    unique(reshape(abs(during),1,prod(size(during))));
    nearzero=ans(1);  if(ans(1)==0)  nearzero=ans(2);  end
    low=min(min(during));
    high=max(max(during));
  end
  if((low>=0) && (high>0))
    tmp=logspace(log10(max(low,nearzero)),log10(high),nbins);
  elseif((low<0) && (high<=0))
    tmp=fliplr(-logspace(log10(max(abs(high),nearzero)),log10(abs(low)),nbins));
  elseif((low<0) && (high>0))
    tmp=[fliplr(-logspace(log10(nearzero),log10(abs(low)),nbins)) ...
        logspace(log10(nearzero),log10(abs(high)),nbins)];
  end
  %set(gca,'xscale','log');
else
  if(notduring)
    tmp=linspace(min(min([during not_during])),max(max([during not_during])),nbins);
  else
    tmp=linspace(min(min(during)),max(max(during)),nbins);
  end
end

if(notduring)
  hist_not_during=hist(not_during',tmp);
  if(size(not_during,1)==1)  hist_not_during=hist_not_during';  end
  hist_not_during.*repmat(([0 diff(tmp)]+[diff(tmp) 0])'/2,1,size(hist_not_during,2));
  hist_not_during=hist_not_during./repmat(sum(ans,1),size(hist_not_during,1),1);
  plot_it(tmp,hist_not_during',style,centraltendency,dispersion,color,1,...
    fid,experiment_list(experiment_value));
end
hist_during=hist(during',tmp);
if(size(during,1)==1)  hist_during=hist_during';  end
hist_during.*repmat(([0 diff(tmp)]+[diff(tmp) 0])'/2,1,size(hist_during,2));
hist_during=hist_during./repmat(sum(ans,1),size(hist_during,1),1);
linewidth=1;  if(notduring)  linewidth=2;  end
plot_it(tmp,hist_during',style,centraltendency,dispersion,color,linewidth,...
    fid,experiment_list(experiment_value));


% --- Executes on button press in FeatureHistogram.
function FeatureHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];

gg=1:length(handles.grouplist);
experimentvalue=handles.experimentvalue;
individual=handles.individualvalue;
if(individual>3)
  tmp=find(cumsum_num_indi_per_exp<(individual-3),1,'last');
  gg=find(cumsum_num_exp_per_group<tmp,1,'last');
  experimentvalue{gg}=tmp-cumsum_num_exp_per_group(gg);
  individual=individual-cumsum_num_indi_per_exp(tmp);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='feature histogram';
%xlabel(get_label(handles.featurelist(handles.featurevalue),feature_units{1}));
tstr=char(strrep(handles.behaviorlist(handles.behaviorvalue),'_','-'));
switch(handles.behaviorlogic)
  case 2
    tstr=[tstr ' AND '];
  case 3
    tstr=[tstr ' AND NOT '];
  case 4
    tstr=[tstr ' OR '];
end
if(handles.behaviorlogic>1)
  tstr=[tstr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
end
%xstr=[xstr ' (%)'];
units=load(fullfile(handles.experimentlist{gg(1)}{experimentvalue{gg(1)}(1)},'perframe',...
    [handles.featurelist{handles.featurevalue} '.mat']),'units');
xstr=get_label(handles.featurelist(handles.featurevalue),units.units);
ystr='normalized';

print_csv_help(fid,handles.type,tstr,xstr,ystr);

guidata(figure('ButtonDownFcn',@ButtonDownFcn_Callback),handles);  hold on;

for g=gg
  set(handles.Status,'string',...
      ['Processing ' num2str(length(handles.experimentvalue{g})) ' experiment(s) in group ' handles.grouplist{g}]);
  drawnow;
  if(~isempty(experimentvalue{g}))
    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
    plot_feature_histogram(experimentvalue{g},handles.experimentlist{g},...
        handles.behaviorvalue,handles.behaviorlist,handles.behaviorlogic,handles.behaviorvalue2,...
        handles.featurevalue,handles.featurelist,individual,...
        handles.sexdata((cumsum_num_exp_per_group(g)+1):cumsum_num_exp_per_group(g+1)),...
        handles.featurehistogram_perwhat,handles.featurehistogram_style,...
        handles.featurehistogram_notduring,handles.featurehistogram_logbinsize,handles.featurehistogram_nbins,...
        handles.prefs_centraltendency,handles.prefs_dispersion,...
        handles.colors{1,1+mod(g-1,length(handles.colors))},fid);
  end
end

fclose(fid);

title(tstr);
xlabel(xstr);
ylabel(ystr);
axis tight;  zoom reset;

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');

%guidata(hObject,handles);


% ---
%function table_data=calculate_interesting_feature_histograms(experiment_value,experiment_list,...
%    behavior_list,feature_list,perwhat)
%
%table_data=zeros(length(behavior_list),length(feature_list),8);
%parfor b=1:length(behavior_list)
%%for b=1:length(behavior_list)
%  k=1;
%  parfor_tmp=zeros(length(feature_list),8);
%  for f=1:length(feature_list)
%    during={};  not_during={};
%    for e=1:length(experiment_value)
%      behavior_data=load(fullfile(experiment_list{experiment_value(e)},[behavior_list{b} '.mat']));
%      feature_data=load(fullfile(experiment_list{experiment_value(e)},'perframe',...
%          [feature_list{f} '.mat']));
%      sexdata={};
%      for i=1:length(feature_data.data)
%        sexdata{i}=ones(1,length(feature_data.data{i}));
%      end
%      [during{e} not_during{e}]=calculate_feature_histogram(behavior_data,1,[],feature_data,sexdata,nan,perwhat);
%    end
%    during=[during{:}];
%    not_during=[not_during{:}];
%    %parfor_tmp(k,:)=[b f (mean(during)-mean(not_during))/sqrt((std(during)^2+std(not_during)^2)/2)];
%    parfor_tmp(k,:)=[b f mean(during) mean(not_during) std(during) std(not_during) length(during) length(not_during)];
%    k=k+1;
%  end
%  table_data(b,:,:)=parfor_tmp;
%  disp([num2str(b) ' of ' num2str(length(behavior_list))]);
%end
%table_data=reshape(table_data,prod(size(table_data))/8,8);


% --- Executes on button press in InterestingFeatureHistograms.
function InterestingFeatureHistograms_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingFeatureHistograms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

if(isempty(handles.interestingfeaturehistograms_cache))
  table_data={};
  bad={};
  for g=1:length(handles.grouplist)
    set(handles.Status,'string',...
        ['Processing ' num2str(length(handles.experimentvalue{g})) ' experiment(s) in group ' handles.grouplist{g}]);
    drawnow;
    if(~isempty(handles.experimentvalue{g}))
%      table_data{end+1}=calculate_interesting_feature_histograms(...
%          handles.experimentvalue{g},handles.experimentlist{g},...
%          handles.behaviorlist,handles.featurelist,handles.featurehistogram_perwhat);
      experiment_value=handles.experimentvalue{g};
      experiment_list=handles.experimentlist{g};
      behavior_list=handles.behaviorlist;
      feature_list=handles.featurelist;
      perwhat=handles.featurehistogram_perwhat;
      parfor_tmp2=zeros(length(behavior_list),length(feature_list),8);
      bad2={};
      parfor b=1:length(behavior_list)
      %for b=1:length(behavior_list)
        behavior_data={};
        for e=1:length(experiment_value)
          behavior_data{e}=load(fullfile(experiment_list{experiment_value(e)},[behavior_list{b} '.mat']));
        end
        k=1;  bad2{b}={};
        parfor_tmp=zeros(length(feature_list),8);
        for f=1:length(feature_list)
          during={};  not_during={};
          for e=1:length(experiment_value)
            feature_data=load(fullfile(experiment_list{experiment_value(e)},'perframe',...
                [feature_list{f} '.mat']));

            if(length(behavior_data{e}.allScores.scores)~=length(feature_data.data))
              [~,foo,~]=fileparts(experiment_list{experiment_value(e)});
              bad2{b}{end+1}=[foo ', ' behavior_list{b} ', ' feature_list{f}];
              continue;
            end

            sexdata={};
            for i=1:length(feature_data.data)
              sexdata{i}=ones(1,length(feature_data.data{i}));
            end
            [during{e} not_during{e}]=calculate_feature_histogram(behavior_data{e},1,[],feature_data,sexdata,nan,perwhat);
          end
          during=[during{:}];
          not_during=[not_during{:}];
          %parfor_tmp(k,:)=[b f (mean(during)-mean(not_during))/sqrt((std(during)^2+std(not_during)^2)/2)];
          parfor_tmp(k,:)=[b f mean(during) mean(not_during) std(during) std(not_during) length(during) length(not_during)];
          k=k+1;
        end
        parfor_tmp2(b,:,:)=parfor_tmp;
        disp([num2str(b) ' of ' num2str(length(behavior_list))]);
      end
      bad=[bad [bad2{:}]];
      table_data{end+1}=reshape(parfor_tmp2,prod(size(parfor_tmp2))/8,8);
    end
  end

  if(~isempty(bad))
    msg{1}=['the following experiments had different numbers of individuals ' ...
        'for the indicated behavior and feature.  they will be ommitted from the analysis.'];
    msg{2}='';
    tmp=length(bad);
    if(tmp>20)
      bad=bad(1:20);
      bad{end+1}='';
      bad{end+1}=['plus another ' num2str(tmp) ' combinations'];
    end
    uiwait(errordlg({msg{:} bad{:}}));
  end

  tmp2=[];
  for g=1:length(handles.grouplist)
    tmp2=[tmp2; repmat(g,size(table_data{g},1),1) nan(size(table_data{g},1),1) table_data{g}(:,1:2) ...
        (table_data{g}(:,3)-table_data{g}(:,4))./sqrt(table_data{g}(:,5).^2+table_data{g}(:,6).^2)...
        table_data{g}(:,7) table_data{g}(:,8)];
    if(g==length(handles.grouplist))  break;  end
    for g2=(g+1):length(handles.grouplist)
      tmp2=[tmp2; repmat(g,size(table_data{g},1),1) repmat(g2,size(table_data{g},1),1) table_data{g}(:,1:2) ...
          (table_data{g2}(:,3)-table_data{g}(:,3))./sqrt(table_data{g2}(:,5).^2+table_data{g}(:,5).^2)...
          table_data{g}(:,7) table_data{g2}(:,7)];
    end
  end
  if(handles.interestingfeaturehistograms_omitnan)
    idx=find(~isnan(tmp2(:,5)));
    tmp2=tmp2(idx,:);
  end
  if(handles.interestingfeaturehistograms_omitinf)
    idx=find(~isinf(tmp2(:,5)));
    tmp2=tmp2(idx,:);
  end
  tmp2=sortrows(tmp2,-5);

  handles.interestingfeaturehistograms_cache=tmp2;
else
  tmp2=handles.interestingfeaturehistograms_cache;
end

tmp=cell(size(tmp2,1),5);
tmp(:,1)=handles.grouplist(tmp2(:,1));
%tmp(:,2)=num2cell(tmp2(:,6));
tmp(:,2)=cellstr(num2str(tmp2(:,6),'%-d'));
idx=~isnan(tmp2(:,2));
tmp(idx,3)=handles.grouplist(tmp2(idx,2));
tmp(:,4)=cellstr(num2str(tmp2(:,7),'%-d'));
tmp(:,5)=handles.behaviorlist(tmp2(:,3));
tmp(:,6)=handles.featurelist(tmp2(:,4));
tmp(:,7)=num2cell(tmp2(:,5));
set(handles.Table,'Data',tmp);
set(handles.Table,'ColumnName',{'Group' 'n' 'Group2' 'n2' 'Behavior' 'Feature' 'd'''});
set(handles.Table,'ColumnWidth',{75 50 75 50 150 100 75});
set(handles.Table,'RowStriping','on','BackgroundColor',[1 1 1; 0.95 0.95 0.95]);

handles.table_data=tmp2;
handles.table='histogram';
guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');


% --- 
function data=calculate_entiretimeseries(behavior_data,feature_data,sexdata,individual);

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  %for i=1:length(feature_data.data)
  %  feature_data.data{i}=strcmp(feature_data.data{i},vals{1});
  %end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

k=1;
data=nan(length(feature_data.data),max(behavior_data.allScores.tEnd));
for i=1:length(feature_data.data)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end
  if(sum(sexdata{i}) < length(sexdata{i})/2)  continue;  end
  tmp=behavior_data.allScores.tEnd(i)-behavior_data.allScores.tStart(i)+1;
  foo=0;
  if tmp==(length(feature_data.data{i})+1)  % n-1 features
    foo=1;
  end
  data(k,behavior_data.allScores.tStart(i):(behavior_data.allScores.tEnd(i)-foo))=...
      feature_data.data{i}(1:(tmp-foo));
  k=k+1;
end


% --- 
function triggered_data=calculate_triggeredtimeseries(behavior_data,behavior_logic,behavior_data2,feature_data,...
    sexdata,individual,timing,windowradius,subtractmean)

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  %for i=1:length(feature_data.data)
  %  feature_data.data{i}=strcmp(feature_data.data{i},vals{1});
  %end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

k=1;
%triggered_data=zeros(length([behavior_data.allScores.t0s{:}]),1+2*windowradius);
triggered_data=[];
for i=1:length(behavior_data.allScores.t0s)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end
  feature_data_padded=[nan(1,windowradius) feature_data.data{i} nan(1,windowradius)];
  sexdata_padded=[nan(1,windowradius) sexdata{i} nan(1,windowradius)];
  if(behavior_logic>1)
    idx2=[behavior_data2.allScores.t0s{i}'-behavior_data2.allScores.tStart(i) ...
       behavior_data2.allScores.t1s{i}'-behavior_data2.allScores.tStart(i)];
  end
  for j=1:length(behavior_data.allScores.t0s{i})  % bout
    switch(timing)
      case(2)
        idx=behavior_data.allScores.t0s{i}(j)-behavior_data.allScores.tStart(i);
      case(3)
        idx=behavior_data.allScores.t1s{i}(j)-behavior_data.allScores.tStart(i);
    end
    if((behavior_logic==2) && (diff(sum(idx2<idx))==0))
      continue;
    end
    if((behavior_logic==3) && (diff(sum(idx2<idx))==-1))
      continue;
    end
    if(sum(sexdata_padded(idx+(0:(2*windowradius))+(3-timing))) < windowradius)  continue;  end
    triggered_data(k,:)=feature_data_padded(idx+(0:(2*windowradius))+(3-timing));
    k=k+1;
  end
end
if(k==1)
  triggered_data(k,:)=nan(1,2*windowradius+1);
end
if(subtractmean)
  triggered_data=triggered_data-...
      repmat(nanmean(triggered_data(:,1:windowradius),2),1,size(triggered_data,2));
end


%---
%function [table_data feature_units range h]=plot_timeseries(experiment_value,experiment_list,...
%function feature_units=plot_feature_timeseries(experiment_value,experiment_list,...
function plot_feature_timeseries(experiment_value,experiment_list,...
    behavior_value,behavior_list,behavior_logic,behavior_value2,feature_value,feature_list,...
    individual,sexdata,timing,style,centraltendency,dispersion,convolutionwidth,subtractmean,windowradius,...
    color,fid)

data=cell(1,length(experiment_value));
bad=zeros(1,length(experiment_value));
parfor e=1:length(experiment_value)
  behavior_data=load(fullfile(experiment_list{experiment_value(e)},...
      [behavior_list{behavior_value} '.mat']));
  if(behavior_logic>1)
    behavior_data2=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list{behavior_value2} '.mat']));
  else
    behavior_data2=[];
  end
  feature_data=load(fullfile(experiment_list{experiment_value(e)},'perframe',...
      [feature_list{feature_value} '.mat']));
  %feature_units{e}=feature_data.units;

  if((length(behavior_data.allScores.scores)~=length(feature_data.data)) || ...
      ((behavior_logic>1) && (length(behavior_data2.allScores.scores)~=length(feature_data.data))))
    bad(e)=1;
    continue;
  end

  tmp2=sexdata{experiment_value(e)};
  for i=1:length(tmp2)
    switch(individual)
      case(2)
      case(3)
        tmp2{i}=~tmp2{i};
      otherwise
        tmp2{i}=ones(1,length(tmp2{i}));
    end
  end
  tmp=nan;  if(individual>3)  tmp=individual-3;  end

  if(timing==1)
    calculate_entiretimeseries(behavior_data,feature_data,tmp2,tmp);
    conv(nanmean(ans,1),ones(1,convolutionwidth),'same');
    data{e}=ans./conv(ones(1,length(ans)),ones(1,convolutionwidth),'same');
  else
    calculate_triggeredtimeseries(behavior_data,behavior_logic,behavior_data2,...
        feature_data,tmp2,tmp,timing,windowradius,subtractmean);
    data{e}=nanmean(ans,1);
  end
end

tmp=find(bad);
if(length(tmp)>0)
  data(tmp)=[];
  msg{1}=['the following experiments had different numbers of individuals ' ...
      'for the selected behavior and feature.  they will be ommitted from the analysis.'];
  msg{2}='';
  for i=1:length(tmp)
    [~,msg{end+1},~]=fileparts(experiment_list{i});
  end
  uiwait(errordlg(msg));
end

if(timing==1)
  max(cellfun(@(x) size(x,2),data));
  cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],data,'uniformoutput',false);
  ydata=cat(1,ans{:});
  xdata=1:size(ydata,2);
  %table_data=[];
else
  ydata=cat(1,data{:});
  xdata=-windowradius:windowradius;
  %table_data=[sqrt(nanmean(ydata(:,1:windowradius).^2,2))...
  %            sqrt(nanmean(ydata(:,(windowradius+1):end).^2,2))];
end

%feature_units=feature_units{1};

plot_it(xdata,ydata,style,centraltendency,dispersion,color,1,fid,experiment_list(experiment_value));


% --- Executes on button press in FeatureTimeSeries.
function FeatureTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];

gg=1:length(handles.grouplist);
experimentvalue=handles.experimentvalue;
individual=handles.individualvalue;
if(individual>3)
  tmp=find(cumsum_num_indi_per_exp<(individual-3),1,'last');
  gg=find(cumsum_num_exp_per_group<tmp,1,'last');
  experimentvalue{gg}=tmp-cumsum_num_exp_per_group(gg);
  individual=individual-cumsum_num_indi_per_exp(tmp);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='feature time series';
tstr='';
xstr='time (frames)';
units=load(fullfile(handles.experimentlist{gg(1)}{experimentvalue{gg(1)}(1)},'perframe',...
    [handles.featurelist{handles.featurevalue} '.mat']),'units');
%ystr=get_label(handles.featurelist(handles.featurevalue),feature_units{1});
ystr=get_label(handles.featurelist(handles.featurevalue),units.units);
if(handles.featuretimeseries_timing>1)
  tstr=char(strrep(handles.behaviorlist(handles.behaviorvalue),'_','-'));
  switch(handles.behaviorlogic)
    case 2
      tstr=[tstr ' AND '];
    case 3
      tstr=[tstr ' AND NOT '];
    case 4
      tstr=[tstr ' OR '];
  end
  if(handles.behaviorlogic>1)
    tstr=[tstr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
  end
  tstr=[tstr ' (%)'];
end

print_csv_help(fid,handles.type,tstr,xstr,ystr);

guidata(figure('ButtonDownFcn',@ButtonDownFcn_Callback),handles);  hold on;

%range=[];
table_data={};  %feature_units={};
for g=gg
  set(handles.Status,'string',...
      ['Processing ' num2str(length(handles.experimentvalue{g})) ' experiment(s) in group ' handles.grouplist{g}]);
  drawnow;
  if(~isempty(experimentvalue{g}))
    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
    %[table_data{end+1} feature_units{end+1} tmp h]=plot_timeseries(experimentvalue{g},handles.experimentlist{g},...
    %feature_units{end+1}=plot_feature_timeseries(experimentvalue{g},handles.experimentlist{g},...
    plot_feature_timeseries(experimentvalue{g},handles.experimentlist{g},...
        handles.behaviorvalue,handles.behaviorlist,handles.behaviorlogic,handles.behaviorvalue2,...
        handles.featurevalue,handles.featurelist,individual,...
        handles.sexdata((cumsum_num_exp_per_group(g)+1):cumsum_num_exp_per_group(g+1)),...
        handles.featuretimeseries_timing,handles.featuretimeseries_style,...
        handles.prefs_centraltendency,handles.prefs_dispersion,handles.prefs_convolutionwidth,...
        handles.featuretimeseries_subtractmean,handles.timeseries_windowradius,...
        handles.colors{1,1+mod(g-1,length(handles.colors))},fid);
%    range=[range; tmp];
  end
end

fclose(fid);

xlabel(xstr);
ylabel(ystr);
%title(tstr);
axis tight;  zoom reset;



%if((handles.timeseries_tight==1) && (min(range(:,1))<max(range(:,2))))
%  v=axis;  axis([v(1) v(2) min(range(:,1)) max(range(:,2))]);
%end

%if(handles.featuretimeseries_timing>1)
%  tmp=cellfun(@(x) [nanmean(x(:,1)) nanstd(x(:,1)) nanmean(x(:,2)) nanstd(x(:,2))],table_data,...
%      'uniformoutput',false);
%  tmp=cat(1,tmp{:});
%
%  {'mean' 'std. dev.' 'mean' 'std. dev.'};
%  tmp2(1,:)=['RMS before,after:' sprintf('%10s ',ans{:})];
%  for i=1:size(tmp,1)
%    tmp2(i+1,:)=['                 ' sprintf('%10.3g ',tmp(i,:))];
%  end
%  v=axis;
%  h=text(v(1),v(4),tmp2,'color',[0 0.5 0],'tag','stats','verticalalignment','top','fontname','fixed');
%  if(handles.stats)
%    set(h,'visible','on');
%  else
%    set(h,'visible','off');
%  end
%end

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');


% ---
function table_data=calculate_interesting_timeseries(experiment_value,experiment_list,...
    behavior_list,feature_list,statistic,windowradius)

table_data=cell(length(behavior_list),length(feature_list),2);
parfor b=1:length(behavior_list)
%for b=1:length(behavior_list)
  parfor_tmp=cell(length(feature_list),2);
  for f=1:length(feature_list)
    for t=2:3  % timing
      parfor_tmp{f,t-1}=[];
      for e=1:length(experiment_value)
        behavior_data=load(fullfile(experiment_list{experiment_value(e)},[behavior_list{b} '.mat']));
        feature_data=load(fullfile(experiment_list{experiment_value(e)},'perframe',...
          [feature_list{f} '.mat']));
        sexdata={};
        for i=1:length(feature_data.data)
          sexdata{i}=ones(1,length(feature_data.data{i}));
        end

        tmp=calculate_triggeredtimeseries(behavior_data,1,[],feature_data,sexdata,nan,t,windowradius,1);
        parfor_tmp{f,t-1}=[parfor_tmp{f,t-1}; tmp];
      end
      parfor_tmp{f,t-1}=[...
          sqrt(nanmean(parfor_tmp{f,t-1}(:,1:windowradius).^2,2))...
          sqrt(nanmean(parfor_tmp{f,t-1}(:,(windowradius+1):end).^2,2))];
    end
  end
  table_data(b,:,:)=parfor_tmp;
  disp([num2str(b) ' of ' num2str(length(behavior_list))]);
end


% --- Executes on button press in InterestingTimeSeries.
function InterestingTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_list=get(handles.BehaviorList,'String');
feature_list=get(handles.FeatureList,'String');
statistic=handles.prefsstat;
windowradius=handles.timeseries_windowradius;

if(isempty(handles.interestingfeaturetimeseries_cache))
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

  handles.interestingfeaturetimeseries_cache=tmp2;
else
  tmp2=handles.interestingfeaturetimeseries_cache;
end

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

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');


% ---
function ret_val=print_stats(data,stat)

switch(stat)
  case(0)
    ret_val=num2str(data,3);
  case(1)
    ret_val=num2str(mean(data),3);
  case(2)
    ret_val=[num2str(mean(data),3) ' (' num2str(std(data),3) ')'];
  case(3)
    ret_val=[num2str(mean(data),3) ' (' num2str(std(data)./sqrt(length(data)),3) ')'];
  case(4)
    ret_val=num2str(median(data),3);
  case(5)
    ret_val=[num2str(median(data),3) ' (' num2str(prctile(data,5),3) '-' num2str(prctile(data,95),3) ')'];
  case(6)
    ret_val=[num2str(median(data),3) ' (' num2str(prctile(data,25),3) '-' num2str(prctile(data,75),3) ')'];
end


% ---
function [ct,dp,dn]=calculate_ct_d(data,centraltendency,dispersion)

switch(centraltendency)
  case 1
    ct=mean(data);
  case 2
    ct=median(data);
  case 3
    ct=mode(data);
end
switch(dispersion)
  case 1
    tmp=std(data);
    dp=ct+tmp;
    dn=ct-tmp;
  case 2
    tmp=std(data)./sqrt(length(data));
    dp=ct+tmp;
    dn=ct-tmp;
  case 3
    dp=prctile(data,95);
    dn=prctile(data,5);
  case 4
    dp=prctile(data,75);
    dn=prctile(data,25);
end


% ---
function [frames_labelled frames_total]=calculate_behavior_barchart(experiment_value,experiment_list,...
    behavior_list,behavior_logic,behavior_value2,individual,sexdata,perwhat)

collated_data=cell(length(experiment_value),length(behavior_list));
bad=zeros(1,length(experiment_value));
parfor e=1:length(experiment_value)
  behavior_data=load(fullfile(experiment_list{experiment_value(e)},[behavior_list{1} '.mat']));
  if(behavior_logic>1)
    behavior_data2=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list{behavior_value2} '.mat']));
  else
    behavior_data2=[];
  end

  if((length(behavior_data.allScores.scores)~=length(sexdata{experiment_value(e)})) || ...
      ((behavior_logic>1) && (length(behavior_data2.allScores.scores)~=length(sexdata{experiment_value(e)}))))
    bad(e)=1;
    continue;
  end

  parfor_tmp=cell(1,length(behavior_list));
  for b=1:length(behavior_list)
    behavior_data=load(fullfile(experiment_list{experiment_value(e)},[behavior_list{b} '.mat']));

    frames_labelled=nan(1,length(behavior_data.allScores.t0s));
    frames_total=nan(1,length(behavior_data.allScores.t0s));
    sex=nan(1,length(behavior_data.allScores.t0s));
    for i=1:length(behavior_data.allScores.t0s)  % individual
      tmp1=false(1,length(sexdata{experiment_value(e)}{i}));
      for j=1:length(behavior_data.allScores.t0s{i})  % bout
        tmp1((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
            -behavior_data.allScores.tStart(i)+1)=true;
      end
      tmp2=[];
      if(behavior_logic>1)
        tmp2=false(1,length(sexdata{experiment_value(e)}{i}));
        for j=1:length(behavior_data2.allScores.t0s{i})  % bout
          tmp2((behavior_data2.allScores.t0s{i}(j):(behavior_data2.allScores.t1s{i}(j)-1))...
              -behavior_data2.allScores.tStart(i)+1)=true;
        end
      end
      partition_idx=[];
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
      sex(i)=sum(sexdata{experiment_value(e)}{i}(1:length(partition_idx))) > (length(partition_idx)/2);
      frames_labelled(i)=sum(partition_idx);
      frames_total(i)=length(partition_idx);
    end

    parfor_tmp{b}={frames_labelled frames_total sex};
  end
  collated_data(e,:)=parfor_tmp;
end

tmp=find(bad);
if(length(tmp)>0)
  collated_data=collated_data(~logical(bad),:);
  msg{1}=['the following experiments had different numbers of individuals ' ...
      'for the selected behavior and sex feature.  they will be ommitted from the analysis.'];
  msg{2}='';
  for i=1:length(tmp)
    [~,msg{end+1},~]=fileparts(experiment_list{i});
  end
  uiwait(errordlg(msg));
end

switch(individual)
  case 1
    frames_labelled=cellfun(@(x) x{1},collated_data,'uniformoutput',false);
    frames_total=cellfun(@(x) x{2},collated_data,'uniformoutput',false);
  case {2,3}
    frames_labelled=cellfun(@(x) x{1}(x{3}==(3-individual)),collated_data,'uniformoutput',false);
    frames_total=cellfun(@(x) x{2}(x{3}==(3-individual)),collated_data,'uniformoutput',false);
  otherwise
    frames_labelled=cellfun(@(x) x{1}(individual-3),collated_data,'uniformoutput',false);
    frames_total=cellfun(@(x) x{2}(individual-3),collated_data,'uniformoutput',false);
end


% --- Executes on button press in BehaviorBarChart.
function BehaviorBarChart_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorBarChart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];

gg=1:length(handles.grouplist);
experimentvalue=handles.experimentvalue;
individual=handles.individualvalue;
if(individual>3)
  tmp=find(cumsum_num_indi_per_exp<(individual-3),1,'last');
  gg=find(cumsum_num_exp_per_group<tmp,1,'last');
  experimentvalue{gg}=tmp-cumsum_num_exp_per_group(gg);
  individual=individual-cumsum_num_indi_per_exp(tmp);
end

frames_labelled={};  frames_total={};  color={};  xticklabels={};
for g=gg
  set(handles.Status,'string',...
      ['Processing ' num2str(length(handles.experimentvalue{g})) ' experiment(s) in group ' handles.grouplist{g}]);
  drawnow;
  if(~isempty(handles.experimentvalue{g}))
    [frames_labelled{end+1} frames_total{end+1}]=calculate_behavior_barchart(...
        experimentvalue{g},handles.experimentlist{g},...
        handles.behaviorlist,handles.behaviorlogic,handles.behaviorvalue2,individual,...
        handles.sexdata((cumsum_num_exp_per_group(g)+1):cumsum_num_exp_per_group(g+1)),...
        handles.behaviorbarchart_perwhat);
    color{end+1}=handles.colors{1,1+mod(gg(g)-1,length(handles.colors))};
    xticklabels{end+1}=handles.grouplist{g};
  end
end

fid=fopen('most_recent_figure.csv','w');

handles.type='behavior bar chart';
guidata(figure('ButtonDownFcn',@ButtonDownFcn_Callback),handles);  hold on;

table_data=cell(1,length(frames_labelled{1}));
for b=1:length(frames_labelled{1})  % behavior
  tstr='';
  ystr=char(strrep(handles.behaviorlist(b),'_','-'));
  switch(handles.behaviorlogic)
    case 2
      ystr=[ystr ' AND '];
    case 3
      ystr=[ystr ' AND NOT '];
    case 4
      ystr=[ystr ' OR '];
  end
  if(handles.behaviorlogic>1)
    ystr=[ystr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
  end
  ystr=[ystr ' (%)'];
  xstr='group';

  print_csv_help(fid,handles.type,tstr,xstr,ystr);

  ceil(sqrt(length(frames_labelled{1})));
  subplot(ceil(length(frames_labelled{1})/ans),ans,b);  hold on;
  k=[];
  switch(handles.behaviorbarchart_perwhat)
    case 1  % per group
      table_data{b}=nan(1,length(frames_labelled));
      for g=1:length(frames_labelled)
        table_data{b}(g)=100*sum([frames_labelled{g}{:,b}])./sum([frames_total{g}{:,b}]);
        bar(g,table_data{b}(g),color{g});
      end
      fprintf(fid,['%% xdata\n']);  fprintf(fid,'%s, ',handles.grouplist{gg});  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, per group\n']);  fprintf(fid,'%g, ',table_data{b});  fprintf(fid,'\n');
    case 2  % per experiment, error bars
      table_data{b}=cell(1,length(frames_labelled));
      for g=1:length(frames_labelled)
        table_data{b}{g}=100*cellfun(@sum,frames_labelled{g}(:,b))./cellfun(@sum,frames_total{g}(:,b));
        [ct{g},dp{g},dn{g}]=...
            calculate_ct_d(table_data{b}{g},handles.prefs_centraltendency,handles.prefs_dispersion);
        errorbarplot(g,ct{g},dp{g}-ct{g},dn{g}-ct{g},color{g});
      end
      fprintf(fid,['%% xdata\n']);  fprintf(fid,'%s, ',handles.grouplist{gg});  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT+D\n']);  fprintf(fid,'%g, ',[dp{:}]);  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT-D\n']);  fprintf(fid,'%g, ',[dn{:}]);  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT\n']);  fprintf(fid,'%g, ',[ct{:}]);  fprintf(fid,'\n');
    case 3  % per fly, grouped
      exp_separators=[];  maxy=0;
      table_data{b}=cell(1,length(frames_labelled));
      for g=1:length(frames_labelled)
        tmp=frames_labelled{g}(:,b);
        tmp2=frames_total{g}(:,b);
        cumsum(cellfun(@length,tmp));
        exp_separators=[exp_separators; ans+sum(k)];
        table_data{b}{g}=100.*[tmp{:}]./[tmp2{:}];
        maxy=max([maxy table_data{b}{g}]);
        bar((1:length(table_data{b}{g}))+sum(k),table_data{b}{g},color{g},'barwidth',1,'edgecolor','none');
        k(end+1)=length(table_data{b}{g});
        fprintf(fid,['%% data, %s\n'],handles.grouplist{g});
        fprintf(fid,'%g, ',table_data{b}{g});
        fprintf(fid,'\n');
      end
      l=exp_separators(1:2:(end-1));
      r=exp_separators(2:2:end);
      h=patch(0.5+[l r r l l]',repmat([0 0 maxy*1.05 maxy*1.05 0]',1,floor(length(exp_separators)/2)),...
          [0.95 0.95 0.95]);
      set(h,'edgecolor','none');
      set(gca,'children',flipud(get(gca,'children')));
      k=round(cumsum(k)-k/2);
    case 4  % per fly, stern-style
      m=0;
      table_data{b}=cell(1,length(frames_labelled));
      for g=1:length(frames_labelled)
        table_data{b}{g}=cell(1,length(frames_labelled{g}(:,b)));
        for e=1:length(frames_labelled{g}(:,b))
          table_data{b}{g}{e}=100.*frames_labelled{g}{e,b}./frames_total{g}{e,b};
          [ct,dp,dn]=calculate_ct_d(table_data{b}{g}{e},handles.prefs_centraltendency,handles.prefs_dispersion);
          plot(m,ct,[color{g} 'o']);
          plot([m m],[dp dn],[color{g} '-']);
          plot(m+(1:length(table_data{b}{g}{e})),table_data{b}{g}{e},[color{g} '.']);
          m=m+16+length(table_data{b}{g}{e});
        end
        [ct,dp,dn]=calculate_ct_d([table_data{b}{g}{:}],handles.prefs_centraltendency,handles.prefs_dispersion);
        plot(m,ct,[color{g} 'o'],'markersize',9);
        plot([m m],[dp dn],[color{g} '-'],'linewidth',3);
        m=m+24;
        k(end+1)=24+16*length(table_data{b}{g})+length([table_data{b}{g}{:}]);
        fprintf(fid,['%% data, %s\n'],handles.grouplist{g});
        fprintf(fid,'%g, ',[table_data{b}{g}{:}]);
        fprintf(fid,'\n');
      end
      k=round(cumsum(k)-k/2);
  end
  if(isempty(k))  k=1:length(frames_labelled);  end
  %title(tstr);
  ylabel(ystr);
  set(gca,'xtick',k,'xticklabel',xticklabels);
  axis tight;  vt=axis;
  axisalmosttight;  vat=axis;
  if(handles.behaviorbarchart_perwhat==3)
    axis([vat(1) vat(2) 0 vt(4)]);
  else
    axis([vat(1) vat(2) 0 vat(4)]);
  end
  fprintf(fid,'\n');
end

if((ismember(handles.behaviorbarchart_perwhat,[2 3])) && (individual<4))
  tmp={'behavior'};
  for b=1:length(table_data)
    tmp{4*b+1,1}=handles.behaviorlist{b};
    for g=1:length(table_data{b})
      if(b==1)  tmp{1,g+1}=handles.grouplist{g};  end
      [~,p,~,~]=kstest(table_data{b}{g});
      tmp{4*b+1,g+1}=p;
    end
    if(b==1)  tmp{1,g+2}='(K-S normal)';  end
    k=1;
    for g=1:length(table_data{b})-1
      for g2=(g+1):length(table_data{b})
        if(b==1)
          tmp{2,k+1}=strcat(handles.grouplist{g},'-',handles.grouplist{g2});
          tmp{3,k+1}=strcat(handles.grouplist{g},'-',handles.grouplist{g2});
        end
        tmp{4*b+2,k+1}=ranksum(table_data{b}{g},table_data{b}{g2});
        [~,tmp{4*b+3,k+1}]=ttest2(table_data{b}{g},table_data{b}{g2});
        k=k+1;
      end
    end
    if(b==1)
      tmp{2,k+1}='(ranksum)';
      tmp{3,k+1}='(t-test)';
    end
  end
  fprintf(fid,'%%');  fprintf(fid,'%s, ',tmp{1,1});      fprintf(fid,'\n');
  fprintf(fid,'%%');  fprintf(fid,'%s, ',tmp{1,2:end});  fprintf(fid,'\n');
  fprintf(fid,'%%');  fprintf(fid,'%s, ',tmp{2,2:end});  fprintf(fid,'\n');
  fprintf(fid,'%%');  fprintf(fid,'%s, ',tmp{3,2:end});  fprintf(fid,'\n\n');
  for i=2:((size(tmp,1)+1)/4)
    fprintf(fid,'%s, ',tmp{4*(i-1)+1,1});      fprintf(fid,'\n');
    fprintf(fid,'%g, ',tmp{4*(i-1)+1,2:end});  fprintf(fid,'\n');
    fprintf(fid,'%g, ',tmp{4*(i-1)+2,2:end});  fprintf(fid,'\n');
    fprintf(fid,'%g, ',tmp{4*(i-1)+3,2:end});  fprintf(fid,'\n\n');
  end
  set(handles.Table,'Data',tmp);
  set(handles.Table,'ColumnWidth','auto');
  set(handles.Table,'ColumnName',{''});
else
  set(handles.Table,'Data',[]);
  set(handles.Table,'ColumnName',{});
end

fclose(fid);

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');


% ---
function plot_behavior_timeseries(experiment_value,experiment_list,...
    behavior_value,behavior_logic,behavior_value2,individual,sexdata,...
    style,centraltendency,dispersion,convolutionwidth,color,fid)

cellfun(@(x) size(x{1},2),sexdata,'uniformoutput',false);
behavior_cumulative=zeros(length(experiment_value),max([ans{:}]));
bad=zeros(1,length(experiment_value));
parfor e=1:length(experiment_value)
  behavior_data=load(fullfile(experiment_list{experiment_value(e)},[behavior_value '.mat']));
  behavior_data2=[];
  if(behavior_logic>1)
    behavior_data2=load(fullfile(experiment_list{experiment_value(e)},[behavior_value2 '.mat']));
  end

  if((length(behavior_data.allScores.scores)~=length(sexdata{experiment_value(e)})) || ...
      ((behavior_logic>1) && (length(behavior_data2.allScores.scores)~=length(sexdata{experiment_value(e)}))))
    bad(e)=1;
    continue;
  end

  parfor_tmp=zeros(1,length(behavior_cumulative(e,:)));
  k=0;
  for i=1:length(behavior_data.allScores.t0s)   % individual
    %if((eventdata.Indices(r,2)>4)&&(i~=ii))  continue;  end
    if((individual>3)&&((individual-3)~=i))  continue;  end
    %partition_idx=false(1,length(sexdata{experiment_value0(e)}{i}));
    tmp1=false(1,length(sexdata{experiment_value(e)}{i}));
    for j=1:length(behavior_data.allScores.t0s{i})  % bout
      tmp1((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
          -behavior_data.allScores.tStart(i)+1)=true;
    end
    tmp2=[];
    if(behavior_logic>1)
      tmp2=false(1,length(sexdata{experiment_value(e)}{i}));
      for j=1:length(behavior_data2.allScores.t0s{i})  % bout
        tmp2((behavior_data2.allScores.t0s{i}(j):(behavior_data2.allScores.t1s{i}(j)-1))...
            -behavior_data2.allScores.tStart(i)+1)=true;
      end
    end
    partition_idx=[];
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
    %switch(eventdata.Indices(r,2))
    switch(individual)
      case(2)
        partition_idx = partition_idx & sexdata{e}{i}(1:length(partition_idx));
      case(3)
        partition_idx = partition_idx & (~sexdata{e}{i}(1:length(partition_idx)));
    end
    idx=find(partition_idx);
    parfor_tmp(idx)=parfor_tmp(idx)+1;
    k=k+1;
  end
  behavior_cumulative(e,:)=parfor_tmp./k;
  behavior_cumulative(e,:)=conv(behavior_cumulative(e,:),ones(1,convolutionwidth),'same')...
      ./conv(ones(1,length(behavior_cumulative(e,:))),ones(1,convolutionwidth),'same');
end

tmp=find(bad);
if(length(tmp)>0)
  behavior_cumulative=behavior_cumulative(~logical(bad),:);
  msg{1}=['the following experiments had different numbers of individuals ' ...
      'for the selected behavior and sex feature.  they will be ommitted from the analysis.'];
  msg{2}='';
  for i=1:length(tmp)
    [~,msg{end+1},~]=fileparts(experiment_list{i});
  end
  uiwait(errordlg(msg));
end

plot_it(1:size(behavior_cumulative,2),100.*behavior_cumulative,style,centraltendency,dispersion,color,1,...
    fid,experiment_list(experiment_value));



% --- Executes on button press in BehaviorTimeSeries.
function BehaviorTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals)];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];

gg=1:length(handles.grouplist);
experimentvalue=handles.experimentvalue;
individual=handles.individualvalue;
if(individual>3)
  tmp=find(cumsum_num_indi_per_exp<(individual-3),1,'last');
  gg=find(cumsum_num_exp_per_group<tmp,1,'last');
  experimentvalue{gg}=tmp-cumsum_num_exp_per_group(gg);
  individual=individual-cumsum_num_indi_per_exp(tmp);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='behavior time series';
tstr='';
xstr='time (frames)';
%ylabel([char(strrep(handles.behaviorlist(handles.behaviorvalue),'_','-')) ' (%)']);
ystr=char(strrep(handles.behaviorlist(handles.behaviorvalue),'_','-'));
switch(handles.behaviorlogic)
  case 2
    ystr=[ystr ' AND '];
  case 3
    ystr=[ystr ' AND NOT '];
  case 4
    ystr=[ystr ' OR '];
end
if(handles.behaviorlogic>1)
  ystr=[ystr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
end
ystr=[ystr ' (%)'];

print_csv_help(fid,handles.type,tstr,xstr,ystr);

guidata(figure('ButtonDownFcn',@ButtonDownFcn_Callback),handles);  hold on;

table_data={};  raw_table_data={};
for g=gg
  set(handles.Status,'string',...
      ['Processing ' num2str(length(handles.experimentvalue{g})) ' experiment(s) in group ' handles.grouplist{g}]);
  drawnow;
  if(~isempty(handles.experimentvalue{g}))
    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);
    plot_behavior_timeseries(experimentvalue{g},handles.experimentlist{g},...
        handles.behaviorlist{handles.behaviorvalue},...
        handles.behaviorlogic,handles.behaviorlist{handles.behaviorvalue2},...
        individual,handles.sexdata((cumsum_num_exp_per_group(g)+1):cumsum_num_exp_per_group(g+1)),...
        handles.behaviortimeseries_style,handles.prefs_centraltendency,handles.prefs_dispersion,...
        handles.prefs_convolutionwidth,handles.colors{1,1+mod(g-1,length(handles.colors))},fid);
  end
end

fclose(fid);

xlabel(xstr);
ylabel(ystr);
axis tight;  zoom reset;

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');

%tmp={};
%for r=1:size(eventdata.Indices,1)
%  if((length(experiment_value)>0) && (length(experiment_value2)>0))
%    b=ceil(eventdata.Indices(r,1)/2);
%    if(mod(eventdata.Indices(r,1),2))
%      experiment_value0=experiment_value;
%      experiment_list0=experiment_list;
%      sexdata=handles.sexdata(1:length(experiment_list));
%      color='r';
%    else
%      experiment_value0=experiment_value2;
%      experiment_list0=experiment_list2;
%      sexdata=handles.sexdata((length(experiment_list)+1):end);
%      color='b';
%    end
%  elseif(length(experiment_value)>0)
%    b=eventdata.Indices(r,1);
%    experiment_value0=experiment_value;
%    experiment_list0=experiment_list;
%      sexdata=handles.sexdata(1:length(experiment_list));
%    color='r';
%  else
%    b=eventdata.Indices(r,1);
%    experiment_value0=experiment_value2;
%    experiment_list0=experiment_list2;
%    sexdata=handles.sexdata((length(experiment_list)+1):end);
%    color='b';
%  end
%
%  if(eventdata.Indices(r,2)>4)
%    e=1;  n=0;
%    while((e<=length(experiment_value0)) && (n<eventdata.Indices(r,2)-4))
%      behavior_data=load(fullfile(experiment_list0{experiment_value0(e)},...
%          [behavior_list{b} '.mat']));
%      e=e+1;
%      n=n+length(behavior_data.allScores.t0s);
%    end
%    if(n<eventdata.Indices(r,2)-4)
%      set(handles.Status,'string','Ready.');
%      return;
%    end
%    ee=e-1;
%    ii=eventdata.Indices(r,2)-4-(n-length(behavior_data.allScores.t0s));
%  else
%    ee=1:length(experiment_value0);
%    ii=nan;
%  end
%
%  handles.behaviorvalue=b;
%  set(handles.BehaviorList,'Value',handles.behaviorvalue);
%  handles.individualvalue=eventdata.Indices(r,2)-1;
%  if(((color=='r') && ...
%      (handles.individualvalue>(3+sum(handles.individuals(1:length(handles.experimentlist)))))) || ...
%     ((color=='b') && ...
%      (handles.individualvalue>(3+sum(handles.individuals((end-length(handles.experimentlist2)+1):end))))))
%    return;
%  end
%  if((color=='b') && (handles.individualvalue>3))
%    handles.individualvalue=handles.individualvalue+sum(handles.individuals(1:length(handles.experimentlist)));
%  end
%  set(handles.IndividualList,'Value',handles.individualvalue);
%
%  cellfun(@(x) size(x{1},2),sexdata,'uniformoutput',false);
%  behavior_cumulative=zeros(length(experiment_value0),max([ans{:}]));
%  k=1;
%  parfor e=ee
%    behavior_data=load(fullfile(experiment_list0{experiment_value0(e)},...
%        [behavior_list{b} '.mat']));
%    parfor_tmp=zeros(1,length(behavior_cumulative(e,:)));
%    for i=1:length(behavior_data.allScores.t0s)   % individual
%      if((eventdata.Indices(r,2)>4)&&(i~=ii))  continue;  end
%      partition_idx=false(1,length(sexdata{experiment_value0(e)}{i}));
%      for j=1:length(behavior_data.allScores.t0s{i})  % bout
%        partition_idx((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
%            -behavior_data.allScores.tStart(i)+1)=true;
%      end
%      switch(eventdata.Indices(r,2))
%        case(3)
%          partition_idx = partition_idx & sexdata{e}{i}(1:length(partition_idx));
%        case(4)
%          partition_idx = partition_idx & (~sexdata{e}{i}(1:length(partition_idx)));
%      end
%      idx=find(partition_idx);
%      parfor_tmp(idx)=parfor_tmp(idx)+1;
%      k=k+1;
%    end
%    behavior_cumulative(e,:)=parfor_tmp;
%  end
%  behavior_cumulative=sum(behavior_cumulative,1)./(k-1).*100;
%  behavior_cumulative=conv(behavior_cumulative,...
%      ones(1,handles.prefs_convolutionwidth)./handles.prefs_convolutionwidth,'same');
%
%  plot(behavior_cumulative,[color '-']);
%
%  if(eventdata.Indices(r,2)<4)
%    data=handles.raw_table_data{eventdata.Indices(r,1),eventdata.Indices(r,2)};
%    tmp{end+1,1}=sprintf('%10.3g ',mean(data));
%    tmp{end  ,2}=sprintf('%10.3g ',std(data));
%    tmp{end  ,3}=sprintf('%10.3g ',std(data)./sqrt(length(data)));
%    tmp{end  ,4}=sprintf('%10.3g ',median(data));
%    tmp{end  ,5}=sprintf('%10.3g ',prctile(data,25));
%    tmp{end  ,6}=sprintf('%10.3g ',prctile(data,75));
%    tmp{end  ,7}=sprintf('%10.3g ',prctile(data,5));
%    tmp{end  ,8}=sprintf('%10.3g ',prctile(data,95));
%  end
%end
%
%xlabel('time (frames)');
%ylabel([char(strrep(behavior_list{b},'_','-')) ' (%)']);
%axis tight;
%
%{'Mean' 'Std.Dev.' 'Std.Err.' 'Median' '25%' '75%' '5%' '95%'};
%tmp2(1,:)=sprintf('%10s ',ans{:});
%for i=1:size(tmp,1)
%  tmp2(i+1,:)=[tmp{i,:}];
%end
%v=axis;
%h=text(v(1),v(4),tmp2,'color',[0 0.5 0],'tag','stats','verticalalignment','top','fontname','fixed');
%if(handles.stats)
%  set(h,'visible','on');
%else
%  set(h,'visible','off');
%end
%
%guidata(hObject,handles);


% --- 
function [bout_lengths sex inter_bout_lengths inter_sex]=...
    calculate_boutstats2(behavior_data,behavior_logic,behavior_data2,sexdata)

bout_lengths=cell(1,length(behavior_data.allScores.t0s));
inter_bout_lengths=cell(1,length(behavior_data.allScores.t0s));
sex=cell(1,length(behavior_data.allScores.t0s));
inter_sex=cell(1,length(behavior_data.allScores.t0s));
for i=1:length(behavior_data.allScores.t0s)  % individual
  tmp1=false(1,length(sexdata{i}));
  for j=1:length(behavior_data.allScores.t0s{i})  % bout
    tmp1((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
        -behavior_data.allScores.tStart(i)+1)=true;
  end
  if(behavior_logic>1)
    tmp2=false(1,length(sexdata{i}));
    for j=1:length(behavior_data2.allScores.t0s{i})  % bout
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
  partition_idx=[0 partition_idx 0];
  start=1+find(~partition_idx(1:(end-1)) &  partition_idx(2:end))-1;
  stop =  find( partition_idx(1:(end-1)) & ~partition_idx(2:end))-1;
  if(length(start)>0)
    bout_lengths{i}=stop-start+1;
    sex{i}=zeros(1,length(bout_lengths{i}));
    for j=1:length(bout_lengths{i})
      sex{i}(j)=sum(sexdata{i}(start(j):stop(j))) > (bout_lengths{i}(j)/2);
    end
    if(length(start)>1)
      inter_bout_lengths{i}=start(2:end)-stop(1:(end-1));
      inter_sex{i}=zeros(1,length(inter_bout_lengths{i}));
      for j=1:length(inter_bout_lengths{i})
        inter_sex{i}(j)=sum(sexdata{i}(stop(j):start(j+1))) > (inter_bout_lengths{i}(j)/2);
      end
    end
  end
end


% ---
function [table_data raw_table_data]=calculate_boutstats(experiment_value,experiment_list,...
    behavior_list,behavior_logic,behavior_value2,behavior_list2,sexdata,stat)

collated_data=cell(length(experiment_value),length(behavior_list));
parfor e=1:length(experiment_value)
  behavior_data=load(fullfile(experiment_list{experiment_value(e)},[behavior_list{1} '.mat']));
  if(behavior_logic>1)
    behavior_data2=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list2{behavior_value2} '.mat']));
  else
    behavior_data2=[];
  end
  parfor_tmp=cell(1,length(behavior_list));
  for b=1:length(behavior_list)
    behavior_data=load(fullfile(experiment_list{experiment_value(e)},[behavior_list{b} '.mat']));

    [bout_lengths sex inter_bout_lengths inter_sex]=...
        calculate_boutstats2(behavior_data,behavior_logic,behavior_data2,sexdata{e});
    parfor_tmp{b}={bout_lengths sex inter_bout_lengths inter_sex};
  end
  collated_data(e,:)=parfor_tmp;
end

table_data={};
raw_table_data={};
for b=1:length(behavior_list)
  bout_lengths=[];  sex=[];  n=[];
  inter_bout_lengths=[];  inter_sex=[];  inter_n=[];
  for e=1:length(experiment_value)
    bout_lengths=[bout_lengths collated_data{e,b}{1}{:}];
    sex=[sex collated_data{e,b}{2}{:}];
    n=[n cellfun(@numel,collated_data{e,b}{1})];
    inter_bout_lengths=[inter_bout_lengths collated_data{e,b}{3}{:}];
    inter_sex=[inter_sex collated_data{e,b}{4}{:}];
    inter_n=[inter_n cellfun(@numel,collated_data{e,b}{3})];
  end
  n=[0 cumsum(n)];
  inter_n=[0 cumsum(inter_n)];
  table_data{b,1}=behavior_list{b};
  raw_table_data{b,2}=bout_lengths;
      table_data{b,2}=print_stats(raw_table_data{b,2},stat);
  raw_table_data{b,3}=inter_bout_lengths;
      table_data{b,3}=print_stats(raw_table_data{b,3},stat);
  raw_table_data{b,4}=bout_lengths(sex==1);
      table_data{b,4}=print_stats(raw_table_data{b,4},stat);
  raw_table_data{b,5}=bout_lengths(sex==0);
      table_data{b,5}=print_stats(raw_table_data{b,5},stat);
  raw_table_data{b,6}=inter_bout_lengths(inter_sex==1);
      table_data{b,6}=print_stats(raw_table_data{b,6},stat);
  raw_table_data{b,7}=inter_bout_lengths(inter_sex==0);
      table_data{b,7}=print_stats(raw_table_data{b,7},stat);
  for i=1:(length(n)-1)
    raw_table_data{b,6+2*i}=bout_lengths((n(i)+1):n(i+1));
        table_data{b,6+2*i}=print_stats(raw_table_data{b,7+i},stat);
    raw_table_data{b,7+2*i}=inter_bout_lengths((inter_n(i)+1):inter_n(i+1));
        table_data{b,7+2*i}=print_stats(raw_table_data{b,8+i},stat);
  end
end


% --- Executes on button press in BoutStats.
function BoutStats_Callback(hObject, eventdata, handles)
% hObject    handle to BoutStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_list=get(handles.BehaviorList,'String');
behavior_logic=get(handles.BehaviorLogic,'Value');
behavior_value2=get(handles.BehaviorList2,'Value');
behavior_list2=get(handles.BehaviorList2,'String');
sexdata=handles.sexdata;

if(length(experiment_value)>0)
  [table_data raw_table_data]=calculate_boutstats(experiment_value,experiment_list,...
      behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      sexdata(experiment_value),handles.prefsstat);
end
if(length(experiment_value2)>0)
  [table_data2 raw_table_data2]=calculate_boutstats(experiment_value2,experiment_list2,...
      behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      sexdata(experiment_value2+length(experiment_list)),handles.prefsstat);
end

if((length(experiment_value)>0) && (length(experiment_value2)>0))
  tmp(1:2:2*size(table_data,1),1:size(table_data, 2))=table_data;
  tmp(2:2:2*size(table_data2,1),1:size(table_data2,2))=table_data2;
  raw_tmp(1:2:2*size(raw_table_data,1),1:size(raw_table_data, 2))=raw_table_data;
  raw_tmp(2:2:2*size(raw_table_data2,1),1:size(raw_table_data2,2))=raw_table_data2;
  ii=max(size(table_data,2),size(table_data2,2));
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','on','BackgroundColor',[1 0.8 0.8; 0.8 0.8 1]);
elseif(length(experiment_value)>0)
  tmp=table_data;
  raw_tmp=raw_table_data;
  ii=size(table_data,2);
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','off','BackgroundColor',[1 0.8 0.8]);
else
  tmp=table_data2;
  raw_tmp=raw_table_data2;
  ii=size(table_data2,2);
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','off','BackgroundColor',[0.8 0.8 1]);
end

handles.raw_table_data=raw_tmp;
tmp={'Behavior' 'Bout Length' 'Inter BL' 'Male BL' 'Male IBL' 'Female BL' 'Female IBL'};
cellstr(num2str((1:((ii-7)/2))','Indi #%d BL'))';
[tmp{8:2:(ii-1)}]=deal(ans{:});
cellstr(num2str((1:((ii-7)/2))','Indi #%d IBL'))';
[tmp{9:2:ii}]=deal(ans{:});
%for i=1:(ii-7)/2
%  tmp{2*i+6}=['Indi #' num2str(i) ' BL'];
%  tmp{2*i+7}=['Indi #' num2str(i) ' IBL'];
%end
set(handles.Table,'ColumnName',tmp);
set(handles.Table,'ColumnWidth',{150 75 75});

handles.table='bout_stats';
guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');


% --- 
function [during not_during raw_during raw_not_during]=...
    calculate_social(behavior_data,behavior_logic,behavior_data2,feature_data)

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

during={};  not_during={};
raw_during={};  raw_not_during={};
for i=1:length(behavior_data.allScores.t0s)  % individual
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
      partition_idx=tmp1;
    case(2)
      partition_idx=tmp1 & tmp2;
    case(3)
      partition_idx=tmp1 & ~tmp2;
    case(4)
      partition_idx=tmp1 | tmp2;
  end
  if length(partition_idx)==(length(feature_data.data{i}+1))  % n-1 features
    partition_idx=partition_idx(1:end-1);
  end
  partition_idx=[0 partition_idx 0];
  start=1+find(~partition_idx(1:(end-1)) &  partition_idx(2:end))-1;
  stop =  find( partition_idx(1:(end-1)) & ~partition_idx(2:end))-1;
  during{i}={};  not_during{i}={};
  raw_during{i}={};  raw_not_during{i}={};
  if(length(start)>0)
    for j=1:length(start)
      raw_during{i}{j}=feature_data.data{i}(start(j):stop(j));
      [m,f,c]=mode(raw_during{i}{j});
      during{i}{j}={c{1} f./length(raw_during{i}{j})*100};
    end
    if(length(start)>1)
      for j=1:(length(start)-1)
        raw_not_during{i}{j}=feature_data.data{i}(stop(j):start(j+1));
        [m,f,c]=mode(raw_not_during{i}{j});
        not_during{i}{j}={c{1} f./length(raw_not_during{i}{j})*100};
      end
    end
  end
end



% ---
function ret_val=print_stats2(arg)

sprintf('%d,',arg{1});
ret_val=[ans(1:(end-1)) ' (' num2str(arg{2},3) '%)'];



% ---
function [table_data raw_table_data]=plot_social(experiment_value,experiment_list,...
    behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
    feature_value,feature_list,color)

during={};  not_during={};
parfor e=1:length(experiment_value)
%for e=1:length(experiment_value)
  behavior_data=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list{behavior_value} '.mat']));
  if(behavior_logic>1)
    behavior_data2=load(fullfile(experiment_list{experiment_value(e)},...
        [behavior_list2{behavior_value2} '.mat']));
  else
    behavior_data2=[];
  end
  feature_data=load(fullfile(experiment_list{experiment_value(e)},'perframe',...
        [feature_list{feature_value} '.mat']));

  [during{e} not_during{e} raw_during{e} raw_not_during{e}]=...
      calculate_social(behavior_data,behavior_logic,behavior_data2,feature_data);
end

k=1;
table_data={};
raw_table_data={};
for e=1:length(experiment_value)
  for i=1:length(during{e})
    table_data{k,1}=['Exp #' color num2str(experiment_value(e)) ', Indi #' num2str(i)];
    if ~isempty(during{e}{i})
      cellfun(@(x) x{1},during{e}{i},'uniformoutput',false);
      raw_table_data{k,2}=cat(1,ans{:});
      [m,f,c]=mode(raw_table_data{k,2});
      sprintf('%d,',c{1});
      table_data{k,2}=[ans(1:(end-1)) ' (' num2str(f/length(raw_table_data{k,2})*100,3) '%)'];
    end
    if ~isempty(not_during{e}{i})
      cellfun(@(x) x{1},not_during{e}{i},'uniformoutput',false);
      raw_table_data{k,3}=cat(1,ans{:});
      [m,f,c]=mode(raw_table_data{k,3});
      sprintf('%d,',c{1});
      table_data{k,3}=[ans(1:(end-1)) ' (' num2str(f/length(raw_table_data{k,3})*100,3) '%)'];
    end
    if ~isempty(during{e}{i})
      tmp=cellfun(@print_stats2,during{e}{i},'uniformoutput',false);
      [table_data{k,4:2:(2+2*length(during{e}{i}))}]=deal(tmp{:});
      [raw_table_data{k,4:2:(2+2*length(during{e}{i}))}]=deal(raw_during{e}{i}{:});
    end
    if ~isempty(not_during{e}{i})
      tmp=cellfun(@print_stats2,not_during{e}{i},'uniformoutput',false);
      [table_data{k,5:2:(3+2*length(not_during{e}{i}))}]=deal(tmp{:});
      [raw_table_data{k,5:2:(3+2*length(not_during{e}{i}))}]=deal(raw_not_during{e}{i}{:});
    end
    k=k+1;
  end
end


% --- Executes on button press in SocialStats.
function SocialStats_Callback(hObject, eventdata, handles)
% hObject    handle to SocialStats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.Status,'string','Thinking...','foregroundcolor','b');  drawnow;
set(handles.figure1,'pointer','watch');

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

if(~strncmp(feature_list(feature_value),'closestfly_',11))
  feature_value=find(cellfun(@(x) strncmp(x,'closestfly_',11),feature_list),1);
  if(isempty(feature_value))
    error('no closestfly feature found');
  end
  set(handles.FeatureList,'Value',feature_value);
end

%axes(handles.Axes);  cla;  hold on;
figure;  hold on;

if(length(experiment_value)>0)
  [table_data raw_table_data]=plot_social(experiment_value,experiment_list,...
      behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      feature_value,feature_list,'R');
end
if(length(experiment_value2)>0)
  [table_data2 raw_table_data2]=plot_social(experiment_value2,experiment_list2,...
      behavior_value,behavior_list,behavior_logic,behavior_value2,behavior_list2,...
      feature_value,feature_list,'B');
end

if((length(experiment_value)>0) && (length(experiment_value2)>0))
  tmp(1:2:2*size(table_data,1),1:size(table_data, 2))=table_data;
  tmp(2:2:2*size(table_data2,1),1:size(table_data2,2))=table_data2;
  raw_tmp(1:2:2*size(raw_table_data,1),1:size(raw_table_data, 2))=raw_table_data;
  raw_tmp(2:2:2*size(raw_table_data2,1),1:size(raw_table_data2,2))=raw_table_data2;
  ii=max(size(table_data,2),size(table_data2,2));
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','on','BackgroundColor',[1 0.8 0.8; 0.8 0.8 1]);
elseif(length(experiment_value)>0)
  tmp=table_data;
  raw_tmp=raw_table_data;
  ii=size(table_data,2);
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','off','BackgroundColor',[1 0.8 0.8]);
else
  tmp=table_data2;
  raw_tmp=raw_table_data2;
  ii=size(table_data2,2);
  set(handles.Table,'Data',tmp);
  set(handles.Table,'RowStriping','off','BackgroundColor',[0.8 0.8 1]);
end

handles.raw_table_data=raw_tmp;
tmp={'Individual' 'Bout Mode' 'Inter Bout Mode'};
cellstr(num2str((1:(((ii-3)/2)+1))','B #%d'))';
[tmp{4:2:ii}]=deal(ans{:});
cellstr(num2str((1: ((ii-3)/2)   )','IB #%d'))';
[tmp{5:2:ii}]=deal(ans{:});
%for i=1:(ii-3)/2
%  tmp{2*i+2}=['B #' num2str(i)];
%  tmp{2*i+3}=['IB #' num2str(i)];
%end
set(handles.Table,'ColumnName',tmp);
set(handles.Table,'ColumnWidth',{100 75 100});

handles.table='social_stats';
guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');  drawnow;
set(handles.figure1,'pointer','arrow');


% --- 
function CellSelectionCallback(hObject,eventdata)

if(size(eventdata.Indices,1)==0)  return;  end

handles=guidata(hObject);

if(strcmp(handles.table,'bout_stats'))
  if(eventdata.Indices(end,2)==1)  return;  end

  %axes(handles.Axes);  cla;  hold on;
  figure;  hold on;
  tmp={};
  for r=1:size(eventdata.Indices,1)
    experiment_value=get(handles.ExperimentList,'Value');
    experiment_value2=get(handles.ExperimentList2,'Value');
    if((length(experiment_value)>0) && (length(experiment_value2)>0))
      handles.behaviorvalue=ceil(eventdata.Indices(r,1)/2);
      if(mod(eventdata.Indices(r,1),2))
        color='r';
      else
        color='b';
      end
    elseif(length(experiment_value)>0)
      handles.behaviorvalue=eventdata.Indices(r,1);
      color='r';
    else
      handles.behaviorvalue=eventdata.Indices(r,1);
      color='b';
    end
    set(handles.BehaviorList,'Value',handles.behaviorvalue);
    handles.individualvalue=floor(eventdata.Indices(r,2)/2);
    if(((color=='r') && ...
        (handles.individualvalue>(3+sum(handles.individuals(1:length(handles.experimentlist)))))) || ...
       ((color=='b') && ...
        (handles.individualvalue>(3+sum(handles.individuals((end-length(handles.experimentlist2)+1):end))))))
      return;
    end
    if((color=='b') && (handles.individualvalue>3))
      handles.individualvalue=handles.individualvalue+sum(handles.individuals(1:length(handles.experimentlist)));
    end
    set(handles.IndividualList,'Value',handles.individualvalue);

    data=handles.raw_table_data{eventdata.Indices(r,1),eventdata.Indices(r,2)};
    [n,x]=hist(data,100);
    plot(x,n./sum(n),[color '-']);

    tmp{r,1}=sprintf('%10.3g ',mean(data));
    tmp{r,2}=sprintf('%10.3g ',std(data));
    tmp{r,3}=sprintf('%10.3g ',std(data)./sqrt(length(data)));
    tmp{r,4}=sprintf('%10.3g ',median(data));
    tmp{r,5}=sprintf('%10.3g ',prctile(data,25));
    tmp{r,6}=sprintf('%10.3g ',prctile(data,75));
    tmp{r,7}=sprintf('%10.3g ',prctile(data,5));
    tmp{r,8}=sprintf('%10.3g ',prctile(data,95));
  end

  xlabel('length (frames)');
  axis tight;

  {'Mean' 'Std.Dev.' 'Std.Err.' 'Median' '25%' '75%' '5%' '95%'};
  tmp2(1,:)=sprintf('%10s ',ans{:});
  for i=1:size(tmp,1)
    tmp2(i+1,:)=[tmp{i,:}];
  end
  v=axis;
  h=text(v(1),v(4),tmp2,'color',[0 0.5 0],'tag','stats','verticalalignment','top','fontname','fixed');
  if(handles.stats)
    set(h,'visible','on');
  else
    set(h,'visible','off');
  end
  
elseif(strcmp(handles.table,'social_stats'))
  if(eventdata.Indices(end,2)==1)  return;  end

  %axes(handles.Axes);  cla;  hold on;
  figure;  hold on;
  tmp={};
  %for r=1:size(eventdata.Indices,1)
    experiment_value=get(handles.ExperimentList,'Value');
    experiment_value2=get(handles.ExperimentList2,'Value');
    if((length(experiment_value)>0) && (length(experiment_value2)>0))
      if(mod(eventdata.Indices(end,1),2))
        color='r';
      else
        color='b';
      end
    elseif(length(experiment_value)>0)
      color='r';
    else
      color='b';
    end
    %handles.individualvalue=floor(eventdata.Indices(r,2)/2);
    %if(((color=='r') && ...
    %    (handles.individualvalue>(3+sum(handles.individuals(1:length(handles.experimentlist)))))) || ...
    %   ((color=='b') && ...
    %    (handles.individualvalue>(3+sum(handles.individuals((end-length(handles.experimentlist2)+1):end))))))
    %  return;
    %end
    %if((color=='b') && (handles.individualvalue>3))
    %  handles.individualvalue=handles.individualvalue+sum(handles.individuals(1:length(handles.experimentlist)));
    %end
    %set(handles.IndividualList,'Value',handles.individualvalue);

    data=handles.raw_table_data{eventdata.Indices(end,1),eventdata.Indices(end,2)};
    tmp=unique(handles.raw_table_data{eventdata.Indices(end,1),2});
    hist(data,min(tmp):max(tmp));
    %plot(x,n./sum(n),[color '-']);

    %tmp{r,1}=sprintf('%10.3g ',mean(data));
    %tmp{r,2}=sprintf('%10.3g ',std(data));
    %tmp{r,3}=sprintf('%10.3g ',std(data)./sqrt(length(data)));
    %tmp{r,4}=sprintf('%10.3g ',median(data));
    %tmp{r,5}=sprintf('%10.3g ',prctile(data,25));
    %tmp{r,6}=sprintf('%10.3g ',prctile(data,75));
    %tmp{r,7}=sprintf('%10.3g ',prctile(data,5));
    %tmp{r,8}=sprintf('%10.3g ',prctile(data,95));
  %end

  %{'Mean' 'Std.Dev.' 'Std.Err.' 'Median' '25%' '75%' '5%' '95%'};
  %tmp2(1,:)=sprintf('%10s ',ans{:});
  %for i=1:size(tmp,1)
  %  tmp2(i+1,:)=[tmp{i,:}];
  %end
  %v=axis;
  %h=text(v(1),v(4),tmp2,'color',[0 0.5 0],'tag','stats','verticalalignment','top','fontname','fixed');
  %if(handles.stats)
  %  set(h,'visible','on');
  %else
  %  set(h,'visible','off');
  %end
  
  xlabel('closest fly (#)');
  axis tight;

elseif(strcmp(handles.table,'timeseries') || strcmp(handles.table,'histogram'))
  handles.behaviorvalue=handles.table_data(eventdata.Indices(end,1),3);
  set(handles.BehaviorList,'Value',handles.behaviorvalue);
  set(handles.BehaviorLogic,'Value',1);
  set(handles.BehaviorList2,'enable','off');
  handles.featurevalue=handles.table_data(eventdata.Indices(end,1),4);
  set(handles.FeatureList,'Value',handles.featurevalue);
  handles.individual=1;
  set(handles.IndividualList,'Value',handles.individual);

  if(strcmp(handles.table,'timeseries'))
    handles.timeseries_timing=handles.table_data(eventdata.Indices(end,1),3)+1;
    FeatureTimeSeries_Callback(hObject, eventdata, handles);

  elseif(strcmp(handles.table,'histogram'))
    handles.featurehistogram_notduring=isnan(handles.table_data(eventdata.Indices(end,1),2));
    menu_featurehistogram_notduring_set(handles.featurehistogram_notduring);
    FeatureHistogram_Callback(hObject, eventdata, handles);
  end
end

guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ---
function menu_featuretimeseries_subtractmean_set(arg)

handles=guidata(gcf);
if(arg)
  set(handles.MenuFeatureTimeSeriesSubtractMean,'Checked','on');
else
  set(handles.MenuFeatureTimeSeriesSubtractMean,'Checked','off');
end


% --------------------------------------------------------------------
function MenuFeatureTimeSeriesSubtractMean_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureTimeSeriesSubtractMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featuretimeseries_subtractmean=~handles.featuretimeseries_subtractmean;
menu_featuretimeseries_subtractmean_set(handles.featuretimeseries_subtractmean);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureTimeSeriesRadius_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureTimeSeriesRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.timeseries_windowradius=str2num(char(inputdlg({'Window radius:'},'',1,...
    {num2str(handles.timeseries_windowradius)})));
guidata(hObject,handles);


% --------------------------------------------------------------------
%function MenuTimeSeriesTight_Callback(hObject, eventdata, handles)
% hObject    handle to MenuTimeSeriesTight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%handles.timeseries_tight=~handles.timeseries_tight;
%if(handles.timeseries_tight)
%  set(handles.MenuTimeSeriesTight,'Checked','on');
%else
%  set(handles.MenuTimeSeriesTight,'Checked','off');
%end
%guidata(hObject,handles);


% ---
function menu_featuretimeseries_timing_set(arg)

handles=guidata(gcf);
set(handles.MenuFeatureTimeSeriesEntireRecording,'Checked','off');
set(handles.MenuFeatureTimeSeriesOnsetTriggered,'Checked','off');
set(handles.MenuFeatureTimeSeriesOffsetTriggered,'Checked','off');
switch(arg)
  case(1), set(handles.MenuFeatureTimeSeriesEntireRecording,'Checked','on');
  case(2), set(handles.MenuFeatureTimeSeriesOnsetTriggered,'Checked','on');
  case(3), set(handles.MenuFeatureTimeSeriesOffsetTriggered,'Checked','on');
end


% --------------------------------------------------------------------
function MenuFeatureTimeSeriesEntireRecording_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureTimeSeriesEntireRecording (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featuretimeseries_timing=1;
menu_featuretimeseries_timing_set(handles.featuretimeseries_timing);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureTimeSeriesOnsetTriggered_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureTimeSeriesOnsetTriggered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featuretimeseries_timing=2;
menu_featuretimeseries_timing_set(handles.featuretimeseries_timing);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureTimeSeriesOffsetTriggered_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureTimeSeriesOffsetTriggered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featuretimeseries_timing=3;
menu_featuretimeseries_timing_set(handles.featuretimeseries_timing);
guidata(hObject,handles);


% ---
function menu_featuretimeseries_style_set(arg)

handles=guidata(gcf);
set(handles.MenuFeatureTimeSeriesCentralTendency,'Checked','off');
set(handles.MenuFeatureTimeSeriesCentralTendencyDispersion,'Checked','off');
set(handles.MenuFeatureTimeSeriesOverlayedPerExpMeans,'Checked','off');
switch(arg)
  case(1), set(handles.MenuFeatureTimeSeriesCentralTendency,'Checked','on');
  case(2), set(handles.MenuFeatureTimeSeriesCentralTendencyDispersion,'Checked','on');
  case(3), set(handles.MenuFeatureTimeSeriesOverlayedPerExpMeans,'Checked','on');
end


% --------------------------------------------------------------------
function MenuFeatureTimeSeriesCentralTendency_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureTimeSeriesCentralTendency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featuretimeseries_style=1;
menu_featuretimeseries_style_set(handles.featuretimeseries_style);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureTimeSeriesCentralTendencyDispersion_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureTimeSeriesCentralTendencyDispersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featuretimeseries_style=2;
menu_featuretimeseries_style_set(handles.featuretimeseries_style);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureTimeSeriesOverlayedPerExpMeans_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureTimeSeriesOverlayedPerExpMeans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featuretimeseries_style=3;
menu_featuretimeseries_style_set(handles.featuretimeseries_style);
guidata(hObject,handles);


% ---
function menu_featurehistogram_perwhat_set(arg)

handles=guidata(gcf);

set(handles.MenuFeatureHistogramPerFrame,'Checked','off');
set(handles.MenuFeatureHistogramMeanPerBout,'Checked','off');
set(handles.MenuFeatureHistogramMedianPerBout,'Checked','off');
set(handles.MenuFeatureHistogramMaxPerBout,'Checked','off');
set(handles.MenuFeatureHistogramMinPerBout,'Checked','off');
set(handles.MenuFeatureHistogramStdPerBout,'Checked','off');
switch(arg)
  case(1), set(handles.MenuFeatureHistogramPerFrame,'Checked','on');
  case(2), set(handles.MenuFeatureHistogramMeanPerBout,'Checked','on');
  case(3), set(handles.MenuFeatureHistogramMedianPerBout,'Checked','on');
  case(4), set(handles.MenuFeatureHistogramMaxPerBout,'Checked','on');
  case(5), set(handles.MenuFeatureHistogramMinPerBout,'Checked','on');
  case(6), set(handles.MenuFeatureHistogramStdPerBout,'Checked','on');
end


% --------------------------------------------------------------------
function MenuFeatureHistogramPerFrame_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramPerFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_perwhat=1;
menu_featurehistogram_perwhat_set(handles.featurehistogram_perwhat);
handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureHistogramMeanPerBout_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramMeanPerBout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_perwhat=2;
menu_featurehistogram_perwhat_set(handles.featurehistogram_perwhat);
handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureHistogramMedianPerBout_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramMedianPerBout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_perwhat=3;
menu_featurehistogram_perwhat_set(handles.featurehistogram_perwhat);
handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureHistogramMaxPerBout_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramMaxPerBout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_perwhat=4;
menu_featurehistogram_perwhat_set(handles.featurehistogram_perwhat);
handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureHistogramMinPerBout_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramMinPerBout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_perwhat=5;
menu_featurehistogram_perwhat_set(handles.featurehistogram_perwhat);
handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureHistogramStdPerBout_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramStdPerBout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_perwhat=6;
menu_featurehistogram_perwhat_set(handles.featurehistogram_perwhat);
handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% ---
function menu_featurehistogram_style_set(arg)

handles=guidata(gcf);

set(handles.MenuFeatureHistogramCentralTendency,'Checked','off');
set(handles.MenuFeatureHistogramCentralTendencyDispersion,'Checked','off');
set(handles.MenuFeatureHistogramOverlayedPerExpMeans,'Checked','off');
switch(arg)
  case(1), set(handles.MenuFeatureHistogramCentralTendency,'Checked','on');
  case(2), set(handles.MenuFeatureHistogramCentralTendencyDispersion,'Checked','on');
  case(3), set(handles.MenuFeatureHistogramOverlayedPerExpMeans,'Checked','on');
end


% --------------------------------------------------------------------
function MenuFeatureHistogramCentralTendency_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramCentralTendency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_style=1;
menu_featurehistogram_style_set(handles.featurehistogram_style);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureHistogramCentralTendencyDispersion_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramCentralTendencyDispersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_style=2;
menu_featurehistogram_style_set(handles.featurehistogram_style);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureHistogramOverlayedPerExpMeans_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramOverlayedPerExpMeans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_style=3;
menu_featurehistogram_style_set(handles.featurehistogram_style);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in Stats.
%function Stats_Callback(hObject, eventdata, handles)
% hObject    handle to Stats (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%
%handles.stats=~handles.stats;
%h=findobj('tag','stats');
%if(handles.stats)
%  set(h,'visible','on');
%  set(handles.Stats,'backgroundcolor',0.4*[1 1 1]);
%else
%  set(h,'visible','off');
%  set(handles.Stats,'backgroundcolor',get(gcf,'color'));
%end
%guidata(hObject, handles);


% --- Executes on button press in CloseAll.
function CloseAll_Callback(hObject, eventdata, handles)
% hObject    handle to CloseAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

questdlg('Close all figures?','','Yes','No','No');
if(strcmp(ans,'No'))  return;  end
hh=findobj('type','figure');
for h=1:length(hh)
  if(hh(h)~=gcf)  close(hh(h));  end
end


%[file,path]=uiputfile('*.txt','Select file to save table to');
%if (file==0)  return;  end
%fid=fopen(fullfile(path,file),'w');
%tmp=get(handles.Table,'columnname');
%tmp=transpose(tmp);
%tmp2=get(handles.Table,'data');
%tmp2(2:end+1,:)=tmp2;
%tmp2(1,:)=tmp;
%tmp2=cellfun(@char,tmp2,'uniformoutput',false);
%for i=1:size(tmp2,1)
%  for j=1:size(tmp2,2)
%    fprintf(fid,'%s\t',char(tmp2(i,j)));
%  end
%  fprintf(fid,'\r\n');
%end
%fclose(fid);


% --- Executes on button press in ExportGraph.
%function ExportGraph_Callback(hObject, eventdata, handles)
% hObject    handle to ExportGraph (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ---
function menu_behaviorbarchart_perwhat_set(arg)

handles=guidata(gcf);

set(handles.MenuBehaviorBarChartPerGroup,'Checked','off');
set(handles.MenuBehaviorBarChartPerExperimentCTD,'Checked','off');
set(handles.MenuBehaviorBarChartPerFlyGrouped,'Checked','off');
set(handles.MenuBehaviorBarChartPerFlyScatter,'Checked','off');
switch(arg)
  case(1), set(handles.MenuBehaviorBarChartPerGroup,'Checked','on');
  case(2), set(handles.MenuBehaviorBarChartPerExperimentCTD,'Checked','on');
  case(3), set(handles.MenuBehaviorBarChartPerFlyGrouped,'Checked','on');
  case(4), set(handles.MenuBehaviorBarChartPerFlyScatter,'Checked','on');
end


% --------------------------------------------------------------------
function MenuBehaviorBarChartPerGroup_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBehaviorBarChartPerGroup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.behaviorbarchart_perwhat=1;
menu_behaviorbarchart_perwhat_set(handles.behaviorbarchart_perwhat);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuBehaviorBarChartPerExperimentCTD_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBehaviorBarChartPerExperimentCTD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.behaviorbarchart_perwhat=2;
menu_behaviorbarchart_perwhat_set(handles.behaviorbarchart_perwhat);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuBehaviorBarChartPerFlyGrouped_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBehaviorBarChartPerFlyGrouped (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.behaviorbarchart_perwhat=3;
menu_behaviorbarchart_perwhat_set(handles.behaviorbarchart_perwhat);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuBehaviorBarChartPerFlyScatter_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBehaviorBarChartPerFlyScatter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.behaviorbarchart_perwhat=4;
menu_behaviorbarchart_perwhat_set(handles.behaviorbarchart_perwhat);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuBehaviorBarChart_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBehaviorBarChart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ---
function menu_interestingfeaturehistograms_omitnan_set(arg)

handles=guidata(gcf);
if(arg)
  set(handles.MenuInterestingFeatureHistogramsOmitNaN,'Checked','on');
else
  set(handles.MenuInterestingFeatureHistogramsOmitNaN,'Checked','off');
end


% --------------------------------------------------------------------
function MenuInterestingFeatureHistogramsOmitNaN_Callback(hObject, eventdata, handles)
% hObject    handle to MenuInterestingFeatureHistogramsOmitNaN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.interestingfeaturehistograms_omitnan=~handles.interestingfeaturehistograms_omitnan;
menu_interestingfeaturehistograms_omitnan_set(handles.interestingfeaturehistograms_omitnan);
handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% ---
function menu_interestingfeaturehistograms_omitinf_set(arg)

handles=guidata(gcf);
if(arg)
  set(handles.MenuInterestingFeatureHistogramsOmitInf,'Checked','on');
else
  set(handles.MenuInterestingFeatureHistogramsOmitInf,'Checked','off');
end


% --------------------------------------------------------------------
function MenuInterestingFeatureHistogramsOmitInf_Callback(hObject, eventdata, handles)
% hObject    handle to MenuInterestingFeatureHistogramsOmitInf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.interestingfeaturehistograms_omitinf=~handles.interestingfeaturehistograms_omitinf;
menu_interestingfeaturehistograms_omitinf_set(handles.interestingfeaturehistograms_omitinf);
handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuInterestingFeatureHistograms_Callback(hObject, eventdata, handles)
% hObject    handle to MenuInterestingFeatureHistograms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ---
function menu_featurehistogram_logbinsize_set(arg)

handles=guidata(gcf);
if(arg)
  set(handles.MenuFeatureHistogramLogBinSize,'Checked','on');
else
  set(handles.MenuFeatureHistogramLogBinSize,'Checked','off');
end


% --------------------------------------------------------------------
function MenuFeatureHistogramLogBinSize_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramLogBinSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_logbinsize=~handles.featurehistogram_logbinsize;
menu_featurehistogram_logbinsize_set(handles.featurehistogram_logbinsize);
guidata(hObject,handles);


% ---
function menu_featurehistogram_notduring_set(arg)

handles=guidata(gcf);
if(arg)
  set(handles.MenuFeatureHistogramNotDuring,'Checked','on');
else
  set(handles.MenuFeatureHistogramNotDuring,'Checked','off');
end


% --------------------------------------------------------------------
function MenuFeatureHistogramNotDuring_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramNotDuring (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_notduring=~handles.featurehistogram_notduring;
menu_featurehistogram_notduring_set(handles.featurehistogram_notduring);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFeatureHistogramNbins_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramNbins (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featurehistogram_nbins=str2num(char(inputdlg({'Number of bins:'},'',1,...
    {num2str(handles.featurehistogram_nbins)})));
guidata(hObject,handles);


% ---
function menu_behaviortimeseries_style_set(arg)

handles=guidata(gcf);

set(handles.MenuBehaviorTimeSeriesCentralTendency,'Checked','off');
set(handles.MenuBehaviorTimeSeriesCentralTendencyDispersion,'Checked','off');
set(handles.MenuBehaviorTimeSeriesOverlayedPerExpMeans,'Checked','off');
switch(arg)
  case(1), set(handles.MenuBehaviorTimeSeriesCentralTendency,'Checked','on');
  case(2), set(handles.MenuBehaviorTimeSeriesCentralTendencyDispersion,'Checked','on');
  case(3), set(handles.MenuBehaviorTimeSeriesOverlayedPerExpMeans,'Checked','on');
end


% --------------------------------------------------------------------
function MenuBehaviorTimeSeriesCentralTendency_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBehaviorTimeSeriesCentralTendency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.behaviortimeseries_style=1;
menu_behaviortimeseries_style_set(handles.behaviortimeseries_style);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuBehaviorTimeSeriesCentralTendencyDispersion_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBehaviorTimeSeriesCentralTendencyDispersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.behaviortimeseries_style=2;
menu_behaviortimeseries_style_set(handles.behaviortimeseries_style);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuBehaviorTimeSeriesOverlayedPerExpMeans_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBehaviorTimeSeriesOverlayedPerExpMeans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.behaviortimeseries_style=3;
menu_behaviortimeseries_style_set(handles.behaviortimeseries_style);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuBehaviorTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to MenuBehaviorTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ---
function [h,h2]=errorbarplot(x,b,dp,dn,color)

h = bar(x,b,'grouped');
set(h,'facecolor',color);
fxdata = get(h, 'XData');
fydata = get(h ,'YData');
if(~iscell(fxdata)) xdata{1,1} = fxdata; else xdata = fxdata; end
if(~iscell(fydata)) ydata{1,1} = fydata; else ydata = fydata; end

% Determine number of bars
sizz = size(b);
nb = sizz(1)*sizz(2);
xb = [];
yb = [];
for i = 1:length(xdata),
    xb = [xb xdata{i,1}];
    yb = [yb ydata{i,1}];
end;

%% To find the center of each bar - need to look at the output vectors xb, yb
%% find where yb is non-zero - for each bar there is a pair of non-zero yb values.
%% The center of these values is the middle of the bar
%
%nz = find(yb);
%for i = 1:nb,
%    center(i) = (xb(nz(i*2))-xb(nz((i*2)-1)))/2 + xb(nz((i*2)-1));
%end;

% To place the error bars - use the following:

%errdata=[.1 .2; .3 .4; .5 .6];

hold on;
h2=errorbar(xb, b, dp, dn);
set(h2(1),'linewidth',1);            % This changes the thickness of the errorbars
set(h2(1),'color','k');              % This changes the color of the errorbars
set(h2(1),'linestyle','none');       % This removes the connecting


% --------------------------------------------------------------------
function MenuPrefs_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ---
function menu_prefscentraltendency_set(arg)

handles=guidata(gcf);

set(handles.MenuPrefsCentralTendencyMean,'Checked','off');
set(handles.MenuPrefsCentralTendencyMedian,'Checked','off');
set(handles.MenuPrefsCentralTendencyMode,'Checked','off');
switch(arg)
  case(1), set(handles.MenuPrefsCentralTendencyMean,'Checked','on');
  case(2), set(handles.MenuPrefsCentralTendencyMedian,'Checked','on');
  case(3), set(handles.MenuPrefsCentralTendencyMode,'Checked','on');
end


% --------------------------------------------------------------------
function MenuPrefsCentralTendency_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsCentralTendency (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuPrefsCentralTendencyMean_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsCentralTendencyMean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_centraltendency=1;
menu_prefscentraltendency_set(handles.prefs_centraltendency);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuPrefsCentralTendencyMedian_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsCentralTendencyMedian (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_centraltendency=2;
menu_prefscentraltendency_set(handles.prefs_centraltendency);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuPrefsCentralTendencyMode_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsCentralTendencyMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_centraltendency=3;
menu_prefscentraltendency_set(handles.prefs_centraltendency);
guidata(hObject,handles);


% ---
function menu_prefsdispersion_set(arg)

handles=guidata(gcf);

set(handles.MenuPrefsDispersionStdDev,'Checked','off');
set(handles.MenuPrefsDispersionStdErr,'Checked','off');
set(handles.MenuPrefsDispersionPrctile5,'Checked','off');
set(handles.MenuPrefsDispersionPrctile25,'Checked','off');
switch(arg)
  case(1), set(handles.MenuPrefsDispersionStdDev,'Checked','on');
  case(2), set(handles.MenuPrefsDispersionStdErr,'Checked','on');
  case(3), set(handles.MenuPrefsDispersionPrctile5,'Checked','on');
  case(4), set(handles.MenuPrefsDispersionPrctile25,'Checked','on');
end


% --------------------------------------------------------------------
function MenuPrefsDispersion_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsDispersion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function MenuPrefsDispersionStdDev_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsDispersionStdDev (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_dispersion=1;
menu_prefsdispersion_set(handles.prefs_dispersion);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuPrefsDispersionStdErr_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsDispersionStdErr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_dispersion=2;
menu_prefsdispersion_set(handles.prefs_dispersion);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuPrefsDispersionPrctile5_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsDispersion5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_dispersion=3;
menu_prefsdispersion_set(handles.prefs_dispersion);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuPrefsDispersionPrctile25_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsDispersion25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_dispersion=4;
menu_prefsdispersion_set(handles.prefs_dispersion);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuPrefsConvolutionWidth_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsConvolutionWidth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_convolutionwidth=str2num(char(inputdlg({'Convolution width:'},'',1,...
    {num2str(handles.prefs_convolutionwidth)})));
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuFileReset_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileReset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

questdlg('Reset configuration to default?','','Yes','No','No');
if(strcmp(ans,'No'))  return;  end
handles=initialize(handles);
update_figure(handles);
set(handles.figure1,'pointer','arrow');
guidata(hObject, handles);


% --------------------------------------------------------------------
function MenuFileLoad_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path]=uigetfile('*.mat','Select configuration file to open');
if(isnumeric(file) && isnumeric(path) && (file==0) && (path==0))  return;  end
handles=load_configuration_file(fullfile(path,file),hObject,eventdata,handles);
update_figure(handles);
guidata(hObject, handles);


% --------------------------------------------------------------------
function MenuFileSave_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path]=uiputfile('*.mat','Select file to save configuration to');
if(isnumeric(file) && isnumeric(path) && (file==0) && (path==0))  return;  end
save(fullfile(path,file),'handles');


% --------------------------------------------------------------------
function MenuFile_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ---
function ButtonDownFcn_Callback(src,evt)

handles=guidata(src);

tmp={};

switch handles.type
  case 'feature histogram'
    tmp{end+1}=['perwhat=' num2str(handles.featurehistogram_perwhat)];
    tmp{end+1}=['logbinsize=' num2str(handles.featurehistogram_logbinsize)];
    tmp{end+1}=['notduring=' num2str(handles.featurehistogram_notduring)];
    tmp{end+1}=['nbins=' num2str(handles.featurehistogram_nbins)];
  case 'feature time series'
    tmp{end+1}=['timing=' num2str(handles.featuretimeseries_timing)];
    tmp{end+1}=['style=' num2str(handles.featuretimeseries_style)];
    tmp{end+1}=['subtractmean=' num2str(handles.featuretimeseries_subtractmean)];
end
tmp{end+1}='';

tmp{end+1}=['CT=' num2str(handles.prefs_centraltendency)];
tmp{end+1}=['D='  num2str(handles.prefs_dispersion)];
tmp{end+1}=['conv. width='  num2str(handles.prefs_convolutionwidth)];
tmp{end+1}='';

for g=1:length(handles.grouplist)
  if(length(handles.experimentvalue{g})==0)  continue;  end
  tmp{end+1}=['group ' handles.grouplist{g}];
  for e=1:length(handles.experimentvalue{g})
    [~,tmp{end+1},~]=fileparts(handles.experimentlist{g}{handles.experimentvalue{g}(e)});
  end
  tmp{end+1}='';
end

h=msgbox(tmp,'Parameters','replace');
