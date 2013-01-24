function varargout = JAABAPlot(varargin)
%JAABAPLOT: Plot the output of JAABA
%
% This program is part of JAABA.
%
% JAABA: The Janelia Automatic Animal Behavior Annotator
% Copyright 2012, Kristin Branson, HHMI Janelia Farm Resarch Campus
% http://jaaba.sourceforge.net/
% bransonk@janelia.hhmi.org
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License (version 3 pasted in LICENSE.txt) for 
% more details.

% Last Modified by GUIDE v2.5 23-Jan-2013 17:53:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JAABAPlot_OpeningFcn, ...
                   'gui_OutputFcn',  @JAABAPlot_OutputFcn, ...
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
handles.classifierlist={};
handles.classifiervalue=1;
handles.configurations={};
handles.scorefiles={};
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
%handles.individuals=1;
handles.sexdata={};
handles.fps=nan;
handles.classify_forcecompute=false;
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
handles.featuretimeseries_windowradius=10;
handles.interestingfeaturehistograms_omitnan=1;
handles.interestingfeaturehistograms_omitinf=1;
handles.interestingfeaturehistograms_absdprime=1;
%handles.timeseries_tight=0;
handles.prefs_centraltendency=1;
handles.prefs_dispersion=1;
handles.prefs_timeseriesxoffset=1;
handles.prefs_convolutionwidth=1000;
%handles.logy=0;
%handles.stats=0;
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
%handles.colors={'r' 'g' 'b' 'c' 'm' 'y' 'k';
%                '#FF0000' '#00FF00' '#0000FF' '#00FFFF' '#FF00FF' '#FFFF00' '#000000'};
handles.colors=[];


% ---
function handles=load_configuration_file(filename,hObject,eventdata,handles)

handles_saved=load(filename);
handles_saved=handles_saved.handles;
handles.grouplist=handles_saved.grouplist;
handles.groupvalue=handles_saved.groupvalue;
handles.experimentlist=handles_saved.experimentlist;
handles.experimentvalue=handles_saved.experimentvalue;
handles.classifierlist=handles_saved.classifierlist;
handles.classifiervalue=handles_saved.classifiervalue;
handles.configurations=handles_saved.configurations;
handles.scorefiles=handles_saved.scorefiles;
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
handles.fps=handles_saved.fps;
handles.classify_forcecompute=handles_saved.classify_forcecompute;
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
handles.featuretimeseries_windowradius=handles_saved.featuretimeseries_windowradius;
handles.interestingfeaturehistograms_omitnan=handles_saved.interestingfeaturehistograms_omitnan;
handles.interestingfeaturehistograms_omitinf=handles_saved.interestingfeaturehistograms_omitinf;
handles.interestingfeaturehistograms_absdprime=handles_saved.interestingfeaturehistograms_absdprime;
%handles.timeseries_tight=handles_saved.timeseries_tight;
handles.prefs_centraltendency=handles_saved.prefs_centraltendency;
handles.prefs_dispersion=handles_saved.prefs_dispersion;
handles.prefs_timeseriesxoffset=handles_saved.prefs_timeseriesxoffset;
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
function handles=update_experiment_data(handles,features,sexdata,individuals)

handlesexperimentlist=[handles.experimentlist{:}];

handlesfeatures=cell(1,length(handlesexperimentlist));
handlessexdata=cell(1,length(handlesexperimentlist));
handlesindividuals=zeros(length(handlesexperimentlist),length(handles.scorefiles));

parfor ge=1:length(handlesexperimentlist)
  if(features)
    tmp=dir(fullfile(handlesexperimentlist{ge},'perframe','*.mat'));
    [handlesfeatures{ge}{1:length(tmp)}]=deal(tmp.name);
    handlesfeatures{ge}=cellfun(@(x) x(1:(end-4)),handlesfeatures{ge},'uniformoutput',false);
  end

  if(sexdata)
    tmp=load(fullfile(handlesexperimentlist{ge},'perframe','sex.mat'));
    cellfun(@(x) strcmp(x,'M'),tmp.data,'uniformoutput',false);
    handlessexdata(ge)={ans};
  end

  if(individuals)
    behavior_data=[];
    parfor_tmp=zeros(1,length(handles.scorefiles));
    for s=1:length(handles.scorefiles)
      classifier=load(handles.classifierlist{s});
      try
        behavior_data=load(fullfile(handlesexperimentlist{ge},handles.scorefiles{s}));
      catch
        behavior_data.allScores.t0s=[];
        behavior_data.timestamp=nan;
      end
      if ((classifier.classifierTS ~= behavior_data.timestamp) || (length(behavior_data.allScores.t0s)==0))
        parfor_tmp(s)=-1;
      else
        parfor_tmp(s)=length(behavior_data.allScores.t0s);
      end
    end
    handlesindividuals(ge,:)=parfor_tmp;
  end
end

if(features)
  handles.features={handlesfeatures{:}};
  handles.featurelist=check_for_diff_and_return_intersection(handles.features);
end

if(sexdata)
  handles.sexdata={handlessexdata{:}};
end

if(individuals)
  handles.individuals=handlesindividuals;
  handles=fillin_individuallist(handles);
end

if((isnan(handles.fps))&&(length(handles.classifierlist)>0))
  classifier=load(handles.classifierlist{1});
  t=load(fullfile(handlesexperimentlist{1},classifier.trxfilename));
  handles.fps=t.trx(1).fps;
end

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];


% ---
function update_figure(handles)

if(isempty(handles.grouplist) || length(handles.experimentlist{handles.groupvalue})==0)
  set(handles.ExperimentList,'enable','off');
  set(handles.ExperimentMove,'enable','off');
else
  set(handles.ExperimentList,'enable','on');
  set(handles.ExperimentMove,'enable','on');
end
if(isempty(handles.grouplist))
  set(handles.GroupList,'enable','off');
  set(handles.ExperimentAdd,'enable','off');
  set(handles.ExperimentDelete,'enable','off');
else
  set(handles.GroupList,'enable','on');
  set(handles.ExperimentAdd,'enable','on');
  set(handles.ExperimentDelete,'enable','on');
end
if(sum(cellfun(@length,handles.experimentlist))==0)
  set(handles.FeatureList,'enable','off');
  set(handles.ClassifierCheck,'enable','off');
else
  set(handles.FeatureList,'enable','on');
  set(handles.ClassifierCheck,'enable','on');
end
if(isempty(handles.classifierlist))
  set(handles.ClassifierList,'enable','off');
  set(handles.ClassifierDelete,'enable','off');
  set(handles.BehaviorList,'enable','off');
  set(handles.BehaviorLogic,'enable','off');
else
  set(handles.ClassifierList,'enable','on');
  set(handles.ClassifierDelete,'enable','on');
  set(handles.BehaviorList,'enable','on');
  set(handles.BehaviorLogic,'enable','on');
end
if(handles.behaviorlogic==1)
  set(handles.BehaviorList2,'enable','off');
else
  set(handles.BehaviorList2,'enable','on');
end
if((sum(cellfun(@length,handles.experimentlist))==0) || (isempty(handles.classifierlist)))
  set(handles.ClassifierClassify,'enable','off');
else
  set(handles.ClassifierClassify,'enable','on');
end
if((sum(cellfun(@length,handles.experimentlist))==0) || (isempty(handles.classifierlist)) || ...
    (sum(sum(handles.individuals==-1))>0) || (sum(sum(diff(handles.individuals,[],2)~=0))>0))
  set(handles.InterestingFeatureHistograms,'enable','off');
  set(handles.BehaviorBarChart,'enable','off');
  set(handles.BehaviorTimeSeries,'enable','off');
  set(handles.FeatureHistogram,'enable','off');
  set(handles.FeatureTimeSeries,'enable','off');
  set(handles.IndividualList,'enable','off');
else
  set(handles.InterestingFeatureHistograms,'enable','on');
  set(handles.BehaviorBarChart,'enable','on');
  set(handles.BehaviorTimeSeries,'enable','on');
  set(handles.FeatureHistogram,'enable','on');
  set(handles.FeatureTimeSeries,'enable','on');
  set(handles.IndividualList,'enable','on');
end

if(isempty(handles.grouplist))
  set(handles.GroupList,'String',{''},'Value',1);
else
  tmp=length(handles.grouplist);
%      {handles.colors{2,1+mod(0:length(handles.grouplist)-1,length(handles.colors))}}',...
  cellstr(strcat(repmat('<html><font color=#"',tmp,1),...
      reshape(dec2hex(round(handles.colors'*255),2)',6,size(handles.colors,1))',...
      repmat('">',tmp,1),handles.grouplist',repmat('</font></html>',tmp,1)));
  set(handles.GroupList,'String',ans,'Value',handles.groupvalue);
end
if(isempty(handles.experimentlist))
  set(handles.ExperimentList,'String',{''},'Value',1);
else
  set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue});
  set(handles.ExperimentList,'Value',handles.experimentvalue{handles.groupvalue});
end
if(isempty(handles.classifierlist))
  set(handles.ClassifierList,'String',{''},'Value',1);
else
  set(handles.ClassifierList,'String',handles.classifierlist);
  set(handles.ClassifierList,'Value',handles.classifiervalue);
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
%set(handles.Table,'Data',[]);

menu_classify_forcecompute_set(handles.classify_forcecompute);
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
menu_interestingfeaturehistograms_absdprime_set(handles.interestingfeaturehistograms_absdprime);
menu_prefscentraltendency_set(handles.prefs_centraltendency);
menu_prefsdispersion_set(handles.prefs_dispersion);
menu_prefstimeseriesxoffset_set(handles.prefs_timeseriesxoffset);


% --- Executes just before JAABAPlot is made visible.
function JAABAPlot_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

if ~isdeployed,
  SetUpJAABAPath;
end

if(exist('matlabpool')==2 && matlabpool('size')==0)
  
  useparallel = true;
  if isdeployed,
    if ispc || (isunix && ~ismac),
      filename = deployedRelative2Global('ParallelComputingConfiguration_Local_Win4.settings');
      if ~exist(filename,'file'),
        useparallel = false;
      else
        setmcruserdata('ParallelProfile','ParallelComputingConfiguration_Local_Win4.settings');
      end
    end
  end
  if useparallel
    matlabpool('open');
  end

end

if isdeployed,
  handles.rcfilename = deployedRelative2Global('most_recent_config.mat');
else
  handles.rcfilename = 'most_recent_config.mat';
end

if(exist(handles.rcfilename)==2)
  try
    handles=load_configuration_file(handles.rcfilename,hObject,eventdata,handles);
  catch ME
    errordlg({'Error loading last used configuration. Using default values.',getReport(ME)},'Error loading last used configuration');
    handles=initialize(handles);
  end
else
  handles=initialize(handles);
end
update_figure(handles);


% Choose default command line output for JAABAPlot
handles.output = hObject;

set(hObject,'CloseRequestFcn',@figure_CloseRequestFcn);
set(handles.Table,'CellSelectionCallback',@CellSelectionCallback);
set(handles.ExperimentList,'Callback',@ExperimentListCallback);
set(handles.ClassifierList,'Callback',@ClassifierListCallback);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes JAABAPlot wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function ExperimentListCallback(hObject,eventdata)

handles=guidata(hObject);
handles.experimentvalue{handles.groupvalue}=get(handles.ExperimentList,'Value');
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


function ClassifierListCallback(hObject,eventdata)

handles=guidata(hObject);
handles.classifiervalue=get(handles.ClassifierList,'Value');
handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];
guidata(hObject,handles);


% ---
function figure_CloseRequestFcn(hObject, eventdata)

handles=guidata(hObject);
try
  save(handles.rcfilename,'handles');
catch ME,
  errordlg({'Error saving last configuration. State will not be saved.',getReport(ME)},'Error saving last configuration');
end
delete(hObject);


% --- Outputs from this function are returned to the command line.
function varargout = JAABAPlot_OutputFcn(hObject, eventdata, handles)
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
function feature_intersection=check_for_diff_and_return_intersection(arg)

if(length(arg)==0)  ret_val=[];  return;  end
if(length(arg)==1)  ret_val=arg{1};  return;  end

feature_union=unique([arg{:}]);
feature_intersection=arg{1};
for i=2:length(arg)
  feature_intersection=intersect(feature_intersection,arg{i});
end
if(numel(feature_intersection)<numel(feature_union))
  char(setdiff(feature_union,feature_intersection));
  [ans repmat(', ',size(ans,1),1)];
  reshape(ans',1,numel(ans));
  uiwait(errordlg(['feature(s) ' ans(1:end-2) ' are/is not in all experiments and so will be ignored.']));
end


% ---
function handles=fillin_individuallist(handles)

if((numel(handles.individuals)==0) || (sum(sum(diff(handles.individuals,[],2)~=0))>0))  return;  end

cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];

tmp=cell(1,3+sum(handles.individuals(:,1)));
tmp(1:3)={'All' 'Male' 'Female'};
k=4;
for e=1:size(handles.individuals,1)
  g=find(cumsum_num_exp_per_group<e,1,'last');
  for i=1:handles.individuals(e,1)
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

newexperiments=uipickfiles('prompt','Select experiment directory','filterspec',directory);
if(~iscell(newexperiments) || (length(newexperiments)==0))  return;  end
if((length(newexperiments)==1)&&(exist(newexperiments{1})==2))
  newexperiments=textread(newexperiments{1},'%s');
end
tmp=ismember(newexperiments,[handles.experimentlist{:}]);
if(sum(tmp)>0)
  msg{1}='The following experiments have already been added:';
  msg{2}='';
  msg(3:(2+sum(tmp)))=newexperiments(tmp);
  uiwait(errordlg(msg));
  newexperiments(tmp)=[];
end

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

[directory,~,~]=fileparts(newexperiments{1});
handles.experimentlist{handles.groupvalue}={handles.experimentlist{handles.groupvalue}{:} newexperiments{:}};
handles.experimentvalue{handles.groupvalue}=1:length(handles.experimentlist{handles.groupvalue});

handlesfeatures=cell(1,length(newexperiments));
handlessexdata=cell(1,length(newexperiments));
handlesindividuals=zeros(length(newexperiments),length(handles.scorefiles));
parfor n=1:length(newexperiments)
%for n=1:length(newexperiments)
  tmp=dir(fullfile(newexperiments{n},'perframe','*.mat'));
  [handlesfeatures{n}{1:length(tmp)}]=deal(tmp.name);
  handlesfeatures{n}=cellfun(@(x) x(1:(end-4)),handlesfeatures{n},'uniformoutput',false);

  sexdata=load(fullfile(newexperiments{n},'perframe','sex.mat'));
  sexdata.data=cellfun(@(x) strcmp(x,'M'),sexdata.data,'uniformoutput',false);
  handlessexdata(n)={sexdata.data};

  behavior_data=[];
  parfor_tmp=zeros(1,length(handles.scorefiles));
  for s=1:length(handles.scorefiles)
    classifier=load(handles.classifierlist{s});
    try
      behavior_data=load(fullfile(newexperiments{n},handles.scorefiles{s}));
    catch
      behavior_data.allScores.t0s=[];
      behavior_data.timestamp=nan;
    end
    if ((classifier.classifierTS ~= behavior_data.timestamp) || (length(behavior_data.allScores.t0s)==0))
      parfor_tmp(s)=-1;
    else
      parfor_tmp(s)=length(behavior_data.allScores.t0s);
    end
  end
  handlesindividuals(n,:)=parfor_tmp;
end

if((isnan(handles.fps))&&(length(handles.classifierlist)>0))
  classifier=load(handles.classifierlist{1});
  t=load(fullfile(newexperiments{1},classifier.trxfilename));
  handles.fps=t.trx(1).fps;
end

%msg{1}='can''t find the configuration file for the following experiments and classifiers:';
%for n=1:length(newexperiments)
%  tmp=find(cellfun(@isempty,handlesbehaviors{n}));
%  if(length(tmp)>0)
%    msg{end+1}='';
%    msg{end+1}=newexperiments{n};
%    for i=1:length(tmp)
%      msg{end+1}=handlesclassifiers{n}{tmp(i)};
%    end
%  end
%end
%if(length(msg)>1)  uiwait(errordlg(msg));  end

handles.features={handles.features{:} handlesfeatures{:}};
handles.sexdata={handles.sexdata{:} handlessexdata{:}};
handles.individuals=[handles.individuals; handlesindividuals];

tmp=length(handles.features);
idx=[1 : (sum(cellfun(@length,handles.experimentlist(1:handles.groupvalue)))-length(newexperiments)) ...
    ((tmp-length(newexperiments)+1) : tmp) ...
    ((tmp-length(newexperiments)-sum(cellfun(@length,handles.experimentlist((handles.groupvalue+1):end)))+1) : ...
        (tmp-length(newexperiments)))];
handles.features=handles.features(idx);
handles.sexdata=handles.sexdata(idx);
handles.individuals=handles.individuals(idx,:);

handles.featurelist=check_for_diff_and_return_intersection(handles.features);
%set(handles.FeatureList,'String',handles.featurelist);
handles=fillin_individuallist(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

update_figure(handles);

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ExperimentDelete.
function ExperimentDelete_Callback(hObject, eventdata, handles)
% hObject    handle to ExperimentDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

if(length(handles.experimentlist)==0)  return;  end

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

idx=handles.experimentvalue{handles.groupvalue};
handles.experimentlist{handles.groupvalue}(idx)=[];
handles.experimentvalue{handles.groupvalue}=1:length(handles.experimentlist{handles.groupvalue});

idx=idx+sum(cellfun(@length,handles.experimentlist(1:(handles.groupvalue-1))));
handles.features(idx)=[];
handles.sexdata(idx)=[];
handles.individuals(idx,:)=[];

if(isempty(handles.experimentlist{handles.groupvalue}))
  handles.colors=handles.colors(setdiff(1:size(handles.colors,1),handles.groupvalue),:);
  handles.experimentlist(handles.groupvalue)=[];
  handles.experimentvalue(handles.groupvalue)=[];
  handles.grouplist(handles.groupvalue)=[];
  handles.groupvalue=max(1,min([handles.groupvalue length(handles.grouplist)]));
end

handles.featurelist=check_for_diff_and_return_intersection(handles.features);
handles=fillin_individuallist(handles);

update_figure(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


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

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

from_group=handles.groupvalue;
idx=handles.experimentvalue{from_group};
idxF=idx+sum(cellfun(@length,handles.experimentlist(1:(from_group-1))));
idxT=sum(cellfun(@length,handles.experimentlist(1:to_group)));

handles.experimentlist{to_group}=...
    {handles.experimentlist{to_group}{:} handles.experimentlist{from_group}{idx}};
handles.experimentvalue{to_group}=1:length(handles.experimentlist{to_group});
handles.experimentlist{from_group}(idx)=[];
handles.experimentvalue{from_group}=1:length(handles.experimentlist{from_group});
%if(isempty(handles.experimentlist{from_group}))
%  set(handles.ExperimentList,'String',{''},'Value',1);
%else
%  set(handles.ExperimentList,'String',handles.experimentlist{from_group});
%  set(handles.ExperimentList,'Value',handles.experimentvalue{from_group});
%end
if(idxF(1)<=idxT)  idxT=idxT-length(idxF);  end
tmp=setdiff(1:length(handles.features),idxF);
tmp=[tmp(1:idxT) idxF tmp((idxT+1):end)];
%handles.behaviors=handles.behaviors(tmp);
handles.features=handles.features(tmp);
handles.sexdata=handles.sexdata(tmp);
handles.individuals=handles.individuals(tmp,:);

handles=fillin_individuallist(handles);

update_figure(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


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
update_figure(handles);
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

inputdlg({'Name:'},'Create new experiment group');
if(isempty(ans))  return;  end
defaultcolors=[1 0 0;  0 0.5 0;  0 0 1;  0 1 1;  1 0 1;  0.749 0.749 0;  0 0 0];
handles.colors(end+1,:)=uisetcolor(defaultcolors(1+mod(length(handles.grouplist),size(defaultcolors,1)),:));
handles.grouplist{end+1}=char(ans);
handles.groupvalue=length(handles.grouplist);
handles.experimentlist{handles.groupvalue}={};
handles.experimentvalue{handles.groupvalue}=[];
%tmp=length(handles.grouplist);
%    {handles.colors{2,1+mod(0:length(handles.grouplist)-1,length(handles.colors))}}',...
%cellstr(strcat(repmat('<html><font color="',tmp,1),...
%    reshape(dec2hex(round(handles.colors'*255),2)',6,size(handles.colors,1))',...
%    repmat('">',tmp,1),handles.grouplist',repmat('</font></html>',tmp,1)));
%set(handles.GroupList,'String',ans,'Value',handles.groupvalue);
%set(handles.ExperimentList,'String',handles.experimentlist{handles.groupvalue},...
%    'Value',handles.experimentvalue{handles.groupvalue});
%set(handles.GroupList,'enable','on');
%set(handles.ExperimentAdd,'enable','on');
%set(handles.ExperimentDelete,'enable','on');
%set(handles.ExperimentMove,'enable','on');
update_figure(handles);
guidata(hObject,handles);


% ---
function handles=classifier_add(handles,newclassifiers)

handlesconfigurations=cell(1,length(newclassifiers));
handlesbehaviorlist=cell(1,length(newclassifiers));
handlesscorefiles=cell(1,length(newclassifiers));
handlesindividuals=zeros(sum(cellfun(@length,handles.experimentlist)),length(newclassifiers));
for c=1:length(newclassifiers)
  classifier=load(newclassifiers{c});
  if(~isfield(classifier,'postprocessparams'))
    uiwait(errordlg(['not a valid classifier file.  skipping ' newclassifiers{c}],''));
    newclassifiers{c}='';
    continue;
  end
  params=[];  params.behaviors.names='';
  try
    handlesconfigurations{c}=classifier.configfilename;
    [~,~,ext] = fileparts(classifier.configfilename);
    if strcmpi(ext,'.xml');
      params=ReadXMLConfigParams(handlesconfigurations{c});
    else
      params = load(handlesconfigurations{c});
    end
  catch
    directory = fileparts(newclassifiers{c});
    fullfile(directory,classifier.configfilename);
    if(exist(ans,'file'))
      [pa,na,ex]=fileparts(ans);
    else
      pa=directory;  na='';  ex='';
    end
    [~,tmp,~]=fileparts(newclassifiers{c});
    [configfile tmp]=uigetfile(pa,['Select configuration file for ' tmp],[na ex]);
    if(isnumeric(configfile) && isnumeric(tmp) && (configfile==0) && (tmp==0))
      uiwait(errordlg(['skipping ' newclassifiers{c}],''));
      newclassifiers{c}='';
      continue;
    else
      try
        handlesconfigurations{c}=fullfile(tmp,configfile);
        if strcmpi(ext,'.xml');
          params=ReadXMLConfigParams(handlesconfigurations{c});
        else
          params = load(handlesconfigurations{c});
        end
      catch
        uiwait(errordlg(['problem loading config file.  skipping ' newclassifiers{c}],''));
        newclassifiers{c}='';
        continue;
      end
    end
  end
  if(~isfield(params.behaviors,'names') || ~isfield(params.file,'scorefilename'))
    uiwait(errordlg(['not a valid config file.  skipping ' newclassifiers{c}],''));
    newclassifiers{c}='';
    continue;
  end

  if iscell(params.behaviors.names),
    handlesbehaviorlist{c} = params.behaviors.names{1};
  else
    handlesbehaviorlist{c}=params.behaviors.names;
  end
  handlesscorefiles{c}=params.file.scorefilename;

  ee=0;  behavior_data=[];  parfor_tmp=zeros(sum(cellfun(@length,handles.experimentlist)),1);
  for g=1:length(handles.grouplist)
    for e=1:length(handles.experimentlist{g})
      try
        behavior_data=load(fullfile(handles.experimentlist{g}{e},handlesscorefiles{c}));
      catch
        behavior_data.allScores.t0s=[];
        behavior_data.timestamp=nan;
      end
      if((classifier.classifierTS ~= behavior_data.timestamp) || (length(behavior_data.allScores.t0s)==0))
        parfor_tmp(ee+e)=-1;
      else
        parfor_tmp(ee+e)=length(behavior_data.allScores.t0s);
      end
    end
    ee=ee+length(handles.experimentlist{g});
  end
  handlesindividuals(:,c)=parfor_tmp;

  if((c==1)&&(isnan(handles.fps)))
    handlesexperimentlist=[handles.experimentlist{:}];
    if(length(handlesexperimentlist)>0)
      t=load(fullfile(handlesexperimentlist{1},classifier.trxfilename));
      handles.fps=t.trx(1).fps;
    end
  end

end

idx=find(~cellfun(@isempty,newclassifiers));
handles.classifierlist={handles.classifierlist{:} newclassifiers{idx}};
handles.configurations={handles.configurations{:} handlesconfigurations{idx}};
handles.behaviorlist={handles.behaviorlist{1:(end-1)} handlesbehaviorlist{idx} 'All'};
handles.scorefiles={handles.scorefiles{:} handlesscorefiles{idx}};
handles.individuals=[handles.individuals handlesindividuals(:,idx)];

handles.classifiervalue=1:length(handles.classifierlist);

handles=fillin_individuallist(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];


% --- Executes on button press in ClassifierAdd.
function ClassifierAdd_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierAdd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

persistent directory
if(isempty(directory))  directory=pwd;  end

tmp=directory;
[newclassifiers directory]=uigetfile(directory,'Select classifier files','multiselect','on');
if(isnumeric(newclassifiers)&&(newclassifiers==0))  directory=tmp; return;  end
if(~iscell(newclassifiers))  newclassifiers={newclassifiers};  end
newclassifiers=cellfun(@(x) fullfile(directory,x),newclassifiers,'uniformoutput',false);

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handles=classifier_add(handles,newclassifiers);
update_figure(handles);

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierDelete.
function ClassifierDelete_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierDelete (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

idx=handles.classifiervalue;
handles.classifierlist(idx)=[];
handles.classifiervalue=1:max(1,length(handles.classifierlist));
handles.configurations(idx)=[];
handles.behaviorlist(idx)=[];
handles.behaviorvalue=max(1,length(handles.behaviorlist));
handles.scorefiles(idx)=[];
handles.individuals(:,idx)=[];

if(length(handles.classifierlist)==0)
  handles.fps=nan;
end

handles=fillin_individuallist(handles);
update_figure(handles);

handles.interestingfeaturehistograms_cache=[];
handles.interestingfeaturetimeseries_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierAuto.
function ClassifierAuto_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierAuto (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

classifiers_found=cell(1,length(handlesexperimentlist));
classifiers_notfound=cell(1,length(handlesexperimentlist));
parfor ge=1:length(handlesexperimentlist)
  tmp=dir(fullfile(handlesexperimentlist{ge},'*.mat'));
  possiblescorefiles=setdiff({tmp.name},handles.scorefiles);
  classifiers_found{ge}={};
  classifiers_notfound{ge}={};
  for p=1:length(possiblescorefiles)
    tmp=load(fullfile(handlesexperimentlist{ge},possiblescorefiles{p}));
    if(~isfield(tmp,'classifierfilename'))
      continue;
    end
    if(exist(tmp.classifierfilename)==2)
      classifiers_found{ge}{end+1}=tmp.classifierfilename;
    else
      classifiers_notfound{ge}{end+1}=tmp.classifierfilename;
    end
  end
end
classifiers_found=unique([classifiers_found{:}]);
classifiers_notfound=unique([classifiers_notfound{:}]);

classifiers_found=setdiff(classifiers_found,handles.classifierlist);
handles=classifier_add(handles,classifiers_found);

if(length(classifiers_notfound)>2)
  uiwait(errordlg({'Could not find these classifiers:' '' [classifiers_notfound{:}]}));
end

update_figure(handles);

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierCheck.
function ClassifierCheck_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierCheck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

table=cell(length(handlesexperimentlist),length(handles.scorefiles));

parfor ge=1:length(handlesexperimentlist)
%for ge=1:length(handlesexperimentlist)
  behavior_data=[];
  parfor_tmp=cell(1,length(handles.scorefiles));
  for c=1:length(handles.classifierlist)
    classifier=load(handles.classifierlist{c});
    try
      behavior_data=load(fullfile(handlesexperimentlist{ge},handles.scorefiles{c}));
    catch ME,
      warning(getReport(ME));
      parfor_tmp{c}='missing';
      continue;
    end
    if (classifier.classifierTS ~= behavior_data.timestamp)
      parfor_tmp{c}=[datestr(classifier.classifierTS) ' ~= ' datestr(behavior_data.timestamp)];
    else
      parfor_tmp{c}='up-to-date';
    end
  end
  table(ge,:)=parfor_tmp;
  tmp=dir(fullfile(handlesexperimentlist{ge},'*.mat'));
  extrascorefiles=setdiff({tmp.name},handles.scorefiles);
  table2{ge}={};
  for s=1:length(extrascorefiles)
    tmp=load(fullfile(handlesexperimentlist{ge},extrascorefiles{s}));
    if(~isfield(tmp,'classifierfilename'))
      continue;
    end
    table2{ge}{end+1}=extrascorefiles{s};
  end
end

[~,tmp,~]=cellfun(@(x) fileparts(x),handlesexperimentlist,'uniformoutput',false);
table=[tmp' table];
[~,tmp,~]=cellfun(@(x) fileparts(x),handles.classifierlist,'uniformoutput',false);
table=[{'' tmp{:}}; table];

tmp=unique([table2{:}]);
if(length(tmp)>0)
  table{1,end+1}='   ';
  table(1,(end+1):(end+length(tmp)))=tmp;
end
for i=1:length(table2)
  for j=1:length(table2{i})
    idx=find(cellfun(@(x) strcmp(x,table2{i}{j}),table(1,:)));
    table{i+1,idx}='extra';
  end
end

if(length(handles.experimentlist)>1)
  %tmp=zeros(1,length(handles.experimentlist)+length([handles.experimentlist{:}]));
  tmp=zeros(1,1+length([handles.experimentlist{:}]));
  cumsum(cellfun(@length,handles.experimentlist));
  tmp(2+ans(1:end-1))=1;
  cumsum(tmp);
  %tmp=ans+(1:(-1+length(handles.experimentlist)+length([handles.experimentlist{:}])));
  tmp=ans+(1:(1+length([handles.experimentlist{:}])));
  tmp2(tmp,:)=[table(:,:)];
  table=tmp2;
end

set(handles.Table,'Data',table);
6*max(cellfun(@length,table),[],1);
set(handles.Table,'ColumnWidth',mat2cell(ans,1,ones(1,length(ans))));
%set(handles.Table,'ColumnWidth','auto');
set(handles.Table,'ColumnName',{''});
set(handles.Table,'RowName',{});

%handles.table_data=table;
%handles.table='classifier';

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in ClassifierClassify.
function ClassifierClassify_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierClassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

%h=waitbar(0,'This will likely take awhile...',...
%    'CreateCancelBtn','fid=fopen(fullfile(tempdir,''cancel.txt''),''w''); fclose(fid);');
h=waitbar(0,'This will likely take awhile...');
fid=fopen(fullfile(tempdir,'progressbar.txt'),'w');
fwrite(fid,length(handles.behaviorlist)*sum(cellfun(@length,handles.experimentvalue)),'uint32');
fclose(fid);
t = timer('TimerFcn',{@progress_bar,h}, 'Period', 3, 'ExecutionMode', 'fixedRate');
start(t);

handlesexperimentlist=[handles.experimentlist{:}];

parfor ge=1:length(handlesexperimentlist)
%for ge=1:length(handlesexperimentlist)
  JAABADetect(handlesexperimentlist{ge},...
      'classifierfiles',handles.classifierlist(handles.classifiervalue),...
      'configfiles',handles.configurations(handles.classifiervalue),...
      'forcecompute',handles.classify_forcecompute);
  fid=fopen(fullfile(tempdir,'progressbar.txt'),'a');  fwrite(fid,1,'uint32');  fclose(fid);
end

stop(t);
delete(t);
delete(h);
delete(fullfile(tempdir,'progressbar.txt'));

%if(exist(fullfile(tempdir,'cancel.txt')))
%  delete(fullfile(tempdir,'cancel.txt'));
%  set(handles.Status,'string','Ready.','foregroundcolor','g');
%  drawnow;
%  set(handles.figure1,'pointer','arrow');
%  return;
%end

handles=update_experiment_data(handles,false,false,true);
update_figure(handles);

handles.interestingfeaturehistograms_cache=[];

guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on selection change in ClassifierList.
function ClassifierList_Callback(hObject, eventdata, handles)
% hObject    handle to ClassifierList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ClassifierList contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ClassifierList


% --- Executes during object creation, after setting all properties.
function ClassifierList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ClassifierList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
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
function h=plot_it(xdata,ydata,style,centraltendency,dispersion,color,linewidth,fid,experimentlist)

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
    h=plot(xdata,data_ct,'color',color,'linewidth',linewidth);
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
    h=plot(xdata,data_ct,'color',color,'linewidth',3*linewidth);
    idx=isnan(data_dp) | isnan(data_dn);
    xdata=xdata(~idx);  data_dp=data_dp(~idx);  data_dn=data_dn(~idx);
    color2=(get(h,'color')+[4 4 4])/5;
    k=1;  m=0;  step=10000;
    while(k<=length(xdata))
      idx=k:min(k+step,length(xdata));
      patch([xdata(idx) fliplr(xdata(idx))],[data_dp(idx) fliplr(data_dn(idx))],color2,'edgecolor','none');
      k=k+step+1;  m=m+1;
    end
    get(gca,'children');  set(gca,'children',circshift(ans,-m));  % send to back
    %get(gca,'children');  set(gca,'children',ans([m+1 1:m (m+2):end]));  % send to back
    fprintf(fid,['%% ydata, ' str_dp '\n']);  fprintf(fid,'%g, ',data_dp);  fprintf(fid,'\n');
    fprintf(fid,['%% ydata, ' str_dn '\n']);  fprintf(fid,'%g, ',data_dn);  fprintf(fid,'\n');
    fprintf(fid,['%% ydata, ' str_ct '\n']);  fprintf(fid,'%g, ',data_ct);  fprintf(fid,'\n');
  case 3
    h=plot(xdata,ydata','color',color,'linewidth',linewidth);
    h=h(1);
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
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

during={};  not_during={};
for i=1:length(behavior_data.allScores.t0s)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end
  tmp1=zeros(1,length(feature_data.data{i}));
  tmp1(behavior_data.allScores.t0s{i}-behavior_data.allScores.tStart(i)+1)=1;
  tmp1(behavior_data.allScores.t1s{i}-behavior_data.allScores.tStart(i)+1)=-1;
  tmp1=logical(cumsum(tmp1(1:length(feature_data.data{i}))));
  
  if(behavior_logic>1)
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
    start=find(~partition_idx(1:(end-1)) &  partition_idx(2:end));
    stop =find( partition_idx(1:(end-1)) & ~partition_idx(2:end));
    during{i}=[];  not_during{i}=[];
    if(length(start)>0)
      for j=1:length(start)
        if(sum(sexdata{i}(start(j):(stop(j)-1))) < ((stop(j)-start(j))/2))  continue;  end
        switch(perwhat)
          case 2
            during{i}(j)=mean(feature_data.data{i}(start(j):(stop(j)-1)));
          case 3
            during{i}(j)=median(feature_data.data{i}(start(j):(stop(j)-1)));
          case 4
            during{i}(j)=max(feature_data.data{i}(start(j):(stop(j)-1)));
          case 5
            during{i}(j)=min(feature_data.data{i}(start(j):(stop(j)-1)));
          case 6
            during{i}(j)=std(feature_data.data{i}(start(j):(stop(j)-1)));
        end
      end
      if(length(start)>1)
        for j=1:(length(start)-1)
          if(sum(sexdata{i}(stop(j):(start(j+1)-1))) < ((start(j+1)-stop(j))/2))  continue;  end
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


% --- Executes on button press in FeatureHistogram.
function FeatureHistogram_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureHistogram (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals(:,1))'];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

ggee=1:length(handlesexperimentlist);
individual=handles.individualvalue;
if(individual>3)
  ggee=find(cumsum_num_indi_per_exp<(individual-3),1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='feature histogram';

h=figure('toolbar','figure');  hold on;
guidata(h,handles);
uicontrol(h,'style','pushbutton','string','?','position',[10 10 20 20],'callback',@figure_info_callback);

bb=handles.behaviorvalue;
if(bb==length(handles.behaviorlist))  bb=1:(bb-1);  end

behavior_logic=handles.behaviorlogic;
score_file2=handles.scorefiles{handles.behaviorvalue2};
feature_value=handles.featurevalue;
feature_list=handles.featurelist;
notduring=handles.featurehistogram_notduring;
nbins=handles.featurehistogram_nbins;
style=handles.featurehistogram_style;
centraltendency=handles.prefs_centraltendency;
dispersion=handles.prefs_dispersion;

h=[];
for b=bb

  tstr=char(strrep(handles.behaviorlist(b),'_','-'));
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
  units=load(fullfile(handlesexperimentlist{ggee(1)},'perframe',...
      [feature_list{feature_value} '.mat']),'units');
  xstr=get_label(feature_list(feature_value),units.units);
  ystr='normalized';

  print_csv_help(fid,handles.type,tstr,xstr,ystr);

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    subplot(ceil(length(bb)/ans),ans,b);
  end
  hold on;

  score_file=handles.scorefiles{b};

  during=cell(1,length(ggee));
  not_during=cell(1,length(ggee));
  %for ge=ggee
  parfor ge=ggee
    if((individual<4)&&(~ismember(ge,selected_exp)))
      during{ge}=nan;
      not_during{ge}=nan;
      continue;
    end

    behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
    if(behavior_logic>1)
      behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
    else
      behavior_data2=[];
    end
    feature_data=load(fullfile(handlesexperimentlist{ge},'perframe',...
          [feature_list{feature_value} '.mat']));

    tmp2=handles.sexdata{ge};
    for i=1:length(tmp2)
      switch(individual)
        case(2)
        case(3)
          tmp2{i}=~tmp2{i};
        otherwise
          tmp2{i}=ones(1,length(tmp2{i}));
      end
    end
    tmploop=nan;  if(individual>3)  tmploop=individual-3;  end

    [during{ge} not_during{ge}]=calculate_feature_histogram(...
        behavior_data,behavior_logic,behavior_data2,feature_data,tmp2,tmploop,...
        handles.featurehistogram_perwhat);
  end

  i=1;
  while i<=length(during)
    if((length(during{i})==1) && isnan(during{i}))
      during(i)=[];
      not_during(i)=[];
    else
      i=i+1;
    end
  end

  max(cellfun(@(x) size(x,2),during));
  cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],during,'uniformoutput',false);
  during=cat(1,ans{:});
  max(cellfun(@(x) size(x,2),not_during));
  cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],not_during,'uniformoutput',false);
  not_during=cat(1,ans{:});

  low=[];  high=[];  nearzero=[];
  if(~isempty(during))
    low=min(min(during));
    high=max(max(during));
    unique(reshape(abs(during),1,prod(size(during))));
    nearzero=ans(1);  if(ans(1)==0)  nearzero=ans(2);  end
  end
  if(notduring && ~isempty(not_during))
    low=min([low min(not_during)]);
    high=max([high max(not_during)]);
    unique(reshape(abs(not_during),1,prod(size(not_during))));
    tmp=ans(1);  if(ans(1)==0)  tmp=ans(2);  end
    nearzero=min(tmp,nearzero);
  end

  if(~isempty(low) && ~isempty(high) && ~isempty(nearzero))
    if(handles.featurehistogram_logbinsize)
      if((low>=0) && (high>0))
        bins=logspace(log10(max(low,nearzero)),log10(high),nbins);
      elseif((low<0) && (high<=0))
        bins=fliplr(-logspace(log10(max(abs(high),nearzero)),log10(abs(low)),nbins));
      elseif((low<0) && (high>0))
        bins=[fliplr(-logspace(log10(nearzero),log10(abs(low)),nbins)) ...
            logspace(log10(nearzero),log10(abs(high)),nbins)];
      end
      plot(bins,zeros(1,length(bins)),'k.');
    else
      bins=linspace(low,high,nbins);
    end
  end

  for g=1:length(handles.grouplist)
    color=handles.colors(g,:);

    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);

    if(individual<4)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end

    if(notduring)
      if(~isempty(not_during(idx,:)))
        hist_not_during=hist(not_during(idx,:)',bins);
        if(size(not_during,1)==1)  hist_not_during=hist_not_during';  end
        hist_not_during.*repmat(([0 diff(bins)]+[diff(bins) 0])'/2,1,size(hist_not_during,2));
        hist_not_during=hist_not_during./repmat(sum(ans,1),size(hist_not_during,1),1);
        plot_it(bins,hist_not_during',style,centraltendency,dispersion,color,1,...
          fid,handlesexperimentlist(idx));
      else
        plot_it(nan,nan,style,centraltendency,dispersion,color,1,...
          fid,handlesexperimentlist(idx));
      end
    end
    if(~isempty(during(idx,:)))
      hist_during=hist(during(idx,:)',bins);
      if(size(hist_during,1)==1)  hist_during=hist_during';  end
      hist_during.*repmat(([0 diff(bins)]+[diff(bins) 0])'/2,1,size(hist_during,2));
      hist_during=hist_during./repmat(sum(ans,1),size(hist_during,1),1);
      linewidth=1;  if(notduring)  linewidth=2;  end
      h(g)=plot_it(bins,hist_during',style,centraltendency,dispersion,color,linewidth,...
          fid,handlesexperimentlist(idx));
    else
      h(g)=plot_it(nan,nan,style,centraltendency,dispersion,color,1,...
          fid,handlesexperimentlist(idx));
    end
  end

  title(tstr);
  xlabel(xstr);
  ylabel(ystr);
  axis tight;  zoom reset;
end

idx=find(h>0);
if(individual<4)
  legend(h(idx),handles.grouplist);
else
  legend(h(idx),handles.individuallist(individual));
end

fclose(fid);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;

%guidata(hObject,handles);


function progress_bar(~,~,h)

fid=fopen(fullfile(tempdir,'progressbar.txt'),'r');
den=fread(fid,1,'uint32');
fseek(fid,0,1);
num=ftell(fid)/4;
waitbar(num/den,h);
fclose(fid);
drawnow;


% --- Executes on button press in InterestingFeatureHistograms.
function InterestingFeatureHistograms_Callback(hObject, eventdata, handles)
% hObject    handle to InterestingFeatureHistograms (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

if(isempty(handles.interestingfeaturehistograms_cache))
  table_data={};

  cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
  mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
  cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
  selected_exp=[ans{:}];
  cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

  handlesexperimentlist=[handles.experimentlist{:}];
  handlesexperimentlist=handlesexperimentlist(selected_exp);

  nexperiments=length(handlesexperimentlist);
  nbehaviors=length(handles.behaviorlist)-1;
  nfeatures=length(handles.featurelist);

  h=waitbar(0,'This will likely take awhile...',...
      'CreateCancelBtn','fid=fopen(fullfile(tempdir,''cancel.txt''),''w''); fclose(fid);');
  fid=fopen(fullfile(tempdir,'progressbar.txt'),'w');
  fwrite(fid,nfeatures*sum(cellfun(@length,handles.experimentvalue)),'uint32');
  fclose(fid);
  t = timer('TimerFcn',{@progress_bar,h}, 'Period', 3, 'ExecutionMode', 'fixedRate');
  start(t);

  table_data=zeros(nexperiments,nbehaviors,nfeatures,6);
  parfor ge=1:nexperiments
  %for ge=1:nexperiments
    behavior_data={};
    for b=1:nbehaviors
      behavior_data{b}=load(fullfile(handlesexperimentlist{ge},handles.scorefiles{b}));
    end
    bad{ge}={};
    parfor_tmp=zeros(nbehaviors,nfeatures,6);
    for f=1:nfeatures
      if(exist(fullfile(tempdir,'cancel.txt')))  break;  end
      feature_data=load(fullfile(handlesexperimentlist{ge},'perframe',...
          [handles.featurelist{f} '.mat']));
      sexdata={};
      for s=1:length(feature_data.data)
        sexdata{s}=ones(1,length(feature_data.data{s}));
      end
      for b=1:nbehaviors
        if(exist(fullfile(tempdir,'cancel.txt')))  break;  end

        [during not_during]=calculate_feature_histogram(behavior_data{b},1,[],...
            feature_data,sexdata,nan,handles.featurehistogram_perwhat);
        parfor_tmp(b,f,:)=[mean(during) mean(not_during) ...
            std(during) std(not_during) length(during) length(not_during)];
      end
      fid=fopen(fullfile(tempdir,'progressbar.txt'),'a');  fwrite(fid,1,'uint32');  fclose(fid);
    end
    table_data(ge,:,:,:)=parfor_tmp;
  end

  stop(t);
  delete(t);
  delete(h);
  delete(fullfile(tempdir,'progressbar.txt'));

  if(exist(fullfile(tempdir,'cancel.txt')))
    delete(fullfile(tempdir,'cancel.txt'));
    set(handles.Status,'string','Ready.','foregroundcolor','g');
    set(handles.figure1,'pointer','arrow');
    drawnow;
    return;
  end

  tmp2=[];
  for g=1:length(handles.grouplist)
    gg=cumsum_num_selexp_per_group(g)+(1:length(handles.experimentlist{g}));
    tmp2=[tmp2; ...
        repmat(g,nbehaviors*nfeatures,1) ...
        reshape(squeeze(sum(table_data(gg,:,:,5))),nbehaviors*nfeatures,1) ...
        nan(nbehaviors*nfeatures,1) ...
        reshape(squeeze(sum(table_data(gg,:,:,6))),nbehaviors*nfeatures,1) ...
        reshape(repmat(1:nbehaviors,nfeatures,1),nbehaviors*nfeatures,1) ...
        reshape(repmat(1:nfeatures,nbehaviors,1)',nbehaviors*nfeatures,1) ...
        reshape(squeeze((mean(table_data(gg,:,:,1))-mean(table_data(gg,:,:,2)))./ ...
          sqrt((mean(table_data(gg,:,:,3).^2)+mean(table_data(gg,:,:,4).^2))/2))', ...
          nbehaviors*nfeatures,1)];

    if(g==length(handles.grouplist))  break;  end
    for g2=(g+1):length(handles.grouplist)
      gg2=cumsum_num_selexp_per_group(g2)+(1:length(handles.experimentlist{g2}));
      tmp2=[tmp2; ...
          repmat(g,nbehaviors*nfeatures,1) ...
          reshape(squeeze(sum(table_data(gg,:,:,5))),nbehaviors*nfeatures,1) ...
          repmat(g2,nbehaviors*nfeatures,1) ...
          reshape(squeeze(sum(table_data(gg2,:,:,5))),nbehaviors*nfeatures,1) ...
          reshape(repmat(1:nbehaviors,nfeatures,1),nbehaviors*nfeatures,1) ...
          reshape(repmat(1:nfeatures,nbehaviors,1)',nbehaviors*nfeatures,1) ...
          reshape(squeeze((mean(table_data(gg,:,:,1))-mean(table_data(gg2,:,:,1)))./ ...
            sqrt((mean(table_data(gg,:,:,3).^2)+mean(table_data(gg2,:,:,3).^2))/2))', ...
            nbehaviors*nfeatures,1)];
    end
  end
  handles.interestingfeaturehistograms_cache=tmp2;
else
  tmp2=handles.interestingfeaturehistograms_cache;
end

if(handles.interestingfeaturehistograms_omitnan)
  idx=find(~isnan(tmp2(:,7)));
  tmp2=tmp2(idx,:);
end
if(handles.interestingfeaturehistograms_omitinf)
  idx=find(~isinf(tmp2(:,7)));
  tmp2=tmp2(idx,:);
end
if(handles.interestingfeaturehistograms_absdprime)
  tmp2(:,7)=abs(tmp2(:,7));
end
tmp2=sortrows(tmp2,-7);

tmp=cell(size(tmp2,1),7);
tmp(:,1)=handles.grouplist(tmp2(:,1));
tmp(:,2)=cellstr(num2str(tmp2(:,2),'%-d'));
idx=~isnan(tmp2(:,3));
tmp(idx,3)=handles.grouplist(tmp2(idx,3));
tmp(:,4)=cellstr(num2str(tmp2(:,4),'%-d'));
tmp(:,5)=handles.behaviorlist(tmp2(:,5));
tmp(:,6)=handles.featurelist(tmp2(:,6));
tmp(:,7)=num2cell(tmp2(:,7));
set(handles.Table,'Data',tmp);
set(handles.Table,'ColumnName',{'Group' 'n' 'Group2' 'n2' 'Behavior' 'Feature' 'd'''});
set(handles.Table,'ColumnWidth',{75 50 75 50 150 100 75});
set(handles.Table,'RowStriping','on','BackgroundColor',[1 1 1; 0.95 0.95 0.95]);

handles.table_data=tmp2;
handles.table='histogram';
guidata(hObject,handles);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- 
function data=calculate_entiretimeseries(behavior_data,feature_data,sexdata,individual,xoffset);

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
  feature_data.data=cellfun(@(x) strcmp(x,vals{1}),feature_data.data,'uniformoutput',false);
end

k=1;
if(xoffset==1)
  data=nan(length(feature_data.data),max(behavior_data.allScores.tEnd));
else
  data=nan(length(feature_data.data),max(behavior_data.allScores.tEnd)-min(behavior_data.allScores.tStart));
end
for i=1:length(feature_data.data)  % individual
  if((~isnan(individual)) && (i~=individual))  continue;  end
  if(sum(sexdata{i}) < length(sexdata{i})/2)  continue;  end
  switch(xoffset)
    case(1)
      behavior_data.allScores.tStart(i);
    case(2)
      1;
    case(3)
      behavior_data.allScores.tStart(i)-min(behavior_data.allScores.tStart)+1;
  end
  data(k,ans:(ans+length(feature_data.data{i})-1))=feature_data.data{i};
  k=k+1;
end


% --- 
function triggered_data=calculate_triggeredtimeseries(behavior_data,behavior_logic,behavior_data2,feature_data,...
    sexdata,individual,timing,windowradius,subtractmean)

if(iscell(feature_data.data{1}))
  vals=unique([feature_data.data{:}]);
  if(length(vals)>2)  error('uhoh');  end
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


% --- Executes on button press in FeatureTimeSeries.
function FeatureTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to FeatureTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user experiment (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals(:,1))'];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

ggee=1:length(handlesexperimentlist);
individual=handles.individualvalue;
if(individual>3)
  ggee=find(cumsum_num_indi_per_exp<(individual-3),1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='feature time series';

h=figure('toolbar','figure');  hold on;
guidata(h,handles);
uicontrol(h,'style','pushbutton','string','?','position',[10 10 20 20],'callback',@figure_info_callback);

bb=handles.behaviorvalue;
if(handles.featuretimeseries_timing==1)  bb=1;  end
if(bb==length(handles.behaviorlist))  bb=1:(bb-1);  end

behavior_logic=handles.behaviorlogic;
score_file2=handles.scorefiles{handles.behaviorvalue2};
feature_value=handles.featurevalue;
feature_list=handles.featurelist;
sexdata=handles.sexdata;
timing=handles.featuretimeseries_timing;
style=handles.featuretimeseries_style;
centraltendency=handles.prefs_centraltendency;
dispersion=handles.prefs_dispersion;
convolutionwidth=handles.prefs_convolutionwidth;
subtractmean=handles.featuretimeseries_subtractmean;
windowradius=handles.featuretimeseries_windowradius;
xoffset=handles.prefs_timeseriesxoffset;

h=[];
for b=bb

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    subplot(ceil(length(bb)/ans),ans,b);
  end
  hold on;

  score_file=handles.scorefiles{b};

  data=cell(1,length(handlesexperimentlist));
  parfor ge=ggee
  %for ge=ggee
    if((individual<4)&&(~ismember(ge,selected_exp)))  continue;  end

    behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
    if(behavior_logic>1)
      behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
    else
      behavior_data2=[];
    end
    feature_data=load(fullfile(handlesexperimentlist{ge},'perframe',...
        [feature_list{feature_value} '.mat']));

    tmp2=sexdata{ge};
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
      calculate_entiretimeseries(behavior_data,feature_data,tmp2,tmp,xoffset);
      conv(nanmean(ans,1),ones(1,convolutionwidth),'same');
      data{ge}=ans./conv(ones(1,length(ans)),ones(1,convolutionwidth),'same');
    else
      calculate_triggeredtimeseries(behavior_data,behavior_logic,behavior_data2,...
          feature_data,tmp2,tmp,timing,windowradius,subtractmean);
      data{ge}=nanmean(ans,1);
    end
  end

  idx=cellfun(@isempty,data);
  data=data(~idx);

  if(timing==1)
    max(cellfun(@(x) size(x,2),data));
    cellfun(@(x) [x nan(size(x,1),ans-size(x,2))],data,'uniformoutput',false);
    ydata=cat(1,ans{:});
    xdata=1:size(ydata,2);
  else
    ydata=cat(1,data{:});
    xdata=-windowradius:windowradius;
  end

  tstr='';
  if(timing>1)
    tstr=char(strrep(handles.behaviorlist(b),'_','-'));
    switch(handles.behaviorlogic)
      case 2
        tstr=[tstr ' AND '];
      case 3
        tstr=[tstr ' AND NOT '];
      case 4
        tstr=[tstr ' OR '];
    end
    if(behavior_logic>1)
      tstr=[tstr char(strrep(handles.behaviorlist(handles.behaviorvalue2),'_','-'))];
    end
  end
  time_base=xdata./handles.fps;
  xstr='time (sec)';
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (min)';
  end
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (hr)';
  end
  if(time_base(end)>120)
    time_base=time_base./24;
    xstr='time (d)';
  end
  units=load(fullfile(handlesexperimentlist{ggee(1)},'perframe',...
      [feature_list{feature_value} '.mat']),'units');
  ystr=get_label(feature_list(feature_value),units.units);

  print_csv_help(fid,handles.type,tstr,xstr,ystr);

  for g=1:length(handles.grouplist)
    color=handles.colors(g,:);

    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);

    if(individual<4)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end

    h(g)=plot_it(time_base,ydata(idx,:),style,centraltendency,dispersion,color,1,fid,handlesexperimentlist(idx));
  end

  xlabel(xstr);
  ylabel(ystr);
  title(tstr);
  axis tight;  zoom reset;

end
idx=find(h>0);
if(individual<4)
  legend(h(idx),handles.grouplist);
else
  legend(h(idx),handles.individuallist(individual+cumsum_num_indi_per_exp(ggee)));
end

fclose(fid);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


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

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

experiment_value=get(handles.ExperimentList,'Value');
experiment_list=get(handles.ExperimentList,'String');
experiment_value2=get(handles.ExperimentList2,'Value');
experiment_list2=get(handles.ExperimentList2,'String');
behavior_list=get(handles.BehaviorList,'String');
feature_list=get(handles.FeatureList,'String');
statistic=handles.prefsstat;
windowradius=handles.featuretimeseries_windowradius;

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

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


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
function [h,h2]=errorbarplot(x,b,dn,dp,color)

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
h2=errorbar(xb, b, dn, dp);
set(h2(1),'linewidth',1);            % This changes the thickness of the errorbars
set(h2(1),'color','k');              % This changes the color of the errorbars
set(h2(1),'linestyle','none');       % This removes the connecting


% --- Executes on button press in BehaviorBarChart.
function BehaviorBarChart_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorBarChart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals(:,1))'];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

ggee=1:length(handlesexperimentlist);
individual=handles.individualvalue;
if(individual>3)
  ggee=find(cumsum_num_indi_per_exp<(individual-3),1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='behavior bar chart';

h=figure('toolbar','figure');  hold on;
guidata(h,handles);
uicontrol(h,'style','pushbutton','string','?','position',[10 10 20 20],'callback',@figure_info_callback);

bb=handles.behaviorvalue;
if(bb==length(handles.behaviorlist))  bb=1:(bb-1);  end

behavior_logic=handles.behaviorlogic;
score_file2=handles.scorefiles{handles.behaviorvalue2};
sexdata=handles.sexdata;
perwhat=handles.behaviorbarchart_perwhat;

h=[];
table_data={};
for b=bb

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

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    subplot(ceil(length(bb)/ans),ans,b);
  end
  hold on;

  score_file=handles.scorefiles{b};

  collated_data=cell(1,length(ggee));
  parfor ge=ggee
  %for ge=ggee
    if((individual<4)&&(~ismember(ge,selected_exp)))  continue;  end

    behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
    if(behavior_logic>1)
      behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
    else
      behavior_data2=[];
    end

    frames_labelled=nan(1,length(behavior_data.allScores.t0s));
    frames_total=nan(1,length(behavior_data.allScores.t0s));
    sex=nan(1,length(behavior_data.allScores.t0s));
    for i=1:length(behavior_data.allScores.t0s)  % individual
      tmp1=false(1,length(sexdata{ge}{i}));
      for j=1:length(behavior_data.allScores.t0s{i})  % bout
        tmp1((behavior_data.allScores.t0s{i}(j):(behavior_data.allScores.t1s{i}(j)-1))...
            -behavior_data.allScores.tStart(i)+1)=true;
      end
      tmp2=[];
      if(behavior_logic>1)
        tmp2=false(1,length(sexdata{ge}{i}));
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
      sex(i)=sum(sexdata{ge}{i}(1:length(partition_idx))) > (length(partition_idx)/2);
      frames_labelled(i)=sum(partition_idx);
      frames_total(i)=length(partition_idx);
    end

    collated_data{ge}={frames_labelled frames_total sex};
  end

  idx=cellfun(@isempty,collated_data);
  collated_data=collated_data(~idx);

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

  exp_separators=[];  maxy=0;  k=[];  m=0;  table_data{end+1}=[];
  for g=1:length(handles.grouplist)
    color=handles.colors(g,:);

    if(individual<4)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end

    xticklabels{g}=handles.grouplist{g};

    switch(handles.behaviorbarchart_perwhat)
      case 1  % per group
        table_data{end}(g)=100*sum([frames_labelled{idx}])./sum([frames_total{idx}]);
        h(g)=bar(g,table_data{end}(g));
        set(h(g),'facecolor',color);
      case 2  % per experiment, error bars
        table_data{end}{g}=100*cellfun(@sum,frames_labelled(idx))./cellfun(@sum,frames_total(idx));
        [ct{g},dp{g},dn{g}]=...
            calculate_ct_d(table_data{end}{g},handles.prefs_centraltendency,handles.prefs_dispersion);
        h(g)=errorbarplot(g,ct{g},ct{g}-dn{g},dp{g}-ct{g},color);
      case 3  % per fly, grouped
        cumsum(cellfun(@length,frames_labelled(idx)))';
        exp_separators=[exp_separators; ans+sum(k)];
        table_data{end}{g}=100.*[frames_labelled{idx}]./[frames_total{idx}];
        maxy=max([maxy table_data{end}{g}]);
        h(g)=bar((1:length(table_data{end}{g}))+sum(k),table_data{end}{g},...
            'barwidth',1,'edgecolor','none');
        set(h(g),'facecolor',color);
        k(end+1)=length(table_data{end}{g});
        fprintf(fid,['%% data, %s\n'],xticklabels{g});
        fprintf(fid,'%g, ',[table_data{end}{g}]);
        fprintf(fid,'\n');
      case 4  % per fly, stern-style
        table_data{end}{g}=cell(1,length(frames_labelled(idx)));
        for e=idx
          table_data{end}{g}{e}=100.*frames_labelled{e}./frames_total{e};
          [ct,dp,dn]=calculate_ct_d(table_data{end}{g}{e},...
              handles.prefs_centraltendency,handles.prefs_dispersion);
          h(g)=plot(m,ct,'o','color',color);
          plot([m m],[dp dn],'-','color',color);
          plot(m+(1:length(table_data{end}{g}{e})),table_data{end}{g}{e},'.','color',color);
          m=m+16+length(table_data{end}{g}{e});
        end
        [ct,dp,dn]=calculate_ct_d([table_data{end}{g}{:}],...
            handles.prefs_centraltendency,handles.prefs_dispersion);
        plot(m,ct,'o','color',color,'markersize',9);
        plot([m m],[dp dn],'-','color',color,'linewidth',3);
        m=m+24;
        k(end+1)=24+16*length(table_data{end}{g})+length([table_data{end}{g}{:}]);
        fprintf(fid,['%% data, %s\n'],xticklabels{g});
        fprintf(fid,'%g, ',[table_data{end}{g}{:}]);
        fprintf(fid,'\n');
    end
  end

  switch(handles.behaviorbarchart_perwhat)
    case 1  % per group
      fprintf(fid,['%% xdata\n']);  fprintf(fid,'%s, ',xticklabels{:});  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, per group\n']);  fprintf(fid,'%g, ',table_data{end});  fprintf(fid,'\n');
    case 2  % per experiment, error bars
      fprintf(fid,['%% xdata\n']);  fprintf(fid,'%s, ',xticklabels{:});  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT+D\n']);  fprintf(fid,'%g, ',[dp{:}]);  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT-D\n']);  fprintf(fid,'%g, ',[dn{:}]);  fprintf(fid,'\n');
      fprintf(fid,['%% ydata, CT\n']);  fprintf(fid,'%g, ',[ct{:}]);  fprintf(fid,'\n');
    case 3  % per fly, grouped
      l=exp_separators(1:2:(end-1));
      r=exp_separators(2:2:end);
      hh=patch(0.5+[l r r l l]',repmat([0 0 maxy*1.05 maxy*1.05 0]',1,floor(length(exp_separators)/2)),...
          [0.95 0.95 0.95]);
      set(hh,'edgecolor','none');
      set(gca,'children',flipud(get(gca,'children')));
      k=round(cumsum(k)-k/2);
    case 4  % per fly, stern-style
      k=round(cumsum(k)-k/2);
  end

  if(isempty(k))  k=1:length(frames_labelled);  end
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
idx=find(h>0);
if(individual<4)
  legend(h(idx),handles.grouplist);
else
  legend(h(idx),handles.individuallist(individual+cumsum_num_indi_per_exp(ggee)));
end

if((ismember(handles.behaviorbarchart_perwhat,[2 3])) && (individual<4))
  tmp={'behavior'};
  for b=1:length(table_data)
    tmp{4*b+1,1}=handles.behaviorlist{bb(b)};
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

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


% --- Executes on button press in BehaviorTimeSeries.
function BehaviorTimeSeries_Callback(hObject, eventdata, handles)
% hObject    handle to BehaviorTimeSeries (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

handlesexperimentlist=[handles.experimentlist{:}];

cumsum_num_indi_per_exp=[0 cumsum(handles.individuals(:,1))'];
cumsum_num_exp_per_group=[0 cumsum(cellfun(@length,handles.experimentlist))];
mat2cell(cumsum_num_exp_per_group(1:end-1),1,ones(1,length(cumsum_num_exp_per_group)-1));
cellfun(@(x,y) x+y,handles.experimentvalue,ans,'uniformoutput',false);
selected_exp=[ans{:}];
cumsum_num_selexp_per_group=[0 cumsum(cellfun(@length,handles.experimentvalue))];

ggee=1:length(handlesexperimentlist);
individual=handles.individualvalue;
if(individual>3)
  ggee=find(cumsum_num_indi_per_exp<(individual-3),1,'last');
  individual=individual-cumsum_num_indi_per_exp(ggee);
end

fid=fopen('most_recent_figure.csv','w');

handles.type='behavior time series';

h=figure('toolbar','figure');  hold on;
guidata(h,handles);
uicontrol(h,'style','pushbutton','string','?','position',[10 10 20 20],'callback',@figure_info_callback);

bb=handles.behaviorvalue;
if(bb==length(handles.behaviorlist))  bb=1:(bb-1);  end

behavior_logic=handles.behaviorlogic;
score_file2=handles.scorefiles{handles.behaviorvalue2};
sexdata=handles.sexdata;
convolutionwidth=handles.prefs_convolutionwidth;
style=handles.behaviortimeseries_style;
centraltendency=handles.prefs_centraltendency;
dispersion=handles.prefs_dispersion;
xoffset=handles.prefs_timeseriesxoffset;

h=[];
for b=bb

  if(length(bb)>1)
    ceil(sqrt(length(bb)));
    subplot(ceil(length(bb)/ans),ans,b);
  end
  hold on;

  score_file=handles.scorefiles{b};

  behavior_cumulative=cell(length(ggee),1);
  parfor gei=1:numel(ggee),
  %for gei=1:numel(ggee),
    
    ge = ggee(gei);
    if((individual<4)&&(~ismember(ge,selected_exp)))  continue;  end

    behavior_data=load(fullfile(handlesexperimentlist{ge},score_file));
    behavior_data2=[];
    if(behavior_logic>1)
      behavior_data2=load(fullfile(handlesexperimentlist{ge},score_file2));
    end

    if(xoffset==1)
      parfor_tmp=zeros(1,max(behavior_data.allScores.tEnd));
    else
      parfor_tmp=zeros(1,max(behavior_data.allScores.tEnd)-min(behavior_data.allScores.tStart));
    end

    k=0;
    for i=1:length(behavior_data.allScores.t0s)   % individual
      if((individual>3)&&((individual-3)~=i))  continue;  end

      tmp1=zeros(1,max(behavior_data.allScores.tEnd));
      tmp1(behavior_data.allScores.t0s{i})=1;
      tmp1(behavior_data.allScores.t1s{i})=-1;
      tmp1=logical(cumsum(tmp1));

      tmp2=[];
      if(behavior_logic>1)
        tmp2=zeros(1,max(behavior_data.allScores.tEnd));
        tmp2(behavior_data2.allScores.t0s{i})=1;
        tmp2(behavior_data2.allScores.t1s{i})=-1;
        tmp2=logical(cumsum(tmp2));
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
      switch(individual)
        case(2)
          partition_idx = partition_idx & sexdata{ge}{i}(1:length(partition_idx));
        case(3)
          partition_idx = partition_idx & (~sexdata{ge}{i}(1:length(partition_idx)));
      end
      
      idx=[];
      switch(xoffset)
        case(1)
          idx = find(partition_idx);
        case(2)
          idx = find(partition_idx) - behavior_data.allScores.tStart(i);
        case(3)
          idx = find(partition_idx) - min(behavior_data.allScores.tStart);
      end
      parfor_tmp(idx)=parfor_tmp(idx)+1;
      k=k+1;
    end
    parfor_tmp./k;
    behavior_cumulative{gei}=conv(ans,ones(1,convolutionwidth),'same')...
        ./conv(ones(1,length(ans)),ones(1,convolutionwidth),'same');
  end
  parfor_tmp_len = min( cellfun(@numel, behavior_cumulative) );
  %parfor_tmp_len = max( cellfun(@numel, behavior_cumulative) );
  for gei = 1:numel(ggee),
    behavior_cumulative{gei} = behavior_cumulative{gei}(1:parfor_tmp_len);
       %[behavior_cumulative{gei},nan(1,parfor_tmp_len-numel(behavior_cumulative{gei}))];
  end
  behavior_cumulative = cell2mat(behavior_cumulative);

  tstr='';
  time_base=(1:size(behavior_cumulative,2))./handles.fps;
  xstr='time (sec)';
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (min)';
  end
  if(time_base(end)>300)
    time_base=time_base./60;
    xstr='time (hr)';
  end
  if(time_base(end)>120)
    time_base=time_base./24;
    xstr='time (d)';
  end
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

  print_csv_help(fid,handles.type,tstr,xstr,ystr);

  for g=1:length(handles.grouplist)
    color=handles.colors(g,:);

    fprintf(fid,['%% group ' handles.grouplist{g} '\n']);

    if(individual<4)
      idx=(cumsum_num_selexp_per_group(g)+1):(cumsum_num_selexp_per_group(g+1));
    else
      find(cumsum_num_exp_per_group<ggee,1,'last');
      if(ans~=g)  continue;  end
      idx=1;
    end

    h(g)=plot_it(time_base,100.*behavior_cumulative(idx,:),...
        style,centraltendency,dispersion,color,1,fid,handlesexperimentlist(idx));
  end

  xlabel(xstr);
  ylabel(ystr);
  axis tight;  zoom reset;

end
idx=find(h>0);
if(individual<4)
  legend(h(idx),handles.grouplist);
else
  legend(h(idx),handles.individuallist(individual+cumsum_num_indi_per_exp(ggee)));
end

fclose(fid);

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


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

set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

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

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


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
ret_val=[ans(1:(end-1)) ' (' num2str(arg{2},3) '%)']; %#ok<COLND>



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


set(handles.Status,'string','Thinking...','foregroundcolor','b');
set(handles.figure1,'pointer','watch');
drawnow;

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

set(handles.Status,'string','Ready.','foregroundcolor','g');
set(handles.figure1,'pointer','arrow');
drawnow;


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
%    %handles.individualvalue=floor(eventdata.Indices(r,2)/2);
%    %if(((color=='r') && ...
%    %    (handles.individualvalue>(3+sum(handles.individuals(1:length(handles.experimentlist)))))) || ...
%    %   ((color=='b') && ...
%    %    (handles.individualvalue>(3+sum(handles.individuals((end-length(handles.experimentlist2)+1):end))))))
%    %  return;
%    %end
%    %if((color=='b') && (handles.individualvalue>3))
%    %  handles.individualvalue=handles.individualvalue+sum(handles.individuals(1:length(handles.experimentlist)));
%    %end
%    %set(handles.IndividualList,'Value',handles.individualvalue);

    data=handles.raw_table_data{eventdata.Indices(end,1),eventdata.Indices(end,2)};
    tmp=unique(handles.raw_table_data{eventdata.Indices(end,1),2});
    hist(data,min(tmp):max(tmp));
%    %plot(x,n./sum(n),[color '-']);

%    %tmp{r,1}=sprintf('%10.3g ',mean(data));
%    %tmp{r,2}=sprintf('%10.3g ',std(data));
%    %tmp{r,3}=sprintf('%10.3g ',std(data)./sqrt(length(data)));
%    %tmp{r,4}=sprintf('%10.3g ',median(data));
%    %tmp{r,5}=sprintf('%10.3g ',prctile(data,25));
%    %tmp{r,6}=sprintf('%10.3g ',prctile(data,75));
%    %tmp{r,7}=sprintf('%10.3g ',prctile(data,5));
%    %tmp{r,8}=sprintf('%10.3g ',prctile(data,95));
%  %end

%  %{'Mean' 'Std.Dev.' 'Std.Err.' 'Median' '25%' '75%' '5%' '95%'};
%  %tmp2(1,:)=sprintf('%10s ',ans{:});
%  %for i=1:size(tmp,1)
%  %  tmp2(i+1,:)=[tmp{i,:}];
%  %end
%  %v=axis;
%  %h=text(v(1),v(4),tmp2,'color',[0 0.5 0],'tag','stats','verticalalignment','top','fontname','fixed');
%  %if(handles.stats)
%  %  set(h,'visible','on');
%  %else
%  %  set(h,'visible','off');
%  %end
  
  xlabel('closest fly (#)');
  axis tight;

%elseif(strcmp(handles.table,'timeseries') || strcmp(handles.table,'histogram'))
%  set(handles.BehaviorList,'Value',handles.behaviorvalue);
%  set(handles.BehaviorLogic,'Value',1);
%  set(handles.BehaviorList2,'enable','off');
%  set(handles.FeatureList,'Value',handles.featurevalue);
%  set(handles.IndividualList,'Value',handles.individual);

%  if(strcmp(handles.table,'timeseries'))
%    handles.timeseries_timing=handles.table_data(eventdata.Indices(end,1),3)+1;
%    FeatureTimeSeries_Callback(hObject, eventdata, handles);

elseif(strcmp(handles.table,'histogram'))
  switch(eventdata.Indices(end,2))
    case {1,3}
      group=handles.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2));
      if(isnan(group))
        questdlg('Remove all during / not-during comparisons from table?','','Yes','No','No');
      else
        questdlg(['Remove all ' handles.grouplist{group} ' groups from table?'],'','Yes','No','No');
      end
      if(strcmp(ans,'No'))  return;  end
      tmp=get(handles.Table,'Data');
      if(isnan(group))
        idx=find((~isnan(handles.table_data(:,1)))&(~isnan(handles.table_data(:,3))));
      else
        idx=find((handles.table_data(:,1)~=group)&(handles.table_data(:,3)~=group));
      end
      set(handles.Table,'Data',tmp(idx,:));
      handles.table_data=handles.table_data(idx,:);
    case {2,4}
      thresh=handles.table_data(eventdata.Indices(end,1),eventdata.Indices(end,2));
      questdlg(['Remove all rows for which n or n2 is less than or equal to ' num2str(thresh) ...
          ' from table?'],'','Yes','No','No');
      if(strcmp(ans,'No'))  return;  end
      tmp=get(handles.Table,'Data');
      idx=find((handles.table_data(:,2)>thresh)&(handles.table_data(:,4)>thresh));
      set(handles.Table,'Data',tmp(idx,:));
      handles.table_data=handles.table_data(idx,:);
    case {5}
      questdlg(['Remove all ' handles.behaviorlist{handles.table_data(eventdata.Indices(end,1),5)} ...
          ' behaviors from table?'],'','Yes','No','No');
      if(strcmp(ans,'No'))  return;  end
      tmp=get(handles.Table,'Data');
      idx=find(handles.table_data(:,5)~=handles.table_data(eventdata.Indices(end,1),5));
      set(handles.Table,'Data',tmp(idx,:));
      handles.table_data=handles.table_data(idx,:);
    case {6}
      questdlg(['Remove all ' handles.featurelist{handles.table_data(eventdata.Indices(end,1),6)} ...
          ' features from table?'],'','Yes','No','No');
      if(strcmp(ans,'No'))  return;  end
      tmp=get(handles.Table,'Data');
      idx=find(handles.table_data(:,6)~=handles.table_data(eventdata.Indices(end,1),6));
      set(handles.Table,'Data',tmp(idx,:));
      handles.table_data=handles.table_data(idx,:);
    case {7}
      handles.behaviorvalue=handles.table_data(eventdata.Indices(end,1),5);
      handles.behaviorlogic=1;
      handles.featurevalue=handles.table_data(eventdata.Indices(end,1),6);
      handles.individual=1;
      handles.featurehistogram_notduring=isnan(handles.table_data(eventdata.Indices(end,1),3));
      menu_featurehistogram_notduring_set(handles.featurehistogram_notduring);
      FeatureHistogram_Callback(hObject, eventdata, handles);
  end
  update_figure(handles);

end

guidata(hObject,handles);


% ---
function menu_classify_forcecompute_set(arg)

handles=guidata(gcf);
if(arg)
  set(handles.MenuClassifyForceCompute,'Checked','on');
else
  set(handles.MenuClassifyForceCompute,'Checked','off');
end


% --------------------------------------------------------------------
function MenuClassifyForceCompute_Callback(hObject, eventdata, handles)
% hObject    handle to MenuClassifyForceCompute (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.classify_forcecompute=~handles.classify_forcecompute;
menu_classify_forcecompute_set(handles.classify_forcecompute);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuClassify_Callback(hObject, eventdata, handles)
% hObject    handle to MenuClassify (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
function MenuFeatureTimeSeriesWindowRadius_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureTimeSeriesWindowRadius (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.featuretimeseries_windowradius=str2num(char(inputdlg({'Window radius:'},'',1,...
    {num2str(handles.featuretimeseries_windowradius)})));
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
set(handles.MenuFeatureHistogramStdDevPerBout,'Checked','off');
switch(arg)
  case(1), set(handles.MenuFeatureHistogramPerFrame,'Checked','on');
  case(2), set(handles.MenuFeatureHistogramMeanPerBout,'Checked','on');
  case(3), set(handles.MenuFeatureHistogramMedianPerBout,'Checked','on');
  case(4), set(handles.MenuFeatureHistogramMaxPerBout,'Checked','on');
  case(5), set(handles.MenuFeatureHistogramMinPerBout,'Checked','on');
  case(6), set(handles.MenuFeatureHistogramStdDevPerBout,'Checked','on');
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
function MenuFeatureHistogramStdDevPerBout_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFeatureHistogramStdDevPerBout (see GCBO)
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
%handles.interestingfeaturehistograms_cache=[];
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
%handles.interestingfeaturehistograms_cache=[];
guidata(hObject,handles);


% ---
function menu_interestingfeaturehistograms_absdprime_set(arg)

handles=guidata(gcf);
if(arg)
  set(handles.MenuInterestingFeatureHistogramsAbsDPrime,'Checked','on');
else
  set(handles.MenuInterestingFeatureHistogramsAbsDPrime,'Checked','off');
end


% --------------------------------------------------------------------
function MenuInterestingFeatureHistogramsAbsDPrime_Callback(hObject, eventdata, handles)
% hObject    handle to MenuInterestingFeatureHistogramsAbsDPrime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.interestingfeaturehistograms_absdprime=~handles.interestingfeaturehistograms_absdprime;
menu_interestingfeaturehistograms_absdprime_set(handles.interestingfeaturehistograms_absdprime);
%handles.interestingfeaturehistograms_cache=[];
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
function MenuPrefsTimeSeriesXOffset_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsTimeSeriesXOffset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% ---
function menu_prefstimeseriesxoffset_set(arg)

handles=guidata(gcf);

set(handles.MenuPrefsTimeSeriesXOffsetNone,'Checked','off');
set(handles.MenuPrefsTimeSeriesXOffsetStart,'Checked','off');
set(handles.MenuPrefsTimeSeriesXOffsetMinStart,'Checked','off');
switch(arg)
  case(1), set(handles.MenuPrefsTimeSeriesXOffsetNone,'Checked','on');
  case(2), set(handles.MenuPrefsTimeSeriesXOffsetStart,'Checked','on');
  case(3), set(handles.MenuPrefsTimeSeriesXOffsetMinStart,'Checked','on');
end


% --------------------------------------------------------------------
function MenuPrefsTimeSeriesXOffsetNone_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsTimeSeriesXOffsetNone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_timeseriesxoffset=1;
menu_prefstimeseriesxoffset_set(handles.prefs_timeseriesxoffset);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuPrefsTimeSeriesXOffsetStart_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsTimeSeriesXOffsetStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_timeseriesxoffset=2;
menu_prefstimeseriesxoffset_set(handles.prefs_timeseriesxoffset);
guidata(hObject,handles);


% --------------------------------------------------------------------
function MenuPrefsTimeSeriesXOffsetMinStart_Callback(hObject, eventdata, handles)
% hObject    handle to MenuPrefsTimeSeriesXOffsetMinStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.prefs_timeseriesxoffset=3;
menu_prefstimeseriesxoffset_set(handles.prefs_timeseriesxoffset);
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
function MenuFileUpdate_Callback(hObject, eventdata, handles)
% hObject    handle to MenuFileUpdate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'pointer','watch');
handles=update_experiment_data(handles,true,true,true);
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
function figure_info_callback(src,evt)

handles=guidata(src);

CT={'Mean' 'Median' 'Mode'};
D={'Std. Dev.' 'Std. Err.' '5%-95%' '25%-75%'};
style={'Central Tendency' 'Central Tendency & Dispersion' 'Overlayed per-Exp Means'};
behaviorbarchart_perwhat={'per Group' 'per Experiment' 'per Fly' 'per Fly'};
featurehistogram_perwhat={'per Frame' 'Mean per Bout' 'Median per Bout' 'Max per Bout' 'Min per Bout' 'Std. Dev. per Bout'};
featuretimeseries_timing={'Entire Recording' 'Onset Triggered' 'Offset Triggered'};
timeseriesxoffset={'none', 'start', 'min(start)'};

tmp={};

switch handles.type
  case 'behavior bar chart'
    tmp{end+1}=['perwhat = ' behaviorbarchart_perwhat{handles.behaviorbarchart_perwhat}];
  case 'behavior time series'
    tmp{end+1}=['style = ' style{handles.behaviortimeseries_style}];
  case 'feature histogram'
    tmp{end+1}=['perwhat = ' featurehistogram_perwhat{handles.featurehistogram_perwhat}];
    tmp{end+1}=['style = ' style{handles.featurehistogram_style}];
    tmp{end+1}=['logbinsize=' num2str(handles.featurehistogram_logbinsize)];
    tmp{end+1}=['notduring=' num2str(handles.featurehistogram_notduring)];
    tmp{end+1}=['nbins=' num2str(handles.featurehistogram_nbins)];
  case 'feature time series'
    tmp{end+1}=['timing = ' featuretimeseries_timing{handles.featuretimeseries_timing}];
    tmp{end+1}=['style = ' style{handles.featuretimeseries_style}];
    tmp{end+1}=['subtractmean=' num2str(handles.featuretimeseries_subtractmean)];
    tmp{end+1}=['windowradius=' num2str(handles.featuretimeseries_windowradius)];
end
tmp{end+1}='';

tmp{end+1}=['central tendency = ' CT{handles.prefs_centraltendency}];
tmp{end+1}=['dispersion = '  D{handles.prefs_dispersion}];
tmp{end+1}=['convolution width = '  num2str(handles.prefs_convolutionwidth) ' frames'];
tmp{end+1}=['time series x-offset = '  timeseriesxoffset{handles.prefs_timeseriesxoffset}];
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



%BOUT STATISTICS -- not converted to multivariate yet
%
%Create a table of bout lengths (BL) and inter-bout lengths (IBL) using the
%"Bout Stats" button.  As for behavior stats, each row is a behavior;  the
%columns break it down by sex and individual;  the experiment groups are
%striped;  logical operators can be applied; and the Prefs contextual menu
%controls statistics.
%
%Selecting a cell plots a histogram, normalized to unit area, of bout and
%inter-bout lengths.  The "LogY" and "Stats" buttons can be used to scale
%the axis and overlay summary statistics, respectively.  Multiple selected
%cells are overlayed.
%
%
%SOCIAL STATISTICS -- not converted to multivariate yet
%
%Create a table of the closest individual using the "Social Stats" button.
%Choose which behavior to analysis using the pull-down menu in the Behavior
%panel, and which metric to use for closest using the pull-down menu in the
%Feature panel.  In the resulting table, each row is an individual, and the
%columns show the mode and percentage of frames of the closest individual
%for each bout, as well as the overall mode of the per-bout modes.
%
%Selecting a cell in the table plots a histogram of the closest indidividual.


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

html_file = 'http://jaaba.sourceforge.net/PlottingResults.html';
if isdeployed,
  [stat,msg] = myweb_nocheck(html_file);
  if stat ~= 0,
    errordlg({'Please see documentation at http://jaaba.sourceforge.net'
      'Error opening webpage within MATLAB:'
      msg});
  end
else
  web('-browser',html_file);
end
