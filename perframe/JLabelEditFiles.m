function varargout = JLabelEditFiles(varargin)
% JLABELEDITFILES MATLAB code for JLabelEditFiles.fig
%      JLABELEDITFILES, by itself, creates a new JLABELEDITFILES or raises the existing
%      singleton*.
%
%      H = JLABELEDITFILES returns the handle to a new JLABELEDITFILES or the handle to
%      the existing singleton*.
%
%      JLABELEDITFILES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JLABELEDITFILES.M with the given input arguments.
%
%      JLABELEDITFILES('Property','Value',...) creates a new JLABELEDITFILES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JLabelEditFiles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JLabelEditFiles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JLabelEditFiles

% Last Modified by GUIDE v2.5 09-Jul-2012 13:42:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JLabelEditFiles_OpeningFcn, ...
                   'gui_OutputFcn',  @JLabelEditFiles_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}) && exist(varargin{1}), %#ok<EXIST>
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before JLabelEditFiles is made visible.
function JLabelEditFiles_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JLabelEditFiles (see VARARGIN)

% have not pushed the done button yet
handles.success = false;

%
[behaviorParamsFile, ...
  disableBehavior, ...
  JLabelHandle,...
  handles.JLabelSplashHandle] = ...
  myparse(varargin,...
  'behaviorParamsFile',fullfile('params','BehaviorList.xml'),...
  'disableBehavior',false,...
  'JLabelHandle',[],...
  'JLabelSplashHandle',[]);

handles.behaviorXMLeditFile = 'xml_edit_tools/behaviorConfigDefaults.xml';
set(handles.popupmode,'String',{'Normal','Advanced','Ground Truthing','Ground Truthing Advanced'});
handles.disableBehavior = disableBehavior;

if disableBehavior,

  handles.success = true;
  % Simply edit files, not editing the behavior part
  handles.JLabelHandle = JLabelHandle;
  handles.data = JLabelHandle.guidata.data;
  if handles.data.IsGTMode()
    if handles.data.IsAdvancedMode()
      set(handles.popupmode,'Value',find(strcmp(get(handles.popupmode,'String'),'Ground Truthing Advanced')));
    else
      set(handles.popupmode,'Value',find(strcmp(get(handles.popupmode,'String'),'Ground Truthing Advanced')));
    end
  elseif handles.data.IsAdvancedMode()
    set(handles.popupmode,'Value',find(strcmp(get(handles.popupmode,'String'),'Advanced')));
  else
    set(handles.popupmode,'Value',find(strcmp(get(handles.popupmode,'String'),'Normal')));
  end
  DisableBehaviorGui(handles);
  InitExperimentsGui(hObject,handles,JLabelHandle.guidata.data);
  handles.needJLabelInit = false;
  boxpos = get(handles.uipanelBehavior,'Position');
  guipos = get(hObject,'Position');
  set(hObject,'Position',[guipos(1:3) boxpos(2)-5]);

  guidata(hObject,handles);
  
else
  
  handles.JLabelHandle = JLabelHandle;
  handles.needJLabelInit = true;
  handles.defaultConfig = GetDefaultConfig();
  handles.projmanager = ProjectManager(behaviorParamsFile);

  set(handles.listbox_behavior,'String',handles.projmanager.GetProjectList());
  
  % If the user supplied a configuration file.
  if ~isempty(handles.JLabelHandle.guidata.configfilename)
    projNum = handles.projmanager.FindProjFromConfigfile(handles.JLabelHandle.guidata.configfilename);
    if ~isempty(projNum),
      set(handles.listbox_behavior,'Value',projNum);
      handles.projmanager.SetCurrentProject(projNum);
    else
      while true
        dlgStr = sprintf('%s\n%s\n%s','A project has not been defined for the configfile',...
          'To add a new project for this configfile enter the name below.',... 
          'To cancel, close this window');
        in = inputdlg(dlgStr);
        if numel(in)<1 || isempty(in{1}),
          break;
        end
        newName = in{1};
        if handles.projmanager.checkExist(newName),
          uiwait(warndlg('A project already exists with that name'));
        else
          handles.projmanager.AddProject(newName,handles.JLabelHandle.guidata.configfilename);
          handles = guidata(hObject);
          break
        end
      end
    end
  end
  
  projNum = handles.projmanager.GetCurrentProject();
  if ~isempty(projNum)
    set(handles.listbox_behavior,'Value',projNum);
  end
  guidata(hObject,handles);
  updateConfigParams(handles);
end


% Add color for mac's.
buttonNames = {'pushbutton_add','pushbutton_remove','pushbutton_load',...
              'pushbutton_loadwoexp','pushbutton_cancel','pushbutton_done',...
              'pushbutton_addlist'};
for buttonNum = 1:length(buttonNames)
  SetButtonImage(handles.(buttonNames{buttonNum}));
end

if ~isempty(handles.JLabelSplashHandle) && ishandle(handles.JLabelSplashHandle),
  delete(handles.JLabelSplashHandle);
end

% UIWAIT makes JLabelEditFiles wait for user response (see UIRESUME)
uiwait(handles.figure_JLabelEditFiles);


function DisableBehaviorGui(handles)
set(handles.listbox_behavior,'enable','off');
set(handles.config_table,'enable','off');
set(handles.pushbutton_addbehavior,'enable','off');
set(handles.pushbutton_removebehavior,'enable','off');
set(handles.pushbutton_chooseperframe,'enable','off');

function defaultConfig = GetDefaultConfig(name)
if nargin<1
  name = 'default';
end
defaultConfig.targets = struct('type','fly');
defaultConfig.behaviors = struct('names',name,'labelcolors',[0.7,0,0,0,0,0.7],'unknowncolor',[0,0,0]);
defaultConfig.file = struct('moviefilename','movie.ufmf',...
      'trxfilename','registered_trx.mat',...
      'labelfilename',sprintf('labeled%s.mat',name),...
      'gt_labelfilename',sprintf('gt_labeled%s.mat',name),...
      'scorefilename',sprintf('scores_%s.mat',name),...
      'perframedir','perframe',...
      'windowfilename','windowfeatures.mat',...
      'rootoutputdir','',...
      'featureparamfilename','',...
      'featureconfigfile',fullfile('params','featureConfig.xml'));
defaultConfig.plot.trx = struct('colormap','jet',...
	'colormap_multiplier','.7');
defaultConfig.plot.labels = struct('colormap','line',...
	'linewidth','3');
defaultConfig.perframe.params = struct('fov',3.1416,'nbodylengths_near',2.5,...
  'thetafil',[0.0625,0.25,0.375,0.25,0.0625]);
defaultConfig.perframe.landmark_params = struct(...
  'arena_center_mm_x',0,'arena_center_mm_y',0,'arena_radius_mm',60,'arena_type','circle');


function InitJLabelGui(handles)

if ~handles.needJLabelInit, return; end
handles.needJLabelInit = false;

% For behavior --
DisableBehaviorGui(handles);
handles.projmanager.WriteProjectManager();
handles.projmanager.SaveAllConfig();
handles.success = true;
% End for behavior --

% Initializes the JLabel gui once the user selects the behavior.
JLabelHandle = handles.JLabelHandle;
JLabelHandle.guidata.configparams = handles.projmanager.GetProjConfig();
JLabelHandle.guidata.configfilename = handles.projmanager.GetCurProjConfigfile();
JLabelHandle = JLabel('GetGUIPositions',JLabelHandle);
JLabelHandle = JLabel('InitializeState',JLabelHandle);
JLabelHandle = JLabel('InitializePlots',JLabelHandle);
handles.data = JLabelHandle.guidata.data;
SetLabelingMode(handles);
handles.JLabelHandle = JLabelHandle;
guidata(handles.figure_JLabelEditFiles,handles);
guidata(JLabelHandle.figure_JLabel,JLabelHandle);
InitExperimentsGui(handles.figure_JLabelEditFiles,handles,handles.data);


function InitExperimentsGui(hObject,handles,varargin)
%
% Remove behavior related parts.

% parse inputs, set defaults
%
if isempty(varargin),
  error('Usage: JLabelEditFiles(configfilename,...)');
end

% if ischar(varargin{1}),
%   configfilename = varargin{1};
%   handles.data = JLabelData(configfilename,varargin{2:end});
% else
%   handles.data = varargin{1};
% end

handles.data.SetStatusFn( @(s) SetStatusEditFiles(handles.figure_JLabelEditFiles,s));
handles.data.SetClearStatusFn( @() ClearStatusEditFiles(handles.figure_JLabelEditFiles));
% initialize listbox
set(handles.editClassifier,'String',handles.data.classifierfilename);
set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',numel(handles.data.expdirs));
set(handles.statusMsg,'String','');

if handles.data.IsGTMode(),
  if handles.data.IsAdvancedMode()
    popupNdx = find(strcmp(get(handles.popupmode,'String'),'Ground Truthing Advanced'));
  else
    popupNdx = find(strcmp(get(handles.popupmode,'String'),'Ground Truthing'));
  end

elseif handles.data.IsAdvancedMode(),
  popupNdx = find(strcmp(get(handles.popupmode,'String'),'Advanced'));
else
  popupNdx = find(strcmp(get(handles.popupmode,'String'),'Normal'));
end
set(handles.popupmode,'Value',popupNdx);
if handles.data.IsModeSet(),
  set(handles.popupmode,'Enable','off');
end

if handles.data.nexps>0
  set(handles.pushbutton_load,'Enable','off');
end

% height of a row in pixels
handles.table_row_height_px = 22;

% sizes for generate buttons
handles.generate_button_border_y_px = 1;
handles.generate_button_border_x_px = 5;
handles.generate_button_width_px = 100;
handles.generate_button_height_px = handles.table_row_height_px - handles.generate_button_border_y_px;

% buttons for generating the files
pos_table = get(handles.uitable_status,'Position');
top_table = pos_table(2)+pos_table(4);
right_table = pos_table(1)+pos_table(3);
for i = 1:numel(handles.data.filetypes),
  if JLabelData.CanGenerateFile(handles.data.filetypes{i}),
    pos = [right_table + handles.generate_button_border_x_px,...
      top_table - (i-1)*handles.table_row_height_px - ...
      handles.generate_button_height_px - handles.generate_button_border_y_px/2,...
      handles.generate_button_width_px,...
      handles.generate_button_height_px];
    buttonColor = get(handles.figure_JLabelEditFiles,'Color');
    handles.pushbutton_generate = uicontrol('Style','pushbutton','Units','Pixels',...
      'Parent',handles.figure_JLabelEditFiles,...
      'String','Generate','Position',pos,...
      'FontUnits','pixels','FontSize',14,...
      'tag',sprintf('pushbutton_generate_%d',i),...
      'Callback',@(hObject,eventdata) JLabelEditFiles('pushbutton_generate_Callback',hObject,eventdata,guidata(hObject),i),...
      'Backgroundcolor',buttonColor);
  end
end

% initialize status table
UpdateStatusTable(handles);

if handles.data.nexps>0,
  set(handles.pushbutton_load,'enable','off');
end

% Update handles structure
guidata(hObject, handles);


function UpdateStatusTable(handles)

v = get(handles.listbox_experiment,'Value');
if isempty(v) || v <= 0 || handles.data.nexps == 0 ,
  set(handles.uitable_status,'Data',{},'Visible','off');
  return;
end
if numel(v) > 1,
  v = v(end);
end

nfiles = numel(handles.data.filetypes);
data = cell([nfiles,2]);
data(:,1) = handles.data.filetypes;
for i = 1:nfiles,
  file = handles.data.filetypes{i};
  [file_exists,timestamp] = handles.data.FileExists(file,v);
  if file_exists,
    timestamp = datestr(timestamp);%,'yyyymmddTHHMMSS');
  end
  if JLabelData.IsRequiredFile(file),
    if file_exists,
      data{i,2} = sprintf('<html><font color="green">%s</font></html>',timestamp);
    else
      data{i,2} = '<html><font color="red">Missing</font></html>';
    end
  else
    if file_exists,
      data{i,2} = timestamp;
    else
      data{i,2} = 'Absent';
    end
  end
end
tableSize = get(handles.uitable_status,'Position');
set(handles.uitable_status,'Data',data,'Visible','on','ColumnWidth',{150 tableSize(3)-155});

% --- Outputs from this function are returned to the command line.
function varargout = JLabelEditFiles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.JLabelHandle;
varargout{2} = handles.success;
delete(hObject);


% --- Executes on selection change in listbox_experiment.
function listbox_experiment_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to listbox_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_experiment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_experiment
UpdateStatusTable(handles);

% --- Executes during object creation, after setting all properties.
function listbox_experiment_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to listbox_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
    SetButtonImage(hObject);
end

% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


oldv = get(handles.listbox_experiment,'Value');

if ~handles.disableBehavior && isempty(handles.projmanager.GetCurrentProject()), 
  uiwait(warndlg('Select a project before adding an experiment'));
  return
end

InitJLabelGui(handles);
handles = guidata(hObject);

% ask user for experiment directory
allexpdirs = uipickfiles('FilterSpec',handles.data.defaultpath,'Prompt','Add experiment directory');
if numel(allexpdirs)==1 && ~ischar(allexpdirs{1}),
  return;
end

set(handles.pushbutton_cancel,'enable','off');
SetLabelingMode(handles);

for ndx = 1:numel(allexpdirs)
  expdir = allexpdirs{ndx};
  if ismember(expdir,handles.data.expdirs),
    uiwait(warndlg(sprintf('Experiment directory %s already added',expdir),'Already added'));
    return;
  end
  
  
  [success,msg] = handles.data.AddExpDir(expdir);
  if ~success,
    uiwait(warndlg(sprintf('Error adding expdir %s: %s',expdir,msg)));
    return;
  end
  
  set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);

end

% update status table
UpdateStatusTable(handles);


function pushbutton_generate_Callback(hObject, eventdata, handles, row)

file = handles.data.filetypes{row};
if ~JLabelData.CanGenerateFile(file),
  return;
end
v = get(handles.listbox_experiment,'Value');
if isempty(v),
  return;
end
if numel(v) > 1,
  v = v(end);
end
expname = handles.data.expnames{v};
switch file,
%   case 'window',
%     [success,msg] = handles.data.GenerateWindowFeaturesFiles(v,true);
%     if ~success,
%       uiwait(warndlg(msg,'Error generating file'));
%       return;
%     end
  case 'perframedir',
    [success,msg] = handles.data.GeneratePerFrameFiles(v);
    if ~success,
      uiwait(warndlg(sprintf('Error generating %s files for %s: %s',file,expname,msg)));
    end
end
UpdateStatusTable(handles);


% --- Executes on button press in pushbutton_remove.
function pushbutton_remove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

v = get(handles.listbox_experiment,'Value');
if isempty(v),
  return;
end
handles.data.RemoveExpDirs(v);
set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);
UpdateStatusTable(handles);
set(handles.pushbutton_cancel,'enable','off');

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InitJLabelGui(handles);
handles = guidata(hObject);
SetLabelingMode(handles);
handles.success = true;
guidata(hObject,handles);
uiresume(handles.figure_JLabelEditFiles);

% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

InitJLabelGui(handles);
handles = guidata(hObject);

if handles.data.nexps>0
  res = questdlg('Really load experiment list from classifier file? All changes will be lost','Really load?','Yes','No','Cancel','Yes');
  if ~strcmpi(res,'Yes'),
    return;
  end
end
[filename,pathname] = uigetfile('*.mat','Classifier mat file'); %,handles.classifierfilename);
if ~ischar(filename),
  return;
end
classifierfilename = fullfile(pathname,filename);
if ~exist(classifierfilename,'file'),
  uiwait(warndlg(sprintf('Classifier mat file %s does not exist',classifierfilename),'Error loading file list'));
end

set(handles.pushbutton_load,'enable','off');
set(handles.pushbutton_loadwoexp,'enable','off');
set(handles.pushbutton_cancel,'enable','off');
SetLabelingMode(handles);

[success,msg] = handles.data.SetClassifierFileName(classifierfilename);
if ~success,
  uiwait(warndlg(msg,'Error loading file list'));
  return;
end
set(handles.editClassifier,'String',classifierfilename);
set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);
% update status table
UpdateStatusTable(handles);
uiwait(warndlg('Done loading the classifier'));


% --- Executes on button press in pushbutton_loadwoexp.
function pushbutton_loadwoexp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_loadwoexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

InitJLabelGui(handles);
handles = guidata(hObject);

[filename,pathname] = uigetfile('*.mat','Classifier mat file'); %,handles.classifierfilename);
if ~ischar(filename),
  return;
end
classifierfilename = fullfile(pathname,filename);
if ~exist(classifierfilename,'file'),
  uiwait(warndlg(sprintf('Classifier mat file %s does not exist',classifierfilename),'Error loading file list'));
end

set(handles.pushbutton_load,'enable','off');
set(handles.pushbutton_loadwoexp,'enable','off');
set(handles.pushbutton_cancel,'enable','off');
set(handles.popupmode,'enable','off');
SetLabelingMode(handles);

[success,msg] = handles.data.SetClassifierFileNameWoExp(classifierfilename);  
if ~success,
  uiwait(waitdlg(msg,'Error loading file list'));
  return;
end
set(handles.editClassifier,'String',classifierfilename);
set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);
% update status table
UpdateStatusTable(handles);


% --- Executes when user attempts to close figure_JLabelEditFiles.
function figure_JLabelEditFiles_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_JLabelEditFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% res = questdlg('Really ignore all changes?','','Yes','No','Cancel','Yes');
% if ~strcmpi(res,'Yes'),
%   return;
% end
uiresume(handles.figure_JLabelEditFiles);


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_JLabelEditFiles_CloseRequestFcn(hObject, eventdata, handles);


function SetStatusEditFiles(hObject,s)
handles = guidata(hObject);
set(handles.statusMsg,'String',s);

function ClearStatusEditFiles(hObject)
handles = guidata(hObject);
set(handles.statusMsg,'String','');


% --- Executes during object creation, after setting all properties.
function editClassifier_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editClassifier (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmode.
function popupmode_Callback(hObject, eventdata, handles)
% hObject    handle to popupmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmode


% --- Executes during object creation, after setting all properties.
function popupmode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_behavior.
function listbox_behavior_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_behavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_behavior contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_behavior

handles.projmanager.SetCurrentProject(get(hObject,'Value'));
updateConfigParams(handles);

% --- Executes during object creation, after setting all properties.
function listbox_behavior_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_behavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_addbehavior.
function pushbutton_addbehavior_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addbehavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[~,newName,configFile,defaultConfigFile] = NewProjectSettings(handles.projmanager.GetProjectList);
if isempty(newName),
  return;
end

handles.projmanager.AddProject(newName,configFile,defaultConfigFile);

if ~updateConfigParams(handles)
  handles.projmanager.RemoveProject(handles.projmanager.GetCurrentProject());
end

set(handles.listbox_behavior,'String',handles.projmanager.GetProjectList());
set(handles.listbox_behavior,'Value',handles.projmanager.GetCurrentProject());
updateConfigParams(handles);

% --- Executes on button press in pushbutton_removebehavior.
function pushbutton_removebehavior_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removebehavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.projmanager.RemoveProject(handles.projmanager.GetCurrentProject());

set(handles.listbox_behavior,'String',handles.projmanager.GetProjectList());
set(handles.listbox_behavior,'Value',handles.projmanager.GetCurrentProject());
updateConfigParams(handles);


function success = updateConfigParams(handles)
[data, success] = handles.projmanager.GetConfigAsTable();
if success
  set(handles.config_table,'Data',data);
else
  uiwait(warndlg('Configuration file is out of whack!'));
end


function SetLabelingMode(handles)
contents = cellstr(get(handles.popupmode,'String'));
curStr = contents{get(handles.popupmode,'Value')};
switch curStr,
  case 'Advanced',
    handles.data.SetAdvancedMode(true);
    handles.data.SetGTMode(false);
  case 'Normal'
    handles.data.SetAdvancedMode(false);
    handles.data.SetGTMode(false);
  case 'Ground Truthing',
    handles.data.SetAdvancedMode(false);
    handles.data.SetGTMode(true);
  case 'Ground Truthing Advanced',
    handles.data.SetAdvancedMode(true);
    handles.data.SetGTMode(true);
end
 
handles.data.SetMode();
set(handles.popupmode,'enable','off');


% --- Executes when entered data in editable cell(s) in config_table.
function config_table_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to config_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

if isempty(eventdata.Indices); return; end
data = get(handles.config_table,'Data');
ndx = eventdata.Indices(1);
if eventdata.Indices(2) == 2
  handles.projmanager.EditConfigValue(data{ndx,1},data{ndx,2});
else
  handles.projmanager.EditConfigName(eventdata.PreviousData,eventdata.NewData);
end



% --- Executes on button press in pushbutton_chooseperframe.
function pushbutton_chooseperframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_chooseperframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.projmanager.GetCurrentProject(), return; end;
[allpf selected missing] = handles.projmanager.GetAllPerframeList();
data = {}; data(:,1) = allpf; data(:,2) = num2cell(selected);
if ~isempty(missing)
  list = missing{1};
  for ndx = 2:numel(missing)
    list = sprintf('%s %s ',list,missing{ndx});
  end
  wstr = sprintf('Perframe feature(s) %s are defined in config file%s\n%s',...
        list,...
        ' but is not defined in featureConfig.',...
        'Ignoring them.');
  uiwait(warndlg(wstr));
end

curFigPos = get(handles.figure_JLabelEditFiles,'Position');

fig = figure('Position',[curFigPos(1)+100 curFigPos(2) + 100 300 500],'Visible','off');
fighandle = guidata(fig);
fighandle.editfileshandle = handles;
fighandle.fig = fig;
t = uitable();
fighandle.t = t;
set(t,'Parent',fig,'Units','normalized');
set(t,'Data',data,'Position',[0.001 0.11 0.998 0.89],...
  'ColumnEditable',[false true]);
set(t,'ColumnWidth',{230 40},'ColumnName',{'Perframe Feature','Select'});

donebutton = uicontrol('parent',fig,'Units','normalized','Style','pushbutton','String','Done','Value',1);
set(donebutton,'Position',[0.05 0.025 0.4 0.075],'Callback',@pushbuttonEditTableCallback);

cancelbutton = uicontrol('parent',fig,'Units','normalized','Style','pushbutton','String','Cancel','Value',2);
set(cancelbutton,'Position',[0.55 0.025 0.4 0.075],'Callback',@pushbuttonEditTableCallback);

set(fig,'Visible','on');
jscrollpane = findjobj(t);
rowHeaderViewport=jscrollpane.getComponent(4);
rowHeader=rowHeaderViewport.getComponent(0);
newWidth=0; 
rowHeaderViewport.setPreferredSize(java.awt.Dimension(newWidth,0));
height=rowHeader.getHeight;
rowHeader.setPreferredSize(java.awt.Dimension(newWidth,height));
rowHeader.setSize(newWidth,height); 
guidata(fig,fighandle);

function pushbuttonEditTableCallback(hObject,eventdata,handles)
% Accept or reject changes to the table to select perframe features.
handles = guidata(hObject);
if strcmp(get(hObject,'String'),'Done') 
  editfileshandle = handles.editfileshandle;
  data = get(handles.t,'Data');
  curpf = {};
  for ndx = 1:size(data,1)
    if data{ndx,2}
      curpf{end+1} = data{ndx,1};
    end
  end
  editfileshandle.projmanager.SetPerframeList(curpf);
end

delete(handles.fig);


% --- Executes on button press in pushbutton_addconfigparam.
function pushbutton_addconfigparam_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addconfigparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
in = inputdlg({'Configuration Parameter Name','Configuration Parameter Value'});
if isempty(in) || numel(in) < 2
  return;
end

if ~handles.projmanager.AddConfig(in{1},in{2})
  uiwait(warndlg(sprintf('Configuration Parameter named %s already exists',in{1})));
  return;
end

updateConfigParams(handles);

% --- Executes on button press in pushbutton_removeconfigparam.
function pushbutton_removeconfigparam_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removeconfigparam (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

data = get(handles.config_table,'Data');

jUIScrollPane = findjobj(handles.config_table);
jUITable = jUIScrollPane.getViewport.getView;
allndx = jUITable.getSelectedRows + 1;
if numel(allndx)==1 && allndx <1, return, end

for ndx = allndx(:)'
  handles.projmanager.RemoveConfig(data{ndx,1});
end
updateConfigParams(handles);


% --- Executes on button press in pushbutton_addlist.
function pushbutton_addlist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

oldv = get(handles.listbox_experiment,'Value');

if ~handles.disableBehavior && isempty(handles.projmanager.GetCurrentProject()), 
  uiwait(warndlg('Select a project before adding an experiment'));
  return
end

InitJLabelGui(handles);
handles = guidata(hObject);

% ask user for experiment directory
[listfile,pname] = uigetfile('*.txt','Add experiments list from text file');
listfile = fullfile(pname,listfile);
if ~ischar(listfile),
  return;
end

fid = fopen(listfile,'r');
if fid<0, 
  uiwait(warndlg(sprintf('Cannot open %s for reading',listfile)));
  return;
end

set(handles.pushbutton_cancel,'enable','off');
SetLabelingMode(handles);

expdir = fgetl(fid);

while(ischar(expdir))
  if ismember(expdir,handles.data.expdirs),
    uiwait(warndlg(sprintf('Experiment directory %s already added',expdir),'Already added'));
    expdir = fgetl(fid);
    continue;
  end
  
  
  [success,msg] = handles.data.AddExpDir(expdir);
  if ~success,
    uiwait(warndlg(sprintf('Error adding expdir %s: %s',expdir,msg)));
    expdir = fgetl(fid);
    continue;
  end
  
  expdir = fgetl(fid);
end

fclose(fid);

set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);

% update status table
UpdateStatusTable(handles);
