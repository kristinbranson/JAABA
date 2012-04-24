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

% Last Modified by GUIDE v2.5 16-Apr-2012 17:02:36

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
  JLabelHandle] = ...
  myparse(varargin,...
  'behaviorParamsFile',fullfile('params','BehaviorList.xml'),...
  'disableBehavior',false,...
  'JLabelHandle',[]);

handles.behaviorXMLeditFile = 'xml_edit_tools/behaviorConfigDefaults.xml';
set(handles.popupmode,'String',{'Normal','Advanced','Ground Truthing'});
handles.disableBehavior = disableBehavior;

if disableBehavior,
  handles.success = true;
  % Simply edit files, not editing the behavior part
  handles.JLabelHandle = JLabelHandle;
  handles.data = JLabelHandle.data;
  if handles.data.IsGTMode()
    set(handles.popupmode,'Value',find(strcmp(get(handles.popupmode,'String'),'Ground Truthing')));
  elseif handles.data.IsAdvancedMode()
    set(handles.popupmode,'Value',find(strcmp(get(handles.popupmode,'String'),'Advanced')));
  else
    set(handles.popupmode,'Value',find(strcmp(get(handles.popupmode,'String'),'Normal')));
  end
  DisableBehaviorGui(handles);
  InitExperimentsGui(hObject,handles,JLabelHandle.data);
  handles.needJLabelInit = false;
  boxpos = get(handles.uipanelBehavior,'Position');
  guipos = get(hObject,'Position');
  set(hObject,'Position',[guipos(1:3) boxpos(2)-5]);

  guidata(hObject,handles);
  
else
  
  handles.JLabelHandle = JLabelHandle;
  handles.needJLabelInit = true;
  handles.defaultConfig = GetDefaultConfig();
  handles.behaviorParamsFile = behaviorParamsFile;
  if exist(behaviorParamsFile,'file'),
    behaviorparams.behaviors = ReadXMLParams(behaviorParamsFile);
  else
    behaviorparams.behaviors = struct;
  end
  
  handles.behaviorparams = behaviorparams;
  behaviorList = fieldnames(behaviorparams.behaviors);
  handles.behaviorList = behaviorList;
  handles.needSave = false(1,numel(behaviorList));
  
  handles.config = cell(1,numel(behaviorList));
  handles.featureList = cell(1,numel(behaviorList));  
  
  for ndx = 1:numel(behaviorList)
    behaviorName = behaviorList{ndx};
    curparams = ReadXMLParams(handles.behaviorparams.behaviors.(behaviorName).configFile);
    if isfield(curparams,'featureparamlist'),
      handles.featureList{ndx} = curparams.featureparamlist;
      curparams = rmfield(curparams,'featureparamlist');
    end
    handles.config{ndx} = curparams;
  end
  
  set(handles.listbox_behavior,'String',behaviorList);
  
  if ~isempty(handles.JLabelHandle.configfilename)
    foundConfig = false;
    for fndx = 1:numel(handles.behaviorList)
      if strcmp(handles.JLabelHandle.configfilename,  ...
          handles.behaviorparams.behaviors.(handles.behaviorList{fndx}).configFile),
        handles.curbehavior = fndx;
        set(handles.listbox_behavior,'Value',handles.curbehavior);
        foundConfig = true;
        break;
      end
    end
    if ~foundConfig
      while true
        dlgStr = sprintf('%s\n%s\n%s','A project has not been defined for the configfile',...
          'To add a new project for this configfile enter the name below.',... 
          'To cancel, close this window');
        in = inputdlg(dlgStr);
        if numel(in)<1 || isempty(in{1}),
          break;
        end
        newName = in{1};
        if any(strcmp(newName,handles.behaviorList)),
          uiwait(warndlg('A project already exists with that name'));
        else
          addBehavior(hObject,handles,newName,handles.JLabelHandle.configfilename);
          handles = guidata(hObject);
          break
        end
      end
    end
  else
    if numel(behaviorList)>0,
      handles.curbehavior = 1;
      set(handles.listbox_behavior,'Value',handles.curbehavior);
    else
      handles.curbehavior = [];
    end
  end
  guidata(hObject,handles);
  updateConfigParams(handles);
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


function SaveBehaviorState(handles)
% Writes all the behavior and params file.
writeConfigFile(handles.behaviorParamsFile,...
handles.behaviorparams.behaviors,'behaviors',[]);
for ndx = 1:numel(handles.behaviorList)
  if ~handles.needSave(ndx), continue, end;
  behaviorName = handles.behaviorList{ndx};
  writeConfigFile(handles.behaviorparams.behaviors.(behaviorName).configFile,...
    handles.config{ndx},'params',handles.featureList{ndx});
end


function InitJLabelGui(handles)

if ~handles.needJLabelInit, return; end
handles.needJLabelInit = false;

% For behavior --
DisableBehaviorGui(handles);
SaveBehaviorState(handles);
handles.success = true;
% End for behavior --

% Initializes the JLabel gui once the user selects the behavior.
JLabelHandle = handles.JLabelHandle;
JLabelHandle.configparams = handles.config{handles.curbehavior};
behaviorName = handles.behaviorList{handles.curbehavior};
JLabelHandle.configfilename = handles.behaviorparams.behaviors.(behaviorName).configFile;
JLabelHandle = JLabel('GetGUIPositions',JLabelHandle);
JLabelHandle = JLabel('InitializeState',JLabelHandle);
JLabelHandle = JLabel('InitializePlots',JLabelHandle);
handles.data = JLabelHandle.data;
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
  popupNdx = find(strcmp(get(handles.popupmode,'String'),'Ground Truthing'));
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

%{
% % status table fields: these should be defined in order so that earlier
% % files do not depend on later files
% handles.status_files = struct('name',{},'file',{},'required',{},'cangenerate',{},'perfly',{},'isoutput',{});
% handles.status_files(end+1).name = 'movie';
% handles.status_files(end).file = params.moviefilename;
% handles.status_files(end).required = true;
% handles.status_files(end).cangenerate = false;
% handles.status_files(end).perfly = false;
% handles.status_files(end).isoutput = false;
% handles.status_files(end+1).name = 'trx';
% handles.status_files(end).file = params.trxfilename;
% handles.status_files(end).required = true;
% handles.status_files(end).cangenerate = false;
% handles.status_files(end).perfly = false;
% handles.status_files(end).isoutput = false;
% handles.status_files(end+1).name = 'label';
% handles.status_files(end).file = params.labelfilename;
% handles.status_files(end).required = false;
% handles.status_files(end).cangenerate = false;
% handles.status_files(end).perfly = false;
% handles.status_files(end).isoutput = true;
% handles.status_files(end+1).name = 'perframedir';
% handles.status_files(end).file = params.perframedir;
% handles.status_files(end).required = true;
% handles.status_files(end).cangenerate = true;
% handles.status_files(end).perfly = false;
% handles.status_files(end).isoutput = false;
% handles.status_files(end+1).name = 'window';
% handles.status_files(end).file = params.windowfilename;
% handles.status_files(end).required = true;
% handles.status_files(end).cangenerate = true;
% handles.status_files(end).perfly = true;
% handles.status_files(end).isoutput = true;
%}

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

buttonNames = {'pushbutton_add','pushbutton_remove','pushbutton_load',...
              'pushbutton_loadwoexp','pushbutton_cancel','pushbutton_done'};
for buttonNum = 1:length(buttonNames)
  SetButtonImage(handles.(buttonNames{buttonNum}));
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
if handles.data.nexps == 0 || isempty(v) || v <= 0,
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

if ~handles.disableBehavior && isempty(handles.curbehavior), 
  uiwait(warndlg('Select a project before adding an experiment'));
  return
end

InitJLabelGui(handles);
handles = guidata(hObject);

% ask user for experiment directory
expdir = uigetdir(handles.data.defaultpath,'Add experiment directory');
if ~ischar(expdir),
  return;
end

if ismember(expdir,handles.data.expdirs),
  uiwait(warndlg(sprintf('Experiment directory %s already added',expdir),'Already added'));
  return;
end

set(handles.pushbutton_cancel,'enable','off');
SetLabelingMode(handles);

[success,msg] = handles.data.AddExpDir(expdir);
if ~success,
  uiwait(warndlg(sprintf('Error adding expdir %s: %s',expdir,msg)));
  return;
end

set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);

% update status table
UpdateStatusTable(handles);

% check for existence of necessary files in this directory
if ~handles.data.filesfixable,
  uiwait(warndlg(sprintf('Experiment %s is missing required files that cannot be generated within this interface. Removing...',expdir),'Bad experiment'));
  % undo
  handles.data.RemoveExpDirs(handles.data.nexps);
  set(handles.listbox_experiment,'Value',oldv);
  UpdateStatusTable(handles);  
end

if handles.data.filesfixable && ~handles.data.allfilesexist,
  res = questdlg(sprintf('Experiment %s is missing required files. Generate now?',expdir),'Generate missing files?','Yes','Cancel','Yes');
  if strcmpi(res,'Yes'),
    [success,msg] = handles.data.GenerateMissingFiles(handles.data.nexps);
    if ~success,
      uiwait(warndlg(sprintf('Error generating missing required files for experiment %s: %s. Removing...',expdir,msg),'Error generating files'));
      % undo
      handles.data.RemoveExpDirs(handles.data.nexps);
      set(handles.listbox_experiment,'Value',oldv);
    end
    
    [success,msg] = handles.data.PreLoadLabeledData();
    if ~success,
      uiwait(warndlg(sprintf('Error computing window data for experiment %s: %s. Removing...',expdir,msg),'Error Computing Window Data'));
      handles.data.RemoveExpDirs(handles.data.nexps);
      set(handles.listbox_experiment,'Value',oldv);
    end
      
  else
    % undo
    handles.data.RemoveExpDirs(handles.data.nexps);
    set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',oldv);
  end
  UpdateStatusTable(handles);
end



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

res = questdlg('Really load experiment list from classifier file? All changes will be lost','Really load?','Yes','No','Cancel','Yes');
if ~strcmpi(res,'Yes'),
  return;
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
  uiwait(waitdlg(msg,'Error loading file list'));
  return;
end
set(handles.editClassifier,'String',classifierfilename);
set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);
% update status table
UpdateStatusTable(handles);

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


% --- Executes on button press in pushbutton_editbehavior.
function list = getPathAndValue(list,pathTillNow,name,param)
  if isstruct(param)
    fnames = fieldnames(param);
    for ndx = 1:numel(fnames)
      list = getPathAndValue(list,[pathTillNow  name '.'],fnames{ndx},param.(fnames{ndx}));
    end
  else
    list{end+1,1} = [pathTillNow name];
    if isnumeric(param)
      q = num2str(param(1));
      for jj = 2:numel(param)
        q = [q ',' num2str(param(jj))];
      end
      list{end,2} = q;
    else
      list{end,2} = param;
    end
  end


% --- Executes on selection change in listbox_behavior.
function listbox_behavior_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_behavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_behavior contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_behavior

handles.curbehavior = get(hObject,'Value');
guidata(hObject,handles);
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

[~,newName,configFile,defaultConfigFile] = NewProjectSettings(handles.behaviorList);
if isempty(newName),
  return;
end

addBehavior(hObject,handles,newName,configFile,defaultConfigFile);

function addBehavior(hObject,handles,newName,configFile,defaultConfigFile)

if nargin <4,
  defaultConfigFile = '';
end

handles.behaviorparams.behaviors.(newName).configFile = configFile;
handles.behaviorList{end+1} = newName;
handles.needSave(end+1) = true;
handles.curbehavior = numel(handles.behaviorList);
set(handles.listbox_behavior,'String',handles.behaviorList,'Value',handles.curbehavior);

fileToRead = '';

if exist(configFile,'file')
  fileToRead = configFile;
  handles.needSave(handles.curbehavior) = false;
elseif ~isempty(defaultConfigFile) && exist(defaultConfigFile,'file')
  fileToRead = defaultConfigFile;
end

if ~isempty(fileToRead)
  curparams = ReadXMLParams(fileToRead);
  if isfield(curparams,'featureparamlist'),
    handles.featureList{handles.curbehavior} = curparams.featureparamlist;
    curparams = rmfield(curparams,'featureparamlist');
  else
    handles.featureList{handles.curbehavior} = [];    
  end
  handles.config{handles.curbehavior} = curparams;
else
    handles.config{handles.curbehavior} = GetDefaultConfig(newName);
    handles.featureList{handles.curbehavior} = [];
    handles.needSave(handles.curbehavior) = true;
end

set(handles.listbox_behavior,'Value',handles.curbehavior);
guidata(hObject,handles);
if ~updateConfigParams(handles)
  pushbutton_removebehavior_Callback(...
    handles.pushbutton_removebehavior,[],handles);
end


% --- Executes on button press in pushbutton_removebehavior.
function pushbutton_removebehavior_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removebehavior (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curname = handles.behaviorList{handles.curbehavior};
handles.behaviorList(handles.curbehavior) = [];
handles.needSave(handles.curbehavior) = [];
handles.behaviorparams.behaviors = rmfield(handles.behaviorparams.behaviors,curname);
handles.config(handles.curbehavior) = [];
handles.featureList(handles.curbehavior) = [];
if numel(handles.featureList)<handles.curbehavior,
  handles.curbehavior = handles.curbehavior -1;
end
if handles.curbehavior==0 && numel(handles.behaviorList)>0,
  handles.curbehavior = 1;
end
if ~handles.curbehavior, handles.curbehavior =[]; end
set(handles.listbox_behavior,'String',handles.behaviorList,'Value',handles.curbehavior);
guidata(hObject,handles);
updateConfigParams(handles);

guidata(hObject,handles);

function success = updateConfigParams(handles)
success = true;
if isempty(handles.curbehavior); 
  set(handles.config_table,'Data',{});
  return; 
end
configparams = handles.config{handles.curbehavior};
curDat = addToList(configparams,{},'');
if any(cellfun(@iscell,curDat(:,2))),
  uiwait(warndlg('Configuration file is out of whack!'));
  success = false;
  return;
end
set(handles.config_table,'Data',curDat);

function list = addToList(curStruct,list,pathTillNow)
if isempty(fieldnames(curStruct)), return; end
fnames = fieldnames(curStruct);
for ndx = 1:numel(fnames)
  if isstruct(curStruct.(fnames{ndx})), 
    list = addToList(curStruct.(fnames{ndx}),list,[pathTillNow fnames{ndx} '.']);
  else
    list{end+1,1} = [pathTillNow fnames{ndx}];
    param = curStruct.(fnames{ndx});
    if isnumeric(param)
      q = num2str(param(1));
      for jj = 2:numel(param)
        q = [q ',' num2str(param(jj))];
      end
      list{end,2} = q;
    else
      list{end,2} = param;
    end
  end
end

function writeConfigFile(fname,topNode,topNodeName,featureList)
docNode = com.mathworks.xml.XMLUtils.createDocument(topNodeName);
toc = docNode.getDocumentElement;
att = fieldnames(topNode);
for ndx = 1:numel(att)
  toc.appendChild(createXMLNode(docNode,att{ndx},topNode.(att{ndx})));
end
if ~isempty(featureList),
  toc.appendChild(createXMLNode(docNode,'featureparamlist',featureList));
end
xmlwrite(fname,docNode);


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
  eval(sprintf('handles.config{handles.curbehavior}.%s = data{ndx,2};',...
      data{ndx,1}));
else
  [fpath,lastfield] = splitext(eventdata.PreviousData);
  if isempty(lastfield)
    handles.config{curbehavior} = rmfield(handles.config{curbehavior},fpath);
  else
    eval(sprintf('handles.config{handles.curbehavior}.%s = rmfield(handles.config{handles.curbehavior}.%s,lastfield(2:end));',...
      fpath,fpath));  
  end
  eval(sprintf('handles.config{handles.curbehavior}.%s = data{ndx,2};',...
      eventdata.NewData));
  
end

handles.needSave(handles.curbehavior) = true;
guidata(hObject,handles);


% --- Executes on button press in pushbutton_chooseperframe.
function pushbutton_chooseperframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_chooseperframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.curbehavior), return; end;
curbehavior = handles.curbehavior;
featureConfig = ReadXMLParams(handles.config{curbehavior}.file.featureconfigfile);
allperframe = fieldnames(featureConfig.perframe);
data = {};
data(:,1) = allperframe;
if isempty(handles.featureList{curbehavior})
  data(:,2) = {true};
else
  data(:,2) = {false};
  curperframe = fieldnames(handles.featureList{curbehavior});
  for ndx = 1:numel(curperframe)
    allndx = find(strcmp(curperframe{ndx},allperframe));
    if isempty(allndx)
      wstr = sprintf('Perframe feature %s is defined in config file%s\n%s',...
        curperframe{ndx},...
        'but is not defined in featureConfig.',...
        'Ignoring it.');
        uiwait(warndlg(wstr));
      continue;
    end
    data{allndx,2} = true;
  end
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
  curbehavior = editfileshandle.curbehavior;
  editfileshandle.featureList{curbehavior} = struct;
  for ndx = 1:size(data,1)
    if data{ndx,2}
      editfileshandle.featureList{curbehavior}.(data{ndx,1}) = struct;
    end
  end
  editfileshandle.needSave(curbehavior) = true;
  guidata(editfileshandle.figure_JLabelEditFiles,editfileshandle);
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

curbehavior = handles.curbehavior;
if isfield(handles.config{curbehavior},in{1}),
  uiwait(warndlg(sprintf('Configuration Parameter named %s already exists',in{1})));
  return;
end

eval(sprintf('handles.config{curbehavior}.%s=in{2};',in{1}));
handles.needSave(curbehavior) = true;
guidata(hObject,handles);
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
curbehavior = handles.curbehavior;

for ndx = allndx(:)'
  [fpath,lastfield] = splitext(data{ndx,1});
  if isempty(lastfield)
    handles.config{curbehavior} = rmfield(handles.config{curbehavior},fpath);
  else
    handles.config{curbehavior}.(fpath) = rmfield(handles.config{curbehavior}.(fpath),lastfield(2:end));  
  end
end
handles.needSave(curbehavior) = true;
guidata(hObject,handles);
updateConfigParams(handles);
