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

% Last Modified by GUIDE v2.5 30-Jun-2011 18:40:04

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

% parse inputs, set defaults
if isempty(varargin),
  error('Usage: JLabelEditFiles(configfilename,...)');
end
if ischar(varargin{1}),
  configfilename = varargin{1};
  handles.data = JLabelData(configfilename,varargin{2:end});
else
  handles.data = varargin{1};
end

handles.data.SetStatusFn( @(s) SetStatusEditFiles(handles.figure_JLabelEditFiles,s));
handles.data.SetClearStatusFn( @() ClearStatusEditFiles(handles.figure_JLabelEditFiles));
% initialize listbox
set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',numel(handles.data.expdirs));
set(handles.statusMsg,'String','');
% height of a row in pixels
handles.table_row_height_px = 22;

% sizes for generate buttons
handles.generate_button_border_y_px = 1;
handles.generate_button_border_x_px = 5;
handles.generate_button_width_px = 100;
handles.generate_button_height_px = handles.table_row_height_px - handles.generate_button_border_y_px;

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
    handles.pushbutton_generate = uicontrol('Style','pushbutton','Units','Pixels',...
      'Parent',handles.figure_JLabelEditFiles,...
      'String','Generate','Position',pos,...
      'FontUnits','pixels','FontSize',14,...
      'tag',sprintf('pushbutton_generate_%d',i),...
      'Callback',@(hObject,eventdata) JLabelEditFiles('pushbutton_generate_Callback',hObject,eventdata,guidata(hObject),i));
  end
end

buttonNames = {'pushbutton_add','pushbutton_remove','pushbutton_load',...
              'pushbutton_cancel','pushbutton_done'};
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

% UIWAIT makes JLabelEditFiles wait for user response (see UIRESUME)
uiwait(handles.figure_JLabelEditFiles);

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
    timestamp = datestr(timestamp,'yyyymmddTHHMMSS');
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
set(handles.uitable_status,'Data',data,'Visible','on');

% --- Outputs from this function are returned to the command line.
function varargout = JLabelEditFiles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.data;
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

% ask user for experiment directory
expdir = uigetdir(handles.data.defaultpath,'Add experiment directory');
if ~ischar(expdir),
  return;
end

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

set(handles.pushbutton_cancel,'enable','off');

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
handles.success = true;
guidata(hObject,handles);
uiresume(handles.figure_JLabelEditFiles);

% --- Executes on button press in pushbutton_load.
function pushbutton_load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
[success,msg] = handles.data.SetClassifierFileName(classifierfilename);
if ~success,
  uiwait(waitdlg(msg,'Error loading file list'));
  return;
end

set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);
% update status table
UpdateStatusTable(handles);
set(handles.pushbutton_load,'enable','off');
set(handles.pushbutton_cancel,'enable','off');


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
