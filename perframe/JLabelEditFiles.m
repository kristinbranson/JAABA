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

% Last Modified by GUIDE v2.5 12-Nov-2012 13:33:54

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


%--------------------------------------------------------------------------
% --- Executes just before JLabelEditFiles is made visible.
function JLabelEditFiles_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JLabelEditFiles (see VARARGIN)

% have not pushed the done button yet
handles.success = false;

% Certain actions trigger the mode popup menu to be disabled, thereby
% fixing the mode.
handles.modeFixed=false;
handles.groundTruthingMode=false;

%
[ disableBehavior, ...
  figureJLabel] = ...
  myparse(varargin,...
          'disableBehavior',false,...
          'figureJLabel',[]);

% % Just hide the mode stuff entirely        
% set(handles.popupmode,'Visible','off');
% set(handles.labeling_mode_text,'Visible','off');
        
set(handles.popupmode,'String', ...
                      {'Normal', ...
                       'Ground Truthing'});
if handles.groundTruthingMode,
  set(handles.popupmode,'value',2);
else
  set(handles.popupmode,'value',1);
end
                     
% disableBehvaior is true iff the part of the UI for choosing the
% behavior/project should be disabled, such as when JLabelEditFiles is
% reached from JLabel's File > Edit files... menu
handles.disableBehavior = disableBehavior;

if disableBehavior,
  % we go here if reached via "Edit Files..."
  handles.success = true;
  % Simply edit files, not editing the behavior part
  handles.figureJLabel = figureJLabel;
  %handles.data = JLabelHandle.guidata.data;
  
  % get the current labeling mode for JLabel "object", set the popup to
  % match
  handles.data = JLabel('getJLabelData',figureJLabel);
  if JLabel('isGroundTruthingMode',figureJLabel)
    set(handles.popupmode,'Value',find(strcmp(get(handles.popupmode,'String'),'Ground Truthing')));
  else
    set(handles.popupmode,'Value',find(strcmp(get(handles.popupmode,'String'),'Labeling')));
  end
  
  % configfilename=handles.JLabelHandle.guidata.configfilename
  configfilename=JLabel('getConfigFileName',figureJLabel);
  set(handles.text_projectfile,'String',configfilename);
  handles.configfilename=configfilename;
  DisableBehaviorGui(handles);
  %InitExperimentsGui(hObject,handles,JLabelHandle.guidata.data);
  InitExperimentsGui(hObject,handles);
  handles.needJLabelInit = false;
  guidata(hObject,handles);
else
  % We go here if reached via "Open old-school files..."
  configfilename=JLabel('getConfigFileName',figureJLabel);
  if ~isempty(configfilename)
    handles.configfilename = configfilename;
    set(handles.text_projectfile,'String',handles.configfilename);
  else
    previousConfigFileName=JLabel('getPreviousConfigFileName',figureJLabel);
    if isempty(previousConfigFileName)
      handles.configfilename = '';
    else
      handles.configfilename = previousConfigFileName;
      set(handles.text_projectfile,'String',handles.configfilename);
    end
  end
  handles.figureJLabel = figureJLabel;
  handles.needJLabelInit = true;
  guidata(hObject,handles);
end

% Add color for mac's.
buttons = findall(hObject,'Style','pushbutton');
for ndx = 1:numel(buttons)
  SetButtonImage(buttons(ndx));
end

%if ~isempty(handles.JLabelSplashHandle) && ishandle(handles.JLabelSplashHandle),
%  delete(handles.JLabelSplashHandle);
%end

% for debugging, make this window non-modal
%set(hObject,'windowstyle','normal');
set(hObject,'windowstyle','modal');

% Don't need this anymore
%uiwait(handles.figure_JLabelEditFiles);

return


%--------------------------------------------------------------------------
function DisableBehaviorGui(handles)
set(handles.text_projectfile,'enable','off');
set(handles.pushbutton_selectproject,'enable','off');
set(handles.pushbutton_edit_project,'enable','off');
set(handles.pushbutton_newproject,'enable','off');


%--------------------------------------------------------------------------
function InitJLabelGui(handles)

if ~handles.needJLabelInit, return; end
set(handles.figure_JLabelEditFiles,'pointer','watch');
handles.needJLabelInit = false;

% For behavior --
DisableBehaviorGui(handles);
handles.success = true;
% End for behavior --

% Notify JLabel of the configfilename, have it initialize itself
% accordingly
figureJLabel=handles.figureJLabel;
JLabel('setConfigFileName', ...
       figureJLabel, ...
       handles.configfilename);
     
% Update our pointer to the JLabelData, which has now changed
handles.data = JLabel('getJLabelData',figureJLabel);

% Disable that selector for choosing labeling mode
handles=fixLabelingMode(handles);

% Notify JLabel of the labeling mode selected in our figure, which will
% propagate that to the JLabelData
JLabel('setGroundTruthingMode',figureJLabel,handles.groundTruthingMode);

% Modify self as appropriate for the new project configuration
InitExperimentsGui(handles.figure_JLabelEditFiles,handles);

% Set the pointer back to normal
set(handles.figure_JLabelEditFiles,'pointer','arrow');

return


% -------------------------------------------------------------------------
function InitExperimentsGui(hObject,handles)
% Remove behavior related parts.

handles.data.SetStatusFn( @(s) SetStatusEditFiles(handles.figure_JLabelEditFiles,s));
handles.data.SetClearStatusFn( @() ClearStatusEditFiles(handles.figure_JLabelEditFiles));
% initialize listbox
set(handles.editClassifier,'String',handles.data.classifierfilename);
set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',numel(handles.data.expdirs));
set(handles.statusMsg,'String','');

%figureJLabel=handles.figureJLabel;
% if JLabel('isGroundTruthingMode',figureJLabel),
%   if JLabel('isAdvancedMode',figureJLabel),
%     popupNdx = find(strcmp(get(handles.popupmode,'String'),'Ground Truthing Advanced'));
%   else
%     popupNdx = find(strcmp(get(handles.popupmode,'String'),'Ground Truthing'));
%   end
% elseif JLabel('isAdvancedMode',figureJLabel),
%   popupNdx = find(strcmp(get(handles.popupmode,'String'),'Advanced'));
% else
%   popupNdx = find(strcmp(get(handles.popupmode,'String'),'Normal'));
% end
% set(handles.popupmode,'Value',popupNdx);
%if handles.modeFixed,
%  set(handles.popupmode,'Enable','off');
%end

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


% -------------------------------------------------------------------------
function UpdateStatusTable(handles)
% Update the table labeled "Experiment details:", based on the currently
% selected experiment and the the information in handles.data.

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
  if handles.data.IsRequiredFile(file),
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


% -------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = JLabelEditFiles_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1}=hObject;
%varargout{1} = handles.JLabelHandle;
%varargout{2} = handles.success;
%delete(hObject);
return


% -------------------------------------------------------------------------
% --- Executes on selection change in listbox_experiment.
function listbox_experiment_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to listbox_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_experiment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_experiment
UpdateStatusTable(handles);


% -------------------------------------------------------------------------
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


%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.disableBehavior,
  if ~isfield(handles,'configfilename') || isempty(handles.configfilename), 
    uiwait(warndlg('Select a project before adding an experiment'));
    return;
  end
end

InitJLabelGui(handles);
handles = guidata(hObject);

% ask user for experiment directory
% defaultdir = fileparts(handles.data.defaultpath);
% allexpdirs = uipickfiles('FilterSpec',defaultdir,...
%   'Prompt','Add experiment directory');
% if (~iscell(allexpdirs) && isnumeric(allexpdirs)) || (numel(allexpdirs)==1 && ~ischar(allexpdirs{1})),
%   return;
% end

defaultdir = fileparts(handles.data.defaultpath);

allexpdirs = uigetdir2(defaultdir,'Add experiment directory');
if ischar(allexpdirs),
  allexpdirs = {allexpdirs};
end
if isempty(allexpdirs) || ~iscell(allexpdirs),
  return;
end

set(handles.pushbutton_cancel,'enable','off');
handles=fixLabelingMode(handles);

for ndx = 1:numel(allexpdirs)
  expdir = allexpdirs{ndx};
  if ismember(expdir,handles.data.expdirs),
    uiwait(warndlg(sprintf('Experiment directory %s already added',expdir),'Already added'));
    return;
  end
  
  SetStatusEditFiles(hObject,sprintf('Adding experiment directory %s',expdir));
  
  [success,msg] = handles.data.AddExpDir(expdir);
  if ~success,
    if iscell(msg)
      uiwait(warndlg(sprintf('Error adding expdir %s: %s',expdir,msg{:})));
    else
      uiwait(warndlg(sprintf('Error adding expdir %s: %s',expdir,msg)));
    end
    ClearStatusEditFiles(hObject);
    
    return;
  end
  
  set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);

end

SetStatusEditFiles(hObject,'Final update to status table...\n');

% update status table
UpdateStatusTable(handles);

ClearStatusEditFiles(hObject);


%--------------------------------------------------------------------------
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


%--------------------------------------------------------------------------
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


%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
InitJLabelGui(handles);
handles = guidata(hObject);
handles=fixLabelingMode(handles);
% add something here to tell JLabel what mode it should be in?
handles.success = true;
guidata(hObject,handles);
%uiresume(handles.figure_JLabelEditFiles);
if handles.disableBehavior
  % Means JLabelEditFiles was invoked via File > Edit files... in JLabel
  % figure
  JLabel('editFilesDone',handles.figureJLabel);
else
  % Means JLabelEditFiles was invoked via File > Import... in JLabel
  % figure
  JLabel('importDone',handles.figureJLabel);
end
delete(gcbf);
return


%--------------------------------------------------------------------------
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

handles=fixLabelingMode(handles);

res = questdlg('Load the labels and the classifier from the file, or just load the classifier?',...
  'Labels?','Load Labels and Classifier','Load Classifier Only','Cancel','Load Labels and Classifier');
if strcmpi(res,'Cancel'), return, end

classifierlabels = strcmpi(res,'Load Labels and Classifier');

[success,msg] = handles.data.SetClassifierFileName(classifierfilename,'classifierlabels',classifierlabels);
if ~success,
  uiwait(warndlg(msg,'Error loading file list'));
  return;
end
set(handles.pushbutton_load,'enable','off');
set(handles.pushbutton_loadwoexp,'enable','off');
set(handles.pushbutton_cancel,'enable','off');

set(handles.editClassifier,'String',classifierfilename);
set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);
% update status table
UpdateStatusTable(handles);
% Don't need the dialog anymore, b/c using the watch cursor
%uiwait(helpdlg('Done loading the classifier'));

return


%--------------------------------------------------------------------------
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
  return;
end

handles=fixLabelingMode(handles);

[success,msg] = handles.data.SetClassifierFileNameWoExp(classifierfilename);  
if ~success,
  uiwait(warndlg([msg ' Error loading file list']));
  return;
end
set(handles.pushbutton_load,'enable','off');
set(handles.pushbutton_loadwoexp,'enable','off');
set(handles.pushbutton_cancel,'enable','off');
set(handles.popupmode,'enable','off');

set(handles.editClassifier,'String',classifierfilename);
set(handles.listbox_experiment,'String',handles.data.expdirs,'Value',handles.data.nexps);
% update status table
UpdateStatusTable(handles);


%--------------------------------------------------------------------------
% --- Executes when user attempts to close figure_JLabelEditFiles.
function figure_JLabelEditFiles_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_JLabelEditFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% res = questdlg('Really ignore all changes?','','Yes','No','Cancel','Yes');
% if ~strcmpi(res,'Yes'),
%   return;
% end
%uiresume(handles.figure_JLabelEditFiles);
if handles.disableBehavior
  % Means JLabelEditFiles was invoked via File > Edit files... in JLabel
  % figure
  JLabel('editFilesDone',handles.figureJLabel);
else
  % Means JLabelEditFiles was invoked via File > Import... in JLabel
  % figure
  % nothing to do
end
delete(gcbf);
return


%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure_JLabelEditFiles_CloseRequestFcn(hObject, eventdata, handles);
return


%--------------------------------------------------------------------------
function SetStatusEditFiles(hObject,s)

handles = guidata(hObject);
set(handles.statusMsg,'String',s);
set(handles.figure_JLabelEditFiles,'pointer','watch');


%--------------------------------------------------------------------------
function ClearStatusEditFiles(hObject)
handles = guidata(hObject);
set(handles.statusMsg,'String','');
set(handles.figure_JLabelEditFiles,'pointer','arrow');


%--------------------------------------------------------------------------
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


%--------------------------------------------------------------------------
% --- Executes on selection change in popupmode.
function popupmode_Callback(hObject, eventdata, handles)
% hObject    handle to popupmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 1 = Labeling mode (i.e. normal)
% 2 = Ground-truthing mode
i=get(hObject,'Value');
if (i==2),
  handles.groundTruthingMode=true;
  set(handles.pushbutton_load,'Enable','off');
else
  handles.groundTruthingMode=false;
  set(handles.pushbutton_load,'Enable','on');
end
guidata(gcbf,handles);
return


%--------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function handles=fixLabelingMode(handles)
%contents = cellstr(get(handles.popupmode,'String'));
%curStr = contents{get(handles.popupmode,'Value')};

%handles.data.setLabelingMode(curStr);
handles.modeFixed=true;
set(handles.popupmode,'enable','off');
guidata(handles.figure_JLabelEditFiles,handles);  % commit the handles

return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_addlist.
function pushbutton_addlist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.disableBehavior,
  if ~isfield(handles,'configfilename') || isempty(handles.configfilename), 
    uiwait(warndlg('Select a project before adding an experiment'));
    return;
  end
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
handles=fixLabelingMode(handles);

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


%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_selectproject.
function pushbutton_selectproject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_selectproject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.configfilename)
  defpath = handles.configfilename;
else
  defpath = pwd;
end
[fpath,dpath] = uigetfile('*.mat','Select Project File',defpath);
if fpath == 0,
  return; 
end

configfilename = fullfile(dpath,fpath);

set(handles.text_projectfile,'String',configfilename);
handles.configfilename = configfilename;
guidata(hObject,handles);


%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_edit_project.
function pushbutton_edit_project_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_edit_project (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles,'configfilename') && ~isempty(handles.configfilename)
  [~,~,ext] = fileparts(handles.configfilename);
  if strcmp(ext,'.xml')
    configParams = ReadXMLParams(xml_file);
    if ~iscell(configParams.behaviors.names),
      configParams.behaviors.names = {configParams.behaviors.names};
    end
  elseif strcmp(ext,'.mat')
    configParams = load(handles.configfilename);
  else
    uiwait(warndlg('Project file has to be either xml or mat file'));
    return;
  end
end   
ProjectSetup('new',false, ...
             'configParams', configParams, ...
             'figureJLabel',handles.figureJLabel, ...
             'figureJLabelEditFiles',handles.figure_JLabelEditFiles);    

return


%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_newproject.
function pushbutton_newproject_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_newproject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ProjectSetup('new',true, ...
             'figureJLabel',handles.figureJLabel, ...
             'figureJLabelEditFiles',handles.figure_JLabelEditFiles);    
return


%--------------------------------------------------------------------------
function projectSetupDone(figureJLabelEditFiles,configParams,new)
% Tells the JLabelEditFiles "object" that the user just clicked "Done" in
% the Project Setup figure.
% configParams are the just-set-up project params
% new is a boolean, true iff the user pressed "New Project" to setup the
%   project

% if no valid config params, just return
if isempty(configParams), return; end;

% get the guidata
handles=guidata(figureJLabelEditFiles);

% make sure we have a file name in configfilename
if new,
  %previousConfigFileName= ...
  %  JLabel('getPreviousConfigFileName',handles.figureJLabel);
  primaryBehaviorName=strtrim(configParams.behaviors.names{1});
  suggestedFileName=sprintf('%s_project.mat',primaryBehaviorName);
  % a new project, need to get a file name
  [fname,pname] = ...
    uiputfile('*.mat', ...
              'Select a location to store the project file',....
              suggestedFileName);
  if fname == 0; return; end;
  configFileName=fullfile(pname,fname);
  if exist(configFileName,'file')
    [didback,msg] = copyfile(configFileName,[configFileName '~']);
    if ~didback,
      warning('Could not create backup of %s: %s',configFileName,msg);  %#ok
    end
  end
else
  configFileName=handles.configfilename;
end

% try to save to $configfilename
try
  save(configFileName,'-struct','configParams');
catch  %#ok
  uiwait(errordlg('Unable to save project file %s.',fname));
  return;
end

% if save worked, commit the configfilename (only matters if new)
handles.configfilename=configFileName;
%JLabel('setProjectParams',handles.figureJLabel,configParams);

% update the UI
set(handles.text_projectfile,'String',configFileName);

% store the guidata
guidata(figureJLabelEditFiles,handles);

return


%--------------------------------------------------------------------------
% function projectSetupCancelled(figureJLabelEditFiles)
% return

%--------------------------------------------------------------------------
