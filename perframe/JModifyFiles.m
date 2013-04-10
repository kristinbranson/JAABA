function varargout = JModifyFiles(varargin)
% JMODIFYFILES MATLAB code for JModifyFiles.fig
%      JMODIFYFILES, by itself, creates a new JMODIFYFILES or raises the existing
%      singleton*.
%
%      H = JMODIFYFILES returns the handle to a new JMODIFYFILES or the handle to
%      the existing singleton*.
%
%      JMODIFYFILES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in JMODIFYFILES.M with the given input arguments.
%
%      JMODIFYFILES('Property','Value',...) creates a new JMODIFYFILES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before JModifyFiles_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to JModifyFiles_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help JModifyFiles

% Last Modified by GUIDE v2.5 21-Jan-2013 00:37:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @JModifyFiles_OpeningFcn, ...
                   'gui_OutputFcn',  @JModifyFiles_OutputFcn, ...
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
% --- Executes just before JModifyFiles is made visible.
function JModifyFiles_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to JModifyFiles (see VARARGIN)

% Need a better name for hObject
figureJModifyFiles=hObject;

% Parse the arguments
figureJLabel = ...
  myparse(varargin,...
          'figureJLabel',[]);
                             
% Store the parent figure
setGuidataField(figureJModifyFiles,'figureJLabel',figureJLabel);

% Store whether the list of experiments has actually changed
setGuidataField(figureJModifyFiles,'listChanged',false);

% Store a direct ref to the data
data=JLabel('getJLabelData',figureJLabel);  % a ref
setGuidataField(figureJModifyFiles,'data',data);

% Add color for Macs.
buttons = findall(figureJModifyFiles,'Style','pushbutton');
for i = 1:numel(buttons)
  %SetButtonImage(buttons(i));
  adjustNonLabelButtonColor(buttons(i));
end

% Make the window modal (But... Does it really to be modal?)
set(figureJModifyFiles,'windowstyle','modal');
% for debugging, make this window non-modal
%set(figureJModifyFiles,'windowstyle','normal');

% Route the status functions
data.SetStatusFn( @(s) SetStatusEditFiles(figureJModifyFiles,s));
data.SetClearStatusFn( @() ClearStatusEditFiles(figureJModifyFiles));

% Initialize some things
listbox_experiment=getGuidataField(figureJModifyFiles,'listbox_experiment');
statusMessageTextbox=getGuidataField(figureJModifyFiles,'statusMessageTextbox');
set(listbox_experiment,'String',data.expdirs,'Value',numel(data.expdirs));
set(statusMessageTextbox,'String','');

% % Specify height of a row in pixels
% table_row_height_px = 22;
% 
% % Specify sizes for generate buttons
% generate_button_border_y_px = 1;
% generate_button_border_x_px = 5;
% generate_button_width_px = 100;
% generate_button_height_px = table_row_height_px - generate_button_border_y_px;
% 
% % Create buttons for generating the files
% uitable_status=getGuidataField(figureJModifyFiles,'uitable_status');
% pos_table = get(uitable_status,'Position');
% top_table = pos_table(2)+pos_table(4);
% right_table = pos_table(1)+pos_table(3);
% for i = 1:numel(data.filetypes),
%   if JLabelData.CanGenerateFile(data.filetypes{i}),
%     pos = [right_table + generate_button_border_x_px,...
%            top_table - (i-1)*table_row_height_px - ...
%            generate_button_height_px - generate_button_border_y_px/2,...
%            generate_button_width_px,...
%            generate_button_height_px];
%     buttonColor = get(figureJModifyFiles,'Color');
%     pushbutton_generate = ...
%       uicontrol('Style','pushbutton', ...
%                 'Units','Pixels',...
%                 'Parent',figureJModifyFiles,...
%                 'String','Generate', ...
%                 'Position',pos, ...
%                 'FontUnits','pixels', ...
%                 'FontSize',14,...
%                 'tag',sprintf('pushbutton_generate_%d',i),...
%                 'Callback',@(hObject,eventdata) JModifyFiles('pushbutton_generate_Callback',hObject,eventdata,guidata(hObject),i),...
%                 'Backgroundcolor',buttonColor);
%     setGuidataField(figureJModifyFiles,'pushbutton_generate',pushbutton_generate);              
%   end
% end

% initialize status table
UpdateStatusTable(figureJModifyFiles);

return


% -------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = JModifyFiles_OutputFcn(hObject, eventdata, handles) 
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
UpdateStatusTable(gcbf);


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function listbox_experiment_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to listbox_experiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), ...
                   get(0,'defaultUicontrolBackgroundColor'))
  set(hObject,'BackgroundColor','white');
  SetButtonImage(hObject);
end


%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_add.
function pushbutton_add_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figureJModifyFiles=gcbf;
data=getGuidataField(figureJModifyFiles,'data');  % ref
defaultdir = fileparts(data.defaultpath);

allexpdirs = uigetdir2(defaultdir,'Add experiment directory');
if ischar(allexpdirs),
  allexpdirs = {allexpdirs};
end
if isempty(allexpdirs) || ~iscell(allexpdirs),
  return;
end

for ndx = 1:numel(allexpdirs)
  expdir = allexpdirs{ndx};
  if ismember(expdir,data.expdirs),
    uiwait(warndlg(sprintf('Experiment directory %s already added',expdir),'Already added'));
    return;
  end
  SetStatusEditFiles(gcbf,sprintf('Adding experiment directory %s',expdir));
  [success,msg] = data.AddExpDir(expdir);
  if ~success,
    if iscell(msg)
      uiwait(warndlg(sprintf('Error adding expdir %s: %s',expdir,msg{:})));
    else
      uiwait(warndlg(sprintf('Error adding expdir %s: %s',expdir,msg)));
    end
    ClearStatusEditFiles(figureJModifyFiles);
    return;
  end
  set(handles.listbox_experiment,'String',data.expdirs,'Value',data.nexps);
  setGuidataField(figureJModifyFiles,'listChanged',true);
end

SetStatusEditFiles(figureJModifyFiles,'Final update to status table...\n');

% update status table
UpdateStatusTable(figureJModifyFiles);

ClearStatusEditFiles(figureJModifyFiles);

return


% %--------------------------------------------------------------------------
% function pushbutton_generate_Callback(hObject, eventdata, handles, row)
% listbox_experiment=getGuidataField(figureJModifyFiles,'listbox_experiment');
% data=getGuidataField(gcbf,'data');  % ref
% file = data.filetypes{row};
% if ~JLabelData.CanGenerateFile(file),
%   return;
% end
% v = get(listbox_experiment,'Value');
% if isempty(v),
%   return;
% end
% if numel(v) > 1,
%   v = v(end);
% end
% expname = data.expnames{v};
% switch file,
%   case 'perframedir',
%     [success,msg] = data.GeneratePerFrameFiles(v);
%     if ~success,
%       uiwait(warndlg(sprintf('Error generating %s files for %s: %s',file,expname,msg)));
%     end
% end
% UpdateStatusTable(hObject);
% return


%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_remove.
function pushbutton_remove_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
listbox_experiment=getGuidataField(gcbf,'listbox_experiment');
data=getGuidataField(gcbf,'data');  % ref
v = get(listbox_experiment,'Value');
if isempty(v),
  return;
end
data.RemoveExpDirs(v);
set(listbox_experiment,'String',data.expdirs,'Value',data.nexps);
setGuidataField(gcbf,'listChanged',true);
UpdateStatusTable(hObject);
return


%--------------------------------------------------------------------------
% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
figureJLabel=getGuidataField(gcbf,'figureJLabel');
listChanged=getGuidataField(gcbf,'listChanged');
JLabel('modifyFilesDone',figureJLabel,listChanged);
delete(gcbf);
return


%--------------------------------------------------------------------------
% --- Executes when user attempts to close figureJModifyFiles.
function figure_JModifyFiles_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figureJModifyFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Do nothing.  The only way out is the "Done" button, to make clear to the
% user that you can't really cancel the things you done in JModifyFiles.
return


%--------------------------------------------------------------------------
function SetStatusEditFiles(figureJModifyFiles,s)
statusMessageTextbox=getGuidataField(figureJModifyFiles, ...
                                     'statusMessageTextbox');
set(statusMessageTextbox,'String',s);
set(figureJModifyFiles,'pointer','watch');
return


%--------------------------------------------------------------------------
function ClearStatusEditFiles(figureJModifyFiles)
statusMessageTextbox=getGuidataField(figureJModifyFiles, ...
                                     'statusMessageTextbox');
set(statusMessageTextbox,'String','');
set(figureJModifyFiles,'pointer','arrow');
return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_addlist.
function pushbutton_addlist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% get the figure handle
figureJModifyFiles=gcbf;

% ask user for experiment directory
[listfile,pname] = uigetfile({'*.txt','Text Files (*.txt)'}, ...
                             'Add experiments list from text file');
if listfile==0 ,
  % means user hit Cancel button
  return
end
listfile = fullfile(pname,listfile);
if ~ischar(listfile),
  return;
end
fid = fopen(listfile,'r');
if fid<0, 
  uiwait(warndlg(sprintf('Cannot open %s for reading',listfile)));
  return;
end

expdir = fgetl(fid);

data=getGuidataField(figureJModifyFiles,'data');  % a ref
while(ischar(expdir))
  if ismember(expdir,data.expdirs),
    uiwait(warndlg(sprintf('Experiment directory %s already added',expdir),'Already added'));
    expdir = fgetl(fid);
    continue;
  end
    
  [success,msg] = data.AddExpDir(expdir);
  if ~success,
    uiwait(warndlg(sprintf('Error adding expdir %s: %s',expdir,msg)));
    expdir = fgetl(fid);
    continue;
  end
  setGuidataField(gcbf,'listChanged',true);
  
  expdir = fgetl(fid);
end

fclose(fid);

listbox_experiment=getGuidataField(figureJModifyFiles,'listbox_experiment');
set(listbox_experiment,'String',data.expdirs,'Value',data.nexps);

% update status table
UpdateStatusTable(figureJModifyFiles);

% Make sure the status is cleared, and that the cursor is normal
ClearStatusEditFiles(figureJModifyFiles);

return


% -------------------------------------------------------------------------
function UpdateStatusTable(figureJModifyFiles)
% Update the table labeled "Experiment details:", based on the currently
% selected experiment and the the information in handles.data.

data=getGuidataField(figureJModifyFiles,'data');  % a ref
listbox_experiment=getGuidataField(figureJModifyFiles,'listbox_experiment');
uitable_status=getGuidataField(figureJModifyFiles,'uitable_status');

v = get(listbox_experiment,'Value');
if isempty(v) || v <= 0 || data.nexps == 0 ,
  set(uitable_status,'Data',{},'Enable','off');
  return;
end
if numel(v) > 1,
  v = v(end);
end

nfiles = numel(data.filetypes); 
tableData = cell([nfiles,2]);
tableData(:,1) = data.filetypes;
for i = 1:nfiles,
  file = data.filetypes{i};
  [file_exists,timestamp] = data.FileExists(file,v);
  if file_exists,
    timestamp = datestr(timestamp);%,'yyyymmddTHHMMSS');
  end
  if data.IsRequiredFile(file),
    if file_exists,
      tableData{i,2} = sprintf('<html><font color="green">%s</font></html>',timestamp);
    else
      tableData{i,2} = '<html><font color="red">Missing</font></html>';
    end
  else
    if file_exists,
      tableData{i,2} = timestamp;
    else
      tableData{i,2} = 'Absent';
    end
  end
end
tableSize = get(uitable_status,'Position');
set(uitable_status,'Data',tableData, ...
                   'Enable','on', ...
                   'ColumnWidth',{150 tableSize(3)-155});
return
