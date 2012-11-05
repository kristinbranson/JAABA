function varargout = ProjectSetup(varargin)
% PROJECTSETUP MATLAB code for ProjectSetup.fig
%      PROJECTSETUP, by itself, creates a new PROJECTSETUP or raises the existing
%      singleton*.
%
%      H = PROJECTSETUP returns the handle to a new PROJECTSETUP or the handle to
%      the existing singleton*.
%
%      PROJECTSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECTSETUP.M with the given input arguments.
%
%      PROJECTSETUP('Property','Value',...) creates a new PROJECTSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProjectSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProjectSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProjectSetup

% Last Modified by GUIDE v2.5 02-Nov-2012 13:22:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProjectSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @ProjectSetup_OutputFcn, ...
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


% --- Executes just before ProjectSetup is made visible.
function ProjectSetup_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProjectSetup (see VARARGIN)

% Choose default command line output for ProjectSetup
set(hObject,'Visible','off');

handles.output = hObject;
handles.outfile = 0;

buttons = findall(hObject,'Style','pushbutton');
for ndx = 1:numel(buttons)
  SetButtonImage(buttons(ndx));
end

buttons = findall(hObject,'Style','togglebutton');
for ndx = 1:numel(buttons)
  SetButtonImage(buttons(ndx));
end

if ismac, % On mac change the foreground color to black.
  allpopups = findall(hObject,'Style','popup');
  set(allpopups,'ForegroundColor',[0 0 0]);
end

handles.mode = 'advanced';
set(handles.togglebutton_advanced,'
curPos = get(handles.figure1,'Position');
tablePos = get(handles.config_table,'Position');
reducedWidth = curPos(3);
reducedHeight = curPos(4) - tablePos(2) - tablePos(4) - 15;
handles.advancedSize = curPos(3:4);
handles.basicSize = [reducedWidth reducedHeight];
updatePosition(handles,'basic');

handles = initFeatureconfig(handles);
handles = initParams(handles);

set(handles.listbox_inputscores,'String',{});

% Update handles structure
guidata(hObject, handles);

updateText(hObject, handles);

set(hObject,'Visible','on');

% UIWAIT makes ProjectSetup wait for user response (see UIRESUME)
uiwait(handles.figure1);


function updatePosition(handles,mode)
if strcmp(handles.mode,mode), return; end
if strcmp(handles.mode,'advanced')
  curPos = get(handles.figure1,'Position');
  newLeft = curPos(1);
  newBottom = curPos(2) + curPos(4) - handles.basicSize(2);
  set(handles.figure1,'Position',[newLeft newBottom handles.basicSize]);
  set(handles.togglebutton_advanced,'Value',0);
else
  curPos = get(handles.figure1,'Position');
  newLeft = curPos(1);
  newBottom = curPos(2) + curPos(4) - handles.advancedSize(2);
  set(handles.figure1,'Position',[newLeft newBottom handles.advancedSize]);
  set(handles.togglebutton_advanced,'Value',1);  
end

function handles = initFeatureconfig(handles)
list = ReadXMLParams('params/featureConfigList.xml');
animal_types = fieldnames(list);

set(handles.featureconfigpopup,'String',animal_types,'Value',1);
handles.params.file.featureconfigfile = handles.list{1};

handles.list = list;
handles.animal_types = animal_types;


function handles = initParams(handles)
handles.params = struct;
handles.params.behaviors = struct;
handles.params.file = struct;
handles.params.file.perframedir = 'perframe';
handles.params.file.clipsdir = 'clips';
vid = fopen('version.txt','r');
vv = textscan(vid,'%s');
fclose(fid);
handles.params.ver = vv{1};
handles.params.scoresinput = struct('classifierfile',{},'ts',{},'scorefilename',{});
handles.params.windowfeatures = struct;

function updateText(handles)
% Copies the parameters in the params structure to the gui.

if isfield(handles.params.behaviors,'names');
  set(handles.editName,'String',handles.params.behaviors.names);
else
  set(handles.editName,'String','');
end

fnames = {'labelfilename','gt_labelfilename','scorefilename',...
  'perframedir','moviefilename','trxfilename','clipsdir'};
boxnames = {'editlabelfilename','editgtlabelfilename','editscorefilename',...
  'editperframedir','editmoviefilename','edittrxfilename','editclipsdir'};

for ndx = 1:numel(fnames)
  curf = fnames{ndx};
  curbox = boxnames{ndx};
  if isfield(handles.params.file,curf);
    set(handles.(curbox),'String',handles.params.file.(curf));
  else
    set(handles.(curbox),'String','');
  end
end

function handles = updateParams(handles)

fnames = {'labelfilename','gt_labelfilename','scorefilename',...
  'perframedir','moviefilename','trxfilename','clipsdir'};
boxnames = {'editlabelfilename','editgtlabelfilename','editscorefilename',...
  'editperframedir','editmoviefilename','edittrxfilename','editclipsdir'};

for ndx = 1:numel(fnames)
  curf = fnames{ndx};
  curbox = boxnames{ndx};
  str = get(handles.(curbox),'String');
  if ~isempty(str)
    handles.params.file.(curf) = str;
  end
end


% --- Outputs from this function are returned to the command line.
function varargout = ProjectSetup_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.outfile;
delete(handles.figure1);

function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editName as text
%        str2double(get(hObject,'String')) returns contents of editName as a double
name = get(hObject,'String');
handles.params.behaviors.names = name;
handles.params.file.labelfilename = sprintf('label_%s',name);
handles.params.file.gt_labelfilename = sprintf('gt_label_%s',name);
handles.params.file.scorefilename = sprintf('scores_%s',name);
updateText(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in featureconfigpopup.
function featureconfigpopup_Callback(hObject, eventdata, handles)
% hObject    handle to featureconfigpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns featureconfigpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from featureconfigpopup

listndx = get(hObject,'Value');
handles.params.file.featureconfigfile = handles.list{value};
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function featureconfigpopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to featureconfigpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editlabelfilename_Callback(hObject, eventdata, handles)
% hObject    handle to editlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editlabelfilename as text
%        str2double(get(hObject,'String')) returns contents of editlabelfilename as a double
handles = updateParams(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function editlabelfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editgtlabelfilename_Callback(hObject, eventdata, handles)
% hObject    handle to editgtlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editgtlabelfilename as text
%        str2double(get(hObject,'String')) returns contents of editgtlabelfilename as a double
handles = updateParams(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editgtlabelfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editgtlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editscorefilename_Callback(hObject, eventdata, handles)
% hObject    handle to editscorefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editscorefilename as text
%        str2double(get(hObject,'String')) returns contents of editscorefilename as a double
handles = updateParams(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editscorefilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editscorefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function editmoviefilename_Callback(hObject, eventdata, handles)
% hObject    handle to editmoviefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmoviefilename as text
%        str2double(get(hObject,'String')) returns contents of editmoviefilename as a double
handles = updateParams(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editmoviefilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editmoviefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edittrxfilename_Callback(hObject, eventdata, handles)
% hObject    handle to edittrxfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittrxfilename as text
%        str2double(get(hObject,'String')) returns contents of edittrxfilename as a double
handles = updateParams(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edittrxfilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edittrxfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editperframedir_Callback(hObject, eventdata, handles)
% hObject    handle to editperframedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editperframedir as text
%        str2double(get(hObject,'String')) returns contents of editperframedir as a double
handles = updateParams(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editperframedir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editperframedir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editclipsdir_Callback(hObject, eventdata, handles)
% hObject    handle to editclipsdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editclipsdir as text
%        str2double(get(hObject,'String')) returns contents of editclipsdir as a double
handles = updateParams(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editclipsdir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editclipsdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on butuiresume(handles.figure1);ton press in pushbutton_setfeatures.
function pushbutton_setfeatures_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setfeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_copyfeatures.
function pushbutton_copyfeatures_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copyfeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox_inputscores.
function listbox_inputscores_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_inputscores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_inputscores contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_inputscores


% --- Executes during object creation, after setting all properties.
function listbox_inputscores_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_inputscores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_addlist.
function pushbutton_addlist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fnames,pname] = uigetfile('*.mat','Classifier whose scores should be used as input');
if fnames == 0; return; end;

cfile = fullfile(pname,fnames);
classifier = load(cfile);
curs = struct;
curs.classifierfile = cfile;
if isfield(classifier,'classifierTS');
  curs.ts = classifier.classifierTS;
else
  uiwait(warndlg('The selected file does not have all the required fields, Check to make sure it is a classifier file'));
  return;
end

if isfield(classifier,'scorefilename')
  curs.scorefilename = scorefilename;
else
  configfile = classifier.configfilename;
  if exist(configfile,'file')
    configparams = ReadXMLParams(configfile);
    if isfield(configparams,'scorefilename');
      curs.scorefilename = configparams.file.scorefilename;
    else
      curs.scorefilename = sprintf('scores_%s.mat',configparams.behaviors.names);
    end
  else
    sname = inputdlg('Cannot find the config file pointed by the classifier. Input the file name used to stored the scores computed by the classifier',...
      'Score file name');
    if isempty(sname), return; end
    curs.scorefilename = sname;
  end
end

handles.params.scoresinput(end+1) = curs;
guidata(hObject,handles);


% --- Executes on button press in pushbutton_removelist.
function pushbutton_removelist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curndx = get(handles.listbox_inputscores,'Value');
if isempty(curndx), return; end
handles.params.scoresinput(curndx) = [];
clist = {handles.params.scoresinput(:).classifierfile};
if isempty(clist),
  set(handles.listbox_inputscores,'String',{});
else
  set(handles.listbox_inputscores,'String',clist,'Value',1);
end
guidata(hObject,handles);

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.outfile = 0;
guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname] = uiputfile('*.mat','Select a location to store the project',....
  sprintf('Project_%s.mat',handles.params.behaviors.name));

if fname == 0; return; end;

handles.outfile = fullfile(pname,fname);
save(handles.outfile,'-struct',handles.params);

guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in pushbutton_copy.
function pushbutton_copy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fnames,pname] = uigetfile('*.mat','Select a project to copy settings from');
if fnames == 0; return; end;

dstr = {'Select the parameters to copy.','Use Control/Command click to select multiple entries'};
[sel,ok] = listdlg('PromptString',dstr,'ListSize',[350 60],'Name','Parameters to copy',...
  'ListString',{'Behavior and File Names','Window Features','Classifier files used as input'});

if ok == 0, return; end

origparams = load(fullfile(pname,fname));

if any(sel==1)
  handles.params.file = origparams.file;
  handles.params.behavior.names = origparams.behavior.names;
  listndx = find(handles.params.file.featureconfigfile,handles.list);
  set(handles.featureconfigpopup,'Value',listndx);
end

if any(sel == 2)
  handles.params.windowfeatures = origparams.params.windowfeatures;
end

if any(sel == 3)
  handles.params.scoresinput = origparams.params.scoresinput;
end

guidata(hObject,handles);


% --- Executes on button press in togglebutton_advanced.
function togglebutton_advanced_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_advanced
if get(hObject,'Value') ==1
  updatePosition(handles,'basic');
else
  updatePosition(handles,'advanced');
end

% --- Executes on button press in pushbutton_perframe.
function pushbutton_perframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_perframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
