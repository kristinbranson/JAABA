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

% Last Modified by GUIDE v2.5 07-Nov-2012 12:56:26

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
function ProjectSetup_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProjectSetup (see VARARGIN)

% Choose default command line output for ProjectSetup
set(hObject,'Visible','off');

handles.output = hObject;
handles.defpath = pwd;

params = [];
[project_params, xml_file, mat_file, visible, outfile] = myparse(varargin,...
  'project_params',[],...
  'xml_file','',...
  'mat_file','',...
  'visible',true,...
  'outfile',0);

handles.outfile = outfile;


if ~isempty(project_params);
  params = project_params;
end
if ~isempty(xml_file),
  if ~exist(xml_file,'file'),
    errordlg('Cannot open the xml file');
    error('Cannot open the xml file');
  end
  params = ReadXMLParams(xml_file);
  [dpath,fpath,~] = fileparts(xml_file);
  handles.defpath = fullfile(dpath,[fpath '.mat']);
  if ~iscell(params.behaviors.names),
    params.behaviors.names = {params.behaviors.names};
  end
end
if ~isempty(mat_file),
  if ~exist(mat_file,'file')
    errordlg('Cannot open the mat file');
    error('Cannot open the mat file');
  end
  handles.defpath = mat_file;
  params = load(mat_file);
end

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

curPos = get(handles.figure1,'Position');
tablePos = get(handles.config_table,'Position');
reducedWidth = tablePos(1)-15;
reducedHeight = curPos(4);
handles.advancedSize = curPos(3:4);
handles.basicSize = [reducedWidth reducedHeight];
handles.mode = 'advanced';
set(handles.togglebutton_advanced,'Value',0);
handles = updatePosition(handles,'basic');

handles = initFeatureconfig(handles);
if ~isempty(params)
  handles.params = params;
  handles = addversion(handles);
  if isfield(params,'windowfeatures')
      handles.params.windowfeatures = params.windowfeatures;
  elseif isfield(params.file,'featureparamfilename')
    [windowfeaturesparams,windowfeaturescellparams,basicFeatureTable,featureWindowSize] = ...
      ReadPerFrameParams(params.file.featureparamfilename,params.file.featureconfigfile);
    handles.params.windowfeatures.windowfeaturesparams = windowfeaturesparams;
    handles.params.windowfeatures.windowfeaturescellparams = windowfeaturescellparams;
    handles.params.windowfeatures.basicFeatureTable = basicFeatureTable;
    handles.params.windowfeatures.featureWindowSize = featureWindowSize;
    handles.params.file = rmfield(handles.params.file,'featureparamfilename');
  else
    uiwait(warndlg(['The selected configuration file does not have any '...
      'window features. Not copying the window features']));
  end
  
  if ~isfield(params.file,'scorefilename')
    name_str = [sprintf('%s_',handles.params.behaviors.names{1:end-1}),handles.params.behaviors.names{end}];
    handles.params.file.scorefilename = sprintf('scores_%s',name_str);
  end
  
  if ~isfield(params,'scoresinput')
    handles.params.scoresinput = struct('classifierfile',{},'ts',{},'scorefilename',{});
  end
else
  handles = initParams(handles);
end

setConfigTable(handles);

% Update handles structure
guidata(hObject, handles);

updateText(handles);

if visible,
  set(hObject,'Visible','on');
end


% UIWAIT makes ProjectSetup wait for user response (see UIRESUME)
if ~ischar(handles.outfile),
  uiwait(handles.figure1);
else
  params2save = handles.params(1); %#ok<NASGU>
  save(handles.outfile,'-struct','params2save');
end



% --- Outputs from this function are returned to the command line.
function varargout = ProjectSetup_OutputFcn(hObject, ~, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles.outfile;
delete(handles.figure1);

function handles = updatePosition(handles,mode)
if strcmp(handles.mode,mode), return; end
if strcmp(handles.mode,'advanced')
  curPos = get(handles.figure1,'Position');
  set(handles.figure1,'Position',[curPos(1:2) handles.basicSize]);
  set(handles.togglebutton_advanced,'Value',0);
else
  curPos = get(handles.figure1,'Position');
  set(handles.figure1,'Position',[curPos(1:2) handles.advancedSize]);
  set(handles.togglebutton_advanced,'Value',1);  
end
handles.mode = mode;

function handles = initFeatureconfig(handles)
if isdeployed,
  filename = deployedRelative2Global('params/featureConfigList.xml');
else
  filename = 'featureConfigList.xml';
end
list = ReadXMLParams(filename);
animal_types = fieldnames(list);

handles.list = list;
handles.animal_types = animal_types;
set(handles.featureconfigpopup,'String',animal_types);


function handles = initParams(handles)
handles.params.behaviors.names = {};
handles.params.file.featureconfigfile = handles.list.(handles.animal_types{1}).file;
handles.params.behaviors.type = handles.list.(handles.animal_types{1}).animal;
handles.params.file.perframedir = 'perframe';
handles.params.file.clipsdir = 'clips';
handles.params.scoresinput = struct('classifierfile',{},'ts',{},'scorefilename',{});
handles.params.windowfeatures = struct;
handles.params.behaviors.labelcolors = [0.7,0,0,0,0,0.7];
handles.params.behaviors.unknowncolor = [0,0,0];
handles.params.trx.colormap = 'jet';
handles.params.trx.colormap_multiplier = 0.7;
handles.params.trx.extra_linestyle = '-';
handles.params.trx.extra_marker = '.';
handles.params.trx.extra_markersize = 12;
handles.params.labels.colormap = 'line';
handles.params.labels.linewidth = 3;
handles.params.file.labelfilename = '';
handles.params.file.gt_labelfilename = '';
handles.params.file.scorefilename = '';
handles.params.file.trxfilename = '';
handles.params.file.moviefilename = '';
handles = addversion(handles);
handles.params.scoresinput = struct('classifierfile',{},'ts',{},'scorefilename',{});

function updateText(handles)
% Copies the parameters in the params structure to the gui.

if isfield(handles.params.behaviors,'names');
  names = handles.params.behaviors.names;
  if numel(names)>0
    namestr = [sprintf('%s,',names{1:end-1}),names{end}];
  else
    namestr = '';
  end
  set(handles.editName,'String',namestr);
else
  set(handles.editName,'String','');
end

fnames = {'labelfilename','gt_labelfilename','scorefilename',...
  'moviefilename','trxfilename'};
boxnames = {'editlabelfilename','editgtlabelfilename','editscorefilename',...
  'editmoviefilename','edittrxfilename'};

for ndx = 1:numel(fnames)
  curf = fnames{ndx};
  curbox = boxnames{ndx};
  if isfield(handles.params.file,curf);
    set(handles.(curbox),'String',handles.params.file.(curf));
  else
    set(handles.(curbox),'String','');
  end
end

alltypes = fieldnames(handles.list);
listndx = [];
for andx = 1:numel(alltypes)
  curfname = handles.list.(alltypes{andx}).file;
  if strcmp(handles.params.file.featureconfigfile,curfname);
    listndx = andx;
  end
end
if ~isempty(listndx),
  set(handles.featureconfigpopup,'Value',listndx);
else
  uiwait(warndlg(['The feature config file used (', handles.params.file.featureconfigfile,...
    ') does not exist.']));
  return;
end

clist = {handles.params.scoresinput(:).classifierfile};
if isempty(clist),
  set(handles.listbox_inputscores,'String',{});
else
  set(handles.listbox_inputscores,'String',clist,'Value',1);
end


function handles = updateParams(handles)

fnames = {'labelfilename','gt_labelfilename','scorefilename',...
  'moviefilename','trxfilename'};
boxnames = {'editlabelfilename','editgtlabelfilename','editscorefilename',...
  'editmoviefilename','edittrxfilename'};

for ndx = 1:numel(fnames)
  curf = fnames{ndx};
  curbox = boxnames{ndx};
  str = strtrim(get(handles.(curbox),'String'));
  if ~isempty(str) && ~IsNiceFileName(str),
      uiwait(warndlg(sprintf(...
          ['The name specified for %s cannot have special characters.'...
          'Please use only alphanumeric characters and underscore'],curf)));
      set(handles.(curbox),'String',handles.params.file.(curf));
      continue;
  end

  
  if ~isempty(str)
    handles.params.file.(curf) = str;
  end
end

setConfigTable(handles);


function setConfigTable(handles)
params = handles.params;
fields2remove = {'featureparamlist','windowfeatures','scoresinput'};
for ndx = 1:numel(fields2remove)
  if isfield(params,fields2remove{ndx}),
    params = rmfield(params,fields2remove{ndx});
  end
  
end
data = GetParamsAsTable(params);
set(handles.config_table,'Data',data);


function data = GetParamsAsTable(configparams)
data = addToList(configparams,{},'');
idx = cellfun(@iscell,data(:,2));
if any(idx),
  for i = find(idx(:)'),
    if all(cellfun(@ischar,data{i,2})),
      data{i,2} = sprintf('%s,',data{i,2}{:});
      if numel(data{i,2})>0 && data{i,2}(end) == ',',
        data{i,2} = data{i,2}(1:end-1);
      end
    end
  end
end

if any(cellfun(@iscell,data(:,2))),
  data = {}; 
  return;
end


function list = addToList(curStruct,list,pathTillNow)
if isempty(fieldnames(curStruct)), return; end
fnames = fieldnames(curStruct);
for ndx = 1:numel(fnames)
  if isstruct(curStruct.(fnames{ndx})),
    list = addToList(curStruct.(fnames{ndx}),list,[pathTillNow fnames{ndx} '.']);
  else
    list{end+1,1} = [pathTillNow fnames{ndx}]; %#ok<AGROW>
    param = curStruct.(fnames{ndx});
    if isnumeric(param)
      q = num2str(param(1));
      for jj = 2:numel(param)
        q = [q ',' num2str(param(jj))]; %#ok<AGROW>
      end
      list{end,2} = q;
    else
      list{end,2} = param;
    end
  end
end

function handles = addversion(handles)
if ~isfield(handles.params,'ver')
  vid = fopen('version.txt','r');
  vv = textscan(vid,'%s');
  fclose(vid);
  handles.params.ver = vv{1};
end

function editName_Callback(hObject, eventdata, handles)
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editName as text
%        str2double(get(hObject,'String')) returns contents of editName as a double
name = get(hObject,'String');
if isempty(regexp(name,'^[a-zA-Z][\w_,]*$','once','start'));
   uiwait(warndlg(['The behavior name cannot have special characters.'...
       'Please use only alphanumeric characters and _']));
   return;
end
    
name = regexp(name,',','split');
name_str = [sprintf('%s_',name{1:end-1}),name{end}];
handles.params.behaviors.names = name;
handles.params.file.labelfilename = sprintf('label_%s.mat',name_str);
handles.params.file.gt_labelfilename = sprintf('gt_label_%s.mat',name_str);
handles.params.file.scorefilename = sprintf('scores_%s.mat',name_str);
updateText(handles);
guidata(hObject,handles);
setConfigTable(handles);


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

contents = cellstr(get(hObject,'String'));
curtype = contents{get(hObject,'Value')};
handles.params.file.featureconfigfile = handles.list.(curtype).file;
handles.params.behaviors.type = handles.list.(curtype).animal;
guidata(hObject,handles);

setConfigTable(handles);


% --- Executes during object creation, after setting all properties.
function featureconfigpopup_CreateFcn(hObject, ~, ~)
% hObject    handle to featureconfigpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editlabelfilename_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to editlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editlabelfilename as text
%        str2double(get(hObject,'String')) returns contents of editlabelfilename as a double
handles = updateParams(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function editlabelfilename_CreateFcn(hObject, eventdata, handles) %#ok<*DEFNU>
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
  uiwait(warndlg(['The selected file does not have all the required fields,'...
    ' Check to make sure it is a classifier file']));
  return;
end

if isfield(classifier,'scorefilename')
  scorefilename = classifier.scorefilename;
  [~,scorefilename,~] = fileparts(scorefilename);
  curs.scorefilename = scorefilename;
else
  configfile = classifier.configfilename;
  if exist(configfile,'file')
    configparams = ReadXMLParams(configfile);
    if isfield(configparams,'scorefilename');
      scorefilename = configparams.file.scorefilename;
      [~,scorefilename,~] = fileparts(scorefilename);
      curs.scorefilename = scorefilename;
    else
      curs.scorefilename = sprintf('scores_%s',configparams.behaviors.names);
    end
  else
    sname = inputdlg(['Cannot find the config file pointed by the classifier.'...
      'Input the file name used to stored the scores computed by the classifier'],...
      'Score file name');
    if isempty(sname), return; end
    curs.scorefilename = sname;
  end
end

handles.params.scoresinput(end+1) = curs;
guidata(hObject,handles);
updateText(handles);


% --- Executes on button press in pushbutton_removelist.
function pushbutton_removelist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curndx = get(handles.listbox_inputscores,'Value');
if isempty(curndx), return; end
handles.params.scoresinput(curndx) = [];
guidata(hObject,handles);
updateText(handles);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.outfile = 0;
guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname,pname] = uiputfile('*.mat','Select a location to store the project',....
  handles.defpath);

if fname == 0; return; end;

handles.outfile = fullfile(pname,fname);

if exist(handles.outfile,'file')
  [didback,msg] = copyfile(handles.outfile,[handles.outfile '~']);
  if ~didback,
    warning('Could not create backup of %s: %s',handles.outfile,msg);
  end
end

params2save = handles.params(1); %#ok<NASGU>
save(handles.outfile,'-struct','params2save');

guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in pushbutton_copy.
function pushbutton_copy_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fnames,pname] = uigetfile('*.mat','Select a project to copy settings from');
if fnames == 0; return; end;

list = {'Target Type','Behavior Name and File Names','Window Features',...
  'List of perframe features','Classifier files used as input',...
  'Advanced Parameters'};

dstr = {'Select the parameters to copy.','Use Control/Command click to select multiple entries'};
[sel,ok] = listdlg('PromptString',dstr,'ListSize',[350 120],'Name','Parameters to copy',...
  'ListString',list);
sellist = list(sel);

if ok == 0, return; end
[~,~,ext] = fileparts(fnames);
if strcmp(ext,'.xml')
  origparams = ReadXMLParams(fullfile(pname,fnames));
else
  origparams = load(fullfile(pname,fnames));
end

if ismember('Target Type',sellist)
  alltypes = fieldnames(handles.list);
  listndx = [];
  for andx = 1:numel(alltypes)
    curfname = handles.list.(alltypes{andx}).file;
    if strcmp(origparams.file.featureconfigfile,curfname); 
      listndx = andx;
    end
  end
  if ~isempty(listndx),
    set(handles.featureconfigpopup,'Value',listndx);
  else
    uiwait(warndlg(['The feature config file used (', origparams.file.featureconfigfile,...
      ') does not exist. Cannot import from the project']));
    return;
  end
  handles.params.file.featureconfigfile = origparams.file.featureconfigfile;
end


if ismember('Behavior Name and File Names',sellist)
  origfeatureconfigfile = handles.params.file.featureconfigfile;
  handles.params.file = origparams.file;
  handles.params.file.featureconfigfile = origfeatureconfigfile;
  if iscell(origparams.behaviors.names),
    handles.params.behaviors.names = origparams.behaviors.names;
  else
    handles.params.behaviors.names = {origparams.behaviors.names};
  end
  if ~isfield(origparams.file,'scorefilename')
     name = handles.params.behaviors.names;
     name_str = [sprintf('%s_',name{1:end-1}),name{end}];
     handles.params.file.scorefilename =  sprintf('scores_%s.mat',name_str);
  end
end

if ismember('Window Features',sellist),
  if isfield(origparams,'windowfeatures')
    
    if ~strcmp(handles.params.file.featureconfigfile,origparams.file.featureconfigfile),
      res = questdlg(['Target type are not the same for the current project '...
        'and the original project. Are you sure you want to import the window features?'],...
        'Import Window features','Yes','No','No');
      if strcmp(res,'Yes')
        handles.params.windowfeatures = origparams.windowfeatures;
      end
    else
      handles.params.windowfeatures = origparams.windowfeatures;
    end
  elseif isfield(origparams.file,'featureparamfilename')
    uiwait(warndlg(['The selected configuration file does not have any '...
      'window features, but points to a file that does have the window features. '...
      'Loading the window features from the referenced file:' origparams.file.featureparamfilename]));
    [windowfeaturesparams,windowfeaturescellparams,basicFeatureTable,featureWindowSize] = ...
      ReadPerFrameParams(origparams.file.featureparamfilename,handles.params.file.featureconfigfile); 

    handles.params.windowfeatures.windowfeaturesparams = windowfeaturesparams;
    handles.params.windowfeatures.windowfeaturescellparams = windowfeaturescellparams;
    handles.params.windowfeatures.basicFeatureTable = basicFeatureTable;
    handles.params.windowfeatures.featureWindowSize = featureWindowSize;
    handles.params.file = rmfield(handles.params.file,'featureparamfilename');    
  else
    uiwait(warndlg(['The selected configuration file does not have any '...
      'window features. Not copying the window features']));
  end
end

if ismember('List of perframe features',sellist)
  
  if isfield(origparams,'featureparamlist')
    
    if ~strcmp(handles.params.file.featureconfigfile,origparams.file.featureconfigfile),
      res = questdlg(['Target type are not the same for the current project '...
        'and the original project. Are you sure you want to import the list of perframe features?'],...
        'Import list of perframe features','Yes','No','No');
      if strcmp(res,'Yes')
        handles.params.featureparamlist = origparams.featureparamlist;
      end
    else
      handles.params.featureparamlist= origparams.featureparamlist;
      
    end
  end
end

if ismember('Classifier files used as input', sellist);
  if isfield(origparams,'scoresinput')
    
    handles.params.scoresinput = origparams.scoresinput;
  end
end

if ismember('Advanced Parameters',sellist),
  adv_params = {'behaviors.labelcolors',...
    'behaviors.unknowncolor',...
    'plot.trx.colormap',...
    'plot.trx.colormap_multiplier',...
    'plot.trx.extra_marker',...
    'plot.trx.extra_markersize',...
    'plot.trx.extra_linestyle',...
    'plot.labels.colormap',...
    'plot.labels.linewidth',...
    'perframe.params.fov',...
    'perframe.params.nbodylengths_near',...
    'perframe.params.thetafil',...
    'perframe.params.max_dnose2ell_anglerange',...
    'perframe.landmark_params.arena_center_mm_x',...
    'perframe.landmark_params.arena_center_mm_y',...
    'perframe.landmark_params.arena_radius_mm',...
    'perframe.landmark_params.arena_type'};
  
  for str = adv_params(:)',
    try %#ok<TRYNC>
      eval(sprintf('handles.params.%s = origparams.%s;',str{1},str{1}));
    end
  end

end

handles = addversion(handles);
guidata(hObject,handles);
updateText(handles);
setConfigTable(handles);

% --- Executes on button press in togglebutton_advanced.
function togglebutton_advanced_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_advanced
if get(hObject,'Value') ==1
  handles = updatePosition(handles,'advanced');
else
  handles = updatePosition(handles,'basic');
end
guidata(hObject,handles);

% --- Executes on button press in pushbutton_perframe.
function pushbutton_perframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_perframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[allpflist,selected,missing] = GetAllPerframeList(handles.params);
if ~isempty(missing),
   list = missing{1};
  for ndx = 2:numel(missing)
    list = sprintf('%s %s ',list,missing{ndx});
  end
  wstr = sprintf('Perframe feature(s) %s are defined in the project file%s\n%s',...
        list,...
        ' but are not defined for the current target type.',...
        'Ignoring them.');
  uiwait(warndlg(wstr));
 
end

[sel,ok] = listdlg('ListString',allpflist,...
  'InitialValue',find(selected),'Name','Selecte perframe features',...
  'PromptString',{'Control/Command click to','select/deselect perframe features'},...
  'ListSize',[250,700]);

if ok,
  tstruct = struct();
  for ndx = sel(:)'
    tstruct.(allpflist{ndx}) = 1;
  end
  handles.params.featureparamlist = tstruct;
  
end

guidata(hObject,handles);



function [allPfList selected missing]= GetAllPerframeList(configparams)
featureconfigfile = configparams.file.featureconfigfile;
params = ReadXMLParams(featureconfigfile);
allPfList = fieldnames(params.perframe);
selected = false(numel(allPfList),1);
missing = {};
if isfield(configparams,'featureparamlist');
  curpf = fieldnames(configparams.featureparamlist);
else
  curpf = {};
end
if ~isempty(curpf),
  for ndx = 1:numel(curpf),
    allndx = find(strcmp(curpf{ndx},allPfList));
    if isempty(allndx)
      missing{end+1} = curpf{ndx};
    else
      selected(allndx) = true;
    end
  end
else
  missing = {};
  selected = true(numel(allPfList,1));
end



% --- Executes on button press in pushbutton_addadvanced.
function pushbutton_addadvanced_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_addadvanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

in = inputdlg({'Configuration Parameter Name','Configuration Parameter Value'});
if isempty(in) || numel(in) < 2
  return;
end

handles = AddConfig(handles,in{1},in{2});
setConfigTable(handles);
updateText(handles);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_removeadvanced.
function pushbutton_removeadvanced_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removeadvanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
data = get(handles.config_table,'Data');

jUIScrollPane = findjobj(handles.config_table);
jUITable = jUIScrollPane.getViewport.getView;
allndx = jUITable.getSelectedRows + 1;
if numel(allndx)==1 && allndx <1, return, end

for ndx = allndx(:)'
  handles = RemoveConfig(handles,data{ndx,1});
end
setConfigTable(handles);
updateText(handles);
guidata(hObject,handles);


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
  handles = EditConfigValue(handles,data{ndx,1},data{ndx,2});
else
  handles = EditConfigName(handles,eventdata.PreviousData,eventdata.NewData);
end
updateText(handles);
guidata(hObject,handles);

function handles = RemoveConfig(handles,name)

[fpath,lastfield] = splitext(name);
if isempty(lastfield)
  handles.params = rmfield(handles.params,fpath);
else
  evalStr = sprintf(...
    'handles.params.%s = rmfield(handles.params.%s,lastfield(2:end));',...
    fpath,fpath);
  eval(evalStr);
end


function handles = EditConfigName(handles,oldName,newName)
eval_str = sprintf(...
  'value = handles.params.%s;',...
  oldName);
eval(eval_str);
handles = AddConfig(handles,newName,value);
handles = RemoveConfig(handles,oldName);

function handles = EditConfigValue(handles,name,value) %#ok<INUSD>
eval_str = sprintf('handles.params.%s = value;',name);
eval(eval_str);

function handles = AddConfig(handles,name,value)

if ischar(value) && ~isempty(str2num(value)) %#ok<ST2NM>
  value = str2num(value); %#ok<NASGU,ST2NM>
end

iname = fliplr(name);
curstruct = handles.params;
while true,
  [iname,lastfield] = splitext(iname);
  if isempty(lastfield)
    fexist = isfield(curstruct,iname);
    break;
  else
    fexist = isfield(curstruct,fliplr(iname(2:end)));
    if ~fexist, break;    end
    curstruct = curstruct.(fliplir(iname(2:end)));
  end
end

eval(sprintf('handles.params.%s = value;',name));


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
uiresume(hObject);
