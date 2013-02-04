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

% Last Modified by GUIDE v2.5 16-Jan-2013 15:37:04

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


% -------------------------------------------------------------------------
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

[new, ...
 figureJLabel, ...
 figureJLabelEditFiles, ...
 configParams] = ...
   myparse(varargin,...
           'new',[], ...
           'figureJLabel',[],...
           'figureJLabelEditFiles',[],...
           'configParams',[]);

handles.new=new;  
  % true iff we are setting up a new project, as opposed to editing an
  % existing one
handles.figureJLabel=figureJLabel;
handles.figureJLabelEditFiles=figureJLabelEditFiles;
handles.configParams = configParams;

% Change a few things so they still work well on Mac
adjustColorsIfMac(hObject);

curPos = get(handles.figureProjectSetup,'Position');
tablePos = get(handles.config_table,'Position');
reducedWidth = tablePos(1)-15;
reducedHeight = curPos(4);
handles.advancedSize = curPos(3:4);
handles.basicSize = [reducedWidth reducedHeight];
handles.mode = 'advanced';
set(handles.togglebutton_advanced,'Value',0);
handles = updatePosition(handles,'basic');

handles = initFeatureconfig(handles);
if ~isempty(configParams)
  handles.configParams = configParams;
  handles = addversion(handles);
  if isfield(configParams,'windowfeatures')
      handles.configParams.windowfeatures = configParams.windowfeatures;
  elseif isfield(configParams.file,'featureparamfilename')
    [windowfeaturesparams,windowfeaturescellparams,basicFeatureTable,featureWindowSize] = ...
      ReadPerFrameParams(configParams.file.featureparamfilename,configParams.file.featureconfigfile);
    handles.configParams.windowfeatures.windowfeaturesparams = windowfeaturesparams;
    handles.configParams.windowfeatures.windowfeaturescellparams = windowfeaturescellparams;
    handles.configParams.windowfeatures.basicFeatureTable = basicFeatureTable;
    handles.configParams.windowfeatures.featureWindowSize = featureWindowSize;
    handles.configParams.file = rmfield(handles.configParams.file,'featureparamfilename');
  else
    uiwait(warndlg(['The selected configuration file does not have any '...
      'window features. Not copying the window features']));
  end
  
  if ~isfield(configParams.file,'scorefilename')
    name_str = [sprintf('%s_',handles.configParams.behaviors.names{1:end-1}),handles.configParams.behaviors.names{end}];
    handles.configParams.file.scorefilename = sprintf('scores_%s',name_str);
  end
  
  if ~isfield(configParams,'scoresinput')
    handles.configParams.scoresinput = struct('classifierfile',{},'ts',{},'scorefilename',{});
  end
else
  handles = initParams(handles);
end

setConfigTable(handles);

% Change the window title to "New..." under some circumstance
if isempty(figureJLabelEditFiles)
  % this means ProjectSetup() was called directly from JLabel
  if new
    set(hObject,'name','New...');
  end
end

% Make invisible some controls if called directly from JLabel
if isempty(figureJLabelEditFiles)
  % this means ProjectSetup() was called directly from JLabel
  set(handles.textLabelFileName,'visible','off');
  set(handles.textGTLabelFileName,'visible','off');
  set(handles.editlabelfilename,'visible','off');
  set(handles.editgtlabelfilename,'visible','off');
end

% Update handles structure
guidata(hObject, handles);

updateText(handles);

set(hObject,'Visible','on');

return


% -------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function ProjectSetup_OutputFcn(hObject, ~, handles)  %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
%varargout{2} = handles.outfile;
%delete(handles.figureProjectSetup);
return


% -------------------------------------------------------------------------
function handles = updatePosition(handles,mode)
if strcmp(handles.mode,mode), return; end
if strcmp(handles.mode,'advanced')
  curPos = get(handles.figureProjectSetup,'Position');
  set(handles.figureProjectSetup,'Position',[curPos(1:2) handles.basicSize]);
  set(handles.togglebutton_advanced,'Value',0);
else
  curPos = get(handles.figureProjectSetup,'Position');
  set(handles.figureProjectSetup,'Position',[curPos(1:2) handles.advancedSize]);
  set(handles.togglebutton_advanced,'Value',1);  
end
handles.mode = mode;
return


% -------------------------------------------------------------------------
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
return


% -------------------------------------------------------------------------
function handles = initParams(handles)
handles.configParams.behaviors.names = {};
handles.configParams.file.featureconfigfile = handles.list.(handles.animal_types{1}).file;
handles.configParams.behaviors.type = handles.list.(handles.animal_types{1}).animal;
handles.configParams.file.perframedir = 'perframe';
handles.configParams.file.clipsdir = 'clips';
handles.configParams.scoresinput = struct('classifierfile',{},'ts',{},'scorefilename',{});
handles.configParams.windowfeatures = struct;
handles.configParams.behaviors.labelcolors = [0.7,0,0,0,0,0.7];
handles.configParams.behaviors.unknowncolor = [0,0,0];
handles.configParams.trx.colormap = 'jet';
handles.configParams.trx.colormap_multiplier = 0.7;
handles.configParams.trx.extra_linestyle = '-';
handles.configParams.trx.extra_marker = '.';
handles.configParams.trx.extra_markersize = 12;
handles.configParams.labels.colormap = 'line';
handles.configParams.labels.linewidth = 3;
handles.configParams.file.labelfilename = '';
handles.configParams.file.gt_labelfilename = '';
handles.configParams.file.scorefilename = '';
handles.configParams.file.trxfilename = '';
handles.configParams.file.moviefilename = '';
handles = addversion(handles);
handles.configParams.scoresinput = struct('classifierfile',{},'ts',{},'scorefilename',{});


% -------------------------------------------------------------------------
function updateText(handles)
% Copies the parameters in the configParams structure to the gui.

if isfield(handles.configParams.behaviors,'names');
  names = handles.configParams.behaviors.names;
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
  if isfield(handles.configParams.file,curf);
    set(handles.(curbox),'String',handles.configParams.file.(curf));
  else
    set(handles.(curbox),'String','');
  end
end

alltypes = fieldnames(handles.list);
listndx = [];
for andx = 1:numel(alltypes)
  curfname = handles.list.(alltypes{andx}).file;
  if strcmp(handles.configParams.file.featureconfigfile,curfname);
    listndx = andx;
  end
end
if ~isempty(listndx),
  set(handles.featureconfigpopup,'Value',listndx);
else
  uiwait(warndlg(['The feature config file used (', handles.configParams.file.featureconfigfile,...
    ') does not exist.']));
  return;
end

clist = {handles.configParams.scoresinput(:).classifierfile};
if isempty(clist),
  set(handles.listbox_inputscores,'String',{});
else
  set(handles.listbox_inputscores,'String',clist,'Value',1);
end


% -------------------------------------------------------------------------
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
      set(handles.(curbox),'String',handles.configParams.file.(curf));
      continue;
  end

  
  if ~isempty(str)
    handles.configParams.file.(curf) = str;
  end
end

setConfigTable(handles);


% -------------------------------------------------------------------------
function setConfigTable(handles)
configParams = handles.configParams;
fields2remove = {'featureparamlist','windowfeatures','scoresinput'};
for ndx = 1:numel(fields2remove)
  if isfield(configParams,fields2remove{ndx}),
    configParams = rmfield(configParams,fields2remove{ndx});
  end
  
end
data = GetParamsAsTable(configParams);
set(handles.config_table,'Data',data);


% -------------------------------------------------------------------------
function data = GetParamsAsTable(configParams)
data = addToList(configParams,{},'');
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


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function handles = addversion(handles)
if ~isfield(handles.configParams,'ver')
  vid = fopen('version.txt','r');
  vv = textscan(vid,'%s');
  fclose(vid);
  handles.configParams.ver = vv{1};
end


% -------------------------------------------------------------------------
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
handles.configParams.behaviors.names = name;
handles.configParams.file.moviefilename = 'movie.ufmf';
handles.configParams.file.trxfilename = 'registered_trx.mat';
handles.configParams.file.labelfilename = sprintf('label_%s.mat',name_str);
handles.configParams.file.gt_labelfilename = sprintf('gt_label_%s.mat',name_str);
handles.configParams.file.scorefilename = sprintf('scores_%s.mat',name_str);
updateText(handles);
guidata(hObject,handles);
setConfigTable(handles);


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
% --- Executes on selection change in featureconfigpopup.
function featureconfigpopup_Callback(hObject, eventdata, handles)
% hObject    handle to featureconfigpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns featureconfigpopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from featureconfigpopup

contents = cellstr(get(hObject,'String'));
curtype = contents{get(hObject,'Value')};
handles.configParams.file.featureconfigfile = handles.list.(curtype).file;
handles.configParams.behaviors.type = handles.list.(curtype).animal;
guidata(hObject,handles);

setConfigTable(handles);


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function editlabelfilename_Callback(hObject, eventdata, handles) %#ok<DEFNU>
% hObject    handle to editlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editlabelfilename as text
%        str2double(get(hObject,'String')) returns contents of editlabelfilename as a double
handles = updateParams(handles);
guidata(hObject,handles);


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function editgtlabelfilename_Callback(hObject, eventdata, handles)
% hObject    handle to editgtlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editgtlabelfilename as text
%        str2double(get(hObject,'String')) returns contents of editgtlabelfilename as a double
handles = updateParams(handles);
guidata(hObject,handles);


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function editscorefilename_Callback(hObject, eventdata, handles)
% hObject    handle to editscorefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editscorefilename as text
%        str2double(get(hObject,'String')) returns contents of editscorefilename as a double
handles = updateParams(handles);
guidata(hObject,handles);


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function editmoviefilename_Callback(hObject, eventdata, handles)
% hObject    handle to editmoviefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmoviefilename as text
%        str2double(get(hObject,'String')) returns contents of editmoviefilename as a double
handles = updateParams(handles);
guidata(hObject,handles);


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function edittrxfilename_Callback(hObject, eventdata, handles)
% hObject    handle to edittrxfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittrxfilename as text
%        str2double(get(hObject,'String')) returns contents of edittrxfilename as a double
handles = updateParams(handles);
guidata(hObject,handles);


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
% --- Executes on butuiresume(handles.figureProjectSetup);ton press in pushbutton_setfeatures.
function pushbutton_setfeatures_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setfeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_copyfeatures.
function pushbutton_copyfeatures_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copyfeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% -------------------------------------------------------------------------
% --- Executes on selection change in listbox_inputscores.
function listbox_inputscores_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_inputscores (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_inputscores contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_inputscores


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
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

handles.configParams.scoresinput(end+1) = curs;
guidata(hObject,handles);
updateText(handles);


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_removelist.
function pushbutton_removelist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_removelist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curndx = get(handles.listbox_inputscores,'Value');
if isempty(curndx), return; end
handles.configParams.scoresinput(curndx) = [];
guidata(hObject,handles);
updateText(handles);


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.outfile = 0;
guidata(hObject,handles);
delete(gcbf);              
%uiresume(handles.figureProjectSetup);
return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, ~, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.figureJLabelEditFiles)
  % this means ProjectSetup() was called directly from JLabel
  configParams=handles.configParams;
  %rmfield(configParams.file,'featureconfigfile');  
    % should uncomment the above eventually, but not now, ALT Jan 24, 2013
  configParams.file=rmfield(configParams.file,'labelfilename');
  configParams.file=rmfield(configParams.file,'gt_labelfilename');
  JLabel('projectSetupDone', ...
         handles.figureJLabel, ...
         configParams, ...
         handles.new);
else
  % this means ProjectSetup() was called from JLabelEditFiles
  JLabelEditFiles('projectSetupDone', ...
                  handles.figureJLabelEditFiles, ...
                  handles.configParams, ...
                  handles.new);
end                
delete(gcbf);              
return


% -------------------------------------------------------------------------
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
  handles.configParams.file.featureconfigfile = origparams.file.featureconfigfile;
end


if ismember('Behavior Name and File Names',sellist)
  origfeatureconfigfile = handles.configParams.file.featureconfigfile;
  handles.configParams.file = origparams.file;
  handles.configParams.file.featureconfigfile = origfeatureconfigfile;
  if iscell(origparams.behaviors.names),
    handles.configParams.behaviors.names = origparams.behaviors.names;
  else
    handles.configParams.behaviors.names = {origparams.behaviors.names};
  end
  if ~isfield(origparams.file,'scorefilename')
     name = handles.configParams.behaviors.names;
     name_str = [sprintf('%s_',name{1:end-1}),name{end}];
     handles.configParams.file.scorefilename =  sprintf('scores_%s.mat',name_str);
  end
end

if ismember('Window Features',sellist),
  if isfield(origparams,'windowfeatures')
    
    if ~strcmp(handles.configParams.file.featureconfigfile,origparams.file.featureconfigfile),
      res = questdlg(['Target type are not the same for the current project '...
        'and the original project. Are you sure you want to import the window features?'],...
        'Import Window features','Yes','No','No');
      if strcmp(res,'Yes')
        handles.configParams.windowfeatures = origparams.windowfeatures;
      end
    else
      handles.configParams.windowfeatures = origparams.windowfeatures;
    end
  elseif isfield(origparams.file,'featureparamfilename')
    uiwait(warndlg(['The selected configuration file does not have any '...
      'window features, but points to a file that does have the window features. '...
      'Loading the window features from the referenced file:' origparams.file.featureparamfilename]));
    [windowfeaturesparams,windowfeaturescellparams,basicFeatureTable,featureWindowSize] = ...
      ReadPerFrameParams(origparams.file.featureparamfilename,handles.configParams.file.featureconfigfile); 

    handles.configParams.windowfeatures.windowfeaturesparams = windowfeaturesparams;
    handles.configParams.windowfeatures.windowfeaturescellparams = windowfeaturescellparams;
    handles.configParams.windowfeatures.basicFeatureTable = basicFeatureTable;
    handles.configParams.windowfeatures.featureWindowSize = featureWindowSize;
    if isfield(handles.configParams.file,'featureparamfilename'),
        handles.configParams.file = rmfield(handles.configParams.file,'featureparamfilename');
    end
  else
    uiwait(warndlg(['The selected configuration file does not have any '...
      'window features. Not copying the window features']));
  end
end

if ismember('List of perframe features',sellist)
  
  if isfield(origparams,'featureparamlist')
    
    if ~strcmp(handles.configParams.file.featureconfigfile,origparams.file.featureconfigfile),
      res = questdlg(['Target type are not the same for the current project '...
        'and the original project. Are you sure you want to import the list of perframe features?'],...
        'Import list of perframe features','Yes','No','No');
      if strcmp(res,'Yes')
        handles.configParams.featureparamlist = origparams.featureparamlist;
      end
    else
      handles.configParams.featureparamlist= origparams.featureparamlist;
      
    end
  end
end

if ismember('Classifier files used as input', sellist);
  if isfield(origparams,'scoresinput')
    
    handles.configParams.scoresinput = origparams.scoresinput;
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
    'perframe.configParams.fov',...
    'perframe.configParams.nbodylengths_near',...
    'perframe.configParams.thetafil',...
    'perframe.configParams.max_dnose2ell_anglerange',...
    'perframe.landmark_params.arena_center_mm_x',...
    'perframe.landmark_params.arena_center_mm_y',...
    'perframe.landmark_params.arena_radius_mm',...
    'perframe.landmark_params.arena_type'};
  
  for str = adv_params(:)',
    try %#ok<TRYNC>
      eval(sprintf('handles.configParams.%s = origparams.%s;',str{1},str{1}));
    end
  end

end

handles = addversion(handles);
guidata(hObject,handles);
updateText(handles);
setConfigTable(handles);


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_perframe.
function pushbutton_perframe_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_perframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[allpflist,selected,missing] = GetAllPerframeList(handles.configParams);
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
  handles.configParams.featureparamlist = tstruct;
  
end

guidata(hObject,handles);


% -------------------------------------------------------------------------
function [allPfList selected missing]= GetAllPerframeList(configparams)
featureconfigfile = configparams.file.featureconfigfile;
configParams = ReadXMLParams(featureconfigfile);
allPfList = fieldnames(configParams.perframe);
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


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function handles = RemoveConfig(handles,name)

[fpath,lastfield] = splitext(name);
if isempty(lastfield)
  handles.configParams = rmfield(handles.configParams,fpath);
else
  evalStr = sprintf(...
    'handles.configParams.%s = rmfield(handles.configParams.%s,lastfield(2:end));',...
    fpath,fpath);
  eval(evalStr);
end


% -------------------------------------------------------------------------
function handles = EditConfigName(handles,oldName,newName)
eval_str = sprintf(...
  'value = handles.configParams.%s;',...
  oldName);
eval(eval_str);
handles = AddConfig(handles,newName,value);
handles = RemoveConfig(handles,oldName);


% -------------------------------------------------------------------------
function handles = EditConfigValue(handles,name,value) %#ok<INUSD>
eval_str = sprintf('handles.configParams.%s = value;',name);
eval(eval_str);


% -------------------------------------------------------------------------
function handles = AddConfig(handles,name,value)

if ischar(value) && ~isempty(str2num(value)) %#ok<ST2NM>
  value = str2num(value); %#ok<NASGU,ST2NM>
end

iname = fliplr(name);
curstruct = handles.configParams;
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

eval(sprintf('handles.configParams.%s = value;',name));


% -------------------------------------------------------------------------
% --- Executes when user attempts to close figureProjectSetup.
function figureProjectSetup_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figureProjectSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
%uiresume(hObject);
pushbutton_cancel_Callback(hObject, eventdata, handles)

% -------------------------------------------------------------------------
