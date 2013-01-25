function varargout = PrepareJAABAData(varargin)
% PREPAREJAABADATA: Convert data from several tracking programs into 
% a format compatible with JAABA. 
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


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PrepareJAABAData_OpeningFcn, ...
                   'gui_OutputFcn',  @PrepareJAABAData_OutputFcn, ...
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


% --- Executes just before PrepareJAABAData is made visible.
function PrepareJAABAData_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PrepareJAABAData (see VARARGIN)

if ~isdeployed,
  SetUpJAABAPath;
end

handles = InitializeData(handles);
handles = InitializeGUI(handles);

buttons = findall(hObject,'Style','pushbutton');
for ndx = 1:numel(buttons)
  %SetButtonImage(buttons(ndx));
  adjustNonLabelButtonColor(buttons(ndx));
end

buttons = findall(hObject,'Style','togglebutton');
for ndx = 1:numel(buttons)
  %SetButtonImage(buttons(ndx));
  adjustNonLabelButtonColor(buttons(ndx));
end

if ismac, % On mac change the foreground color to black.
  allpopups = findall(hObject,'Style','popup');
  set(allpopups,'ForegroundColor',[0 0 0]);
end

guidata(hObject,handles);

% UIWAIT makes PrepareJAABAData wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PrepareJAABAData_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout = {hObject};

function handles = InitializeGUI(handles)

handles = CreateInputDataTypeRadioButtons(handles);

% set value for input data type
handles.InputDataTypeIndex = find(strcmpi(handles.InputDataTypeNames,handles.InputDataType),1);
if isempty(handles.InputDataTypeIndex),
  handles.InputDataTypeIndex = 1;
  handles.InputDataType = handles.InputDataTypeNames{1};
end
set(handles.radiobutton_InputDataTypes(handles.InputDataTypeIndex),'Value',1);

handles = UpdateInputDataType(handles);

handles = UpdateOutputFilesTable(handles);

% options
set(handles.checkbox_softlink,'Value',handles.SoftLinkFiles);
set(handles.checkbox_fliplr,'Value',handles.fliplr);
set(handles.checkbox_flipud,'Value',handles.flipud);
set(handles.checkbox_transpose,'Value',handles.dotransposeimage);
set(handles.edit_fps,'String',num2str(handles.fps));
set(handles.edit_pxpermm,'String',num2str(handles.pxpermm));
set(handles.checkbox_OverRideFPS,'Value',handles.OverRideFPS);
set(handles.checkbox_OverRideArena,'Value',handles.OverRideArena);
set(handles.edit_CropFirstFrame,'String',num2str(handles.CropFirstFrame));
set(handles.edit_CropEndFrame,'String',num2str(handles.CropEndFrame));

% arena parameters
UpdateArenaParameters(handles);
%set(handles.checkbox_ComputeArenaParameters,'Enable','off');

function UpdateArenaSize(handles)

if strcmpi(handles.ArenaType,'None'),
  set([handles.text_centerx,handles.edit_centerx,...
    handles.text_centery,handles.edit_centery,...
    handles.text_arenasize1,handles.edit_arenasize1,...
    handles.text_arenasize2,handles.edit_arenasize2],'Visible','off');
  %set(handles.checkbox_ComputeArenaParameters,'Visible','off');
  
elseif strcmpi(handles.ArenaType,'Circle'),
  set([handles.text_centerx,handles.edit_centerx,...
    handles.text_centery,handles.edit_centery,...
    handles.text_arenasize1,handles.edit_arenasize1],'Visible',handles.visible_edit_arena);
  set([handles.text_arenasize2,handles.edit_arenasize2],'Visible','off');
  set(handles.text_arenasize1,'String','Arena radius (px)');
  set(handles.edit_arenasize1,'String',num2str(handles.ArenaRadius));

else
  set([handles.text_centerx,handles.edit_centerx,...
    handles.text_centery,handles.edit_centery,...
    handles.text_arenasize1,handles.edit_arenasize1,...
    handles.text_arenasize2,handles.edit_arenasize2],'Visible',handles.visible_edit_arena);
  set(handles.text_arenasize1,'String','Arena width (px)');
  set(handles.edit_arenasize1,'String',num2str(handles.ArenaWidth));
  set(handles.edit_arenasize2,'String',num2str(handles.ArenaHeight));
end

if strcmpi(handles.visible_edit_pxpermm,'off') && strcmpi(handles.ArenaType,'None'),
  set(handles.pushbutton_ComputeArenaParameters,'Enable','off');
else
  set(handles.pushbutton_ComputeArenaParameters,'Enable','on');
end


function handles = UpdateInputDataType(handles)

InputDataType = handles.InputDataTypes.(handles.InputDataType);
InputFiles = handles.InputFiles.(handles.InputDataType);
ninputfiles = numel(InputDataType.files);
data = cell(ninputfiles+1,3);
data{1,1} = 'Video file';
data{1,2} = handles.InputVideoFile;
data{1,3} = ~InputDataType.videorequired || CheckInputFile(handles,1,handles.InputVideoFile);
for i = 1:ninputfiles,
  data{i+1,1} = InputDataType.files(i).name;
  data{i+1,2} = InputFiles{i};
  data{i+1,3} = ~InputDataType.files(i).required;
  data{i+1,3} = ~InputDataType.files(i).required || CheckInputFile(handles,i+1,InputFiles{i});
end
set(handles.uitable_InputFiles,'Data',data);

hedit_arena = [handles.text_arenatype,handles.popupmenu_arenatype,...
  handles.text_centerx,handles.edit_centerx,...
  handles.text_centery,handles.edit_centery,...
  handles.text_arenasize1,handles.edit_arenasize1,...
  handles.text_arenasize2,handles.edit_arenasize2];

handles.visible_edit_fps = 'on';
handles.visible_override_fps = 'on';
handles.visible_edit_pxpermm = 'on';
handles.visible_edit_arena = 'on';
handles.visible_override_arena = 'on';
handles.visible_compute_arena = 'on';
handles.visible_read_arena = 'on';
handles.visible_read_fps = 'on';

if strcmpi(InputDataType.writefps,'no'),
  handles.visible_edit_fps = 'off';
end
if strcmpi(InputDataType.writefps,'no') || strcmpi(InputDataType.readfps,'no'),
  handles.visible_override_fps = 'off';
end
if strcmpi(InputDataType.writepxpermm,'no'),
  handles.visible_edit_pxpermm = 'off';
end
if strcmpi(InputDataType.writearena,'no'),
  handles.visible_edit_arena = 'off';
end
if (strcmpi(InputDataType.writearena,'no')&&strcmpi(InputDataType.writepxpermm,'no')) || ...
    (strcmpi(InputDataType.readarena,'no')&&strcmpi(InputDataType.readpxpermm,'no')) || ...
    (strcmpi(InputDataType.readarena,'no')&&strcmpi(InputDataType.writepxpermm,'no')) || ...
    (strcmpi(InputDataType.writearena,'no')&&strcmpi(InputDataType.readpxpermm,'no')),
  handles.visible_override_arena = 'off';
end
if (strcmpi(InputDataType.writearena,'no')&&strcmpi(InputDataType.writepxpermm,'no')),
  handles.visible_compute_arena = 'off';
end
if (strcmpi(InputDataType.writearena,'no')&&strcmpi(InputDataType.writepxpermm,'no')) || ...
    (strcmpi(InputDataType.readarena,'no')&&strcmpi(InputDataType.readpxpermm,'no')) || ...
    (strcmpi(InputDataType.readarena,'no')&&strcmpi(InputDataType.writepxpermm,'no')) || ...
    (strcmpi(InputDataType.writearena,'no')&&strcmpi(InputDataType.readpxpermm,'no')),
  handles.visible_read_arena = 'off';
end
if strcmpi(InputDataType.writefps,'no'),
  handles.visible_read_fps = 'off';
end

set([handles.text_fps,handles.edit_fps],'Visible',handles.visible_edit_fps);
set([handles.text_OverRideFPS,handles.checkbox_OverRideFPS],'Visible',handles.visible_override_fps);
set([handles.text_pxpermm,handles.edit_pxpermm],'Visible',handles.visible_edit_pxpermm);
set([handles.text_arenatype,handles.popupmenu_arenatype,...
  handles.text_centerx,handles.edit_centerx,...
  handles.text_centery,handles.edit_centery,...
  handles.text_arenasize1,handles.edit_arenasize1,...
  handles.text_arenasize2,handles.edit_arenasize2],...
  'Visible',handles.visible_edit_arena);
set([handles.text_OverRideArena,handles.checkbox_OverRideArena],'Visible',handles.visible_override_arena);
set(handles.pushbutton_ComputeArenaParameters,'Visible',handles.visible_compute_arena);
set(handles.pushbutton_ReadArenaParameters,'Visible',handles.visible_read_arena);
set(handles.pushbutton_ReadFPS,'Visible',handles.visible_read_fps);

function handles = UpdateOutputFilesTable(handles)

data = get(handles.uitable_OutputFiles,'data');
data{1,2} = handles.ExperimentDirectory;
data{2,2} = handles.moviefilestr;
data{3,2} = handles.trxfilestr;
data{4,2} = handles.perframedirstr;
set(handles.uitable_OutputFiles,'data',data);

function handles = CreateInputDataTypeRadioButtons(handles)

pos1 = get(handles.radiobutton_InputDataType1,'Position');
xoff = pos1(1);
yoff = 1-pos1(2)-pos1(4);
width = pos1(3);

% where should the buttons go
height = (1-2*yoff)/handles.NInputDataTypes;
y = linspace(1-yoff-height,yoff,handles.NInputDataTypes);

fns = {'BackgroundColor','Callback','CData','Enable','FontAngle','FontName','FontSize','FontUnits','FontWeight',...
  'ForegroundColor','HorizontalAlignment','TooltipString','Units','UserData','Visible'};

props = get(handles.radiobutton_InputDataType1);
props = rmfield(props,setdiff(fieldnames(props),fns));

handles.radiobutton_InputDataTypes = nan(1,handles.NInputDataTypes);

% first radiobutton already made
set(handles.radiobutton_InputDataType1,...
  'String',handles.InputDataTypes.(handles.InputDataTypeNames{1}).name,...
  'Position',[xoff,y(1),width,height],...
  'Value',0,...
  'UserData',1);
handles.radiobutton_InputDataTypes(1) = handles.radiobutton_InputDataType1;

for i = 2:handles.NInputDataTypes,
  handles.radiobutton_InputDataTypes(i) = uicontrol('Parent',handles.uipanel_InputDataType,'Style','radiobutton');
  set(handles.radiobutton_InputDataTypes(i),props);
  set(handles.radiobutton_InputDataTypes(i),...
    'String',handles.InputDataTypes.(handles.InputDataTypeNames{i}).name,...
    'Position',[xoff,y(i),width,height],...
    'Value',0,...
    'UserData',i);
end



function handles = InitializeData(handles)

% types of inputs that can be converted
handles.InputDataTypes = SetInputDataTypeParams();
handles.NInputDataTypes = numel(fieldnames(handles.InputDataTypes));
handles.InputDataTypeNames = fieldnames(handles.InputDataTypes);


% rc file
if isdeployed,
  handles.rcfilename = deployedRelative2Global('.PrepareJAABAData_rc.mat');
else
  p = fileparts(mfilename('fullpath'));
  handles.rcfilename = fullfile(p,'.PrepareJAABAData_rc.mat');
end

% set default values
handles = SetDefaultValues(handles);

% allowed extensions for videos
handles.VideoFileExts = {'*.ufmf';'*.avi';'*.fmf';'*.sbfmf';'*.seq';'*.mov';'*.mmf';'*.mp4';'*.*'};

% arena file types
handles.ArenaTypes = get(handles.popupmenu_arenatype,'String');

function handles = SetDefaultValues(handles)

handles.config_fns = {'InputDataType','SoftLinkFiles','flipud','fliplr','dotransposeimage',...
  'fps','pxpermm','OverRideFPS','OverRideArena','ArenaType','ArenaCenterX','ArenaCenterY',...
  'ArenaRadius','ArenaWidth','ArenaHeight',...
  'inputdir','ExperimentDirectory','moviefilestr',...
  'trxfilestr','perframedirstr','inputfilestrs','inputmoviefilestr'};

handles.InputDataType = 'Ctrax';
handles.SoftLinkFiles = false;
handles.flipud = false;
handles.fliplr = false;
handles.dotransposeimage = false;
handles.fps = 30;
handles.OverRideFPS = false;
handles.OverRideArena = false;
handles.pxpermm = 1;
handles.ArenaType = 'None';
handles.ArenaCenterX = 0;
handles.ArenaCenterY = 0;
handles.ArenaRadius = 123;
handles.ArenaWidth = 123;
handles.ArenaHeight = 123;
handles.CropFirstFrame = 1;
handles.CropEndFrame = inf;

handles.inputdir = '';
defaultExpDir = sprintf('Exp%s',datestr(now,'yyyymmddTHHMMSS'));
p = pwd;
handles.ExperimentDirectory = fullfile(p,defaultExpDir);
handles.moviefilestr = '';
handles.trxfilestr = 'trx.mat';
handles.perframedirstr = 'perframe';

handles.inputfilestrs = struct;
fns1 = fieldnames(handles.InputDataTypes);
for i = 1:numel(fns1),
  handles.inputfilestrs.(fns1{i})= cell(1,numel(handles.InputDataTypes.(fns1{i}).files));
  handles.InputFiles.(fns1{i}) = cell(1,numel(handles.InputDataTypes.(fns1{i}).files));
end
handles.inputmoviefilestr = '';
handles.InputVideoFile = fullfile(handles.inputdir,handles.inputmoviefilestr);

if exist(handles.rcfilename,'file'),
  [handles,success,msg] = LoadConfiguration(handles,handles.rcfilename);
  if ~success,
    warndlg(msg);
  end
end

% experiment directory default
p = fileparts(handles.ExperimentDirectory);
handles.ExperimentDirectory = fullfile(p,defaultExpDir);

% moviefilestr default
if isempty(handles.moviefilestr),
  if ~isempty(handles.inputmoviefilestr),
    [~,~,ext] = fileparts(handles.inputmoviefilestr);
    handles.moviefilestr = ['movie',ext];
  else
    handles.moviefilestr = 'movie.ufmf';
  end
end

% --- Executes on button press in checkbox_softlink.
function checkbox_softlink_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_softlink (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_softlink
handles.SoftLinkFiles = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in checkbox_fliplr.
function checkbox_fliplr_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_fliplr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_fliplr
handles.fliplr = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes on button press in checkbox_flipud.
function checkbox_flipud_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_flipud (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_flipud
handles.flipud = get(hObject,'Value');
guidata(hObject,handles);



function edit_pxpermm_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pxpermm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pxpermm as text
%        str2double(get(hObject,'String')) returns contents of edit_pxpermm as a double
pxpermm = str2double(get(hObject,'String'));
if isnan(pxpermm) || pxpermm <= 0,
  warndlg('pxpermm must be a positive number','Bad pxpermm value');
  set(hObject,'String',num2str(handles.pxpermm));
  return;
end
handles.pxpermm = pxpermm;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_pxpermm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pxpermm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_fps_Callback(hObject, eventdata, handles)
% hObject    handle to edit_fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_fps as text
%        str2double(get(hObject,'String')) returns contents of edit_fps as a double
fps = str2double(get(hObject,'String'));
if isnan(fps) || fps <= 0,
  warndlg('fps must be a positive number','Bad fps value');
  set(hObject,'String',num2str(handles.fps));
  return;
end
handles.fps = fps;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_fps_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_fps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_centerx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_centerx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_centerx as text
%        str2double(get(hObject,'String')) returns contents of edit_centerx as a double
CenterX = str2double(get(hObject,'String'));
if isnan(CenterX),
  warndlg('Arena center x-coordinate must be a number','Bad ArenaCenterX value');
  set(hObject,'String',num2str(handles.ArenaCenterX));
  return;
end
handles.ArenaCenterX = CenterX;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_centerx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_centerx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_centery_Callback(hObject, eventdata, handles)
% hObject    handle to edit_centery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_centery as text
%        str2double(get(hObject,'String')) returns contents of edit_centery as a double
CenterY = str2double(get(hObject,'String'));
if isnan(CenterY),
  warndlg('Arena center y-coordinate must be a number','Bad ArenaCenterY value');
  set(hObject,'String',num2str(handles.ArenaCenterY));
  return;
end
handles.ArenaCenterY = CenterY;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_centery_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_centery (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_arenatype.
function popupmenu_arenatype_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_arenatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_arenatype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_arenatype
v = get(hObject,'Value');
handles.ArenaType = handles.ArenaTypes{v};
UpdateArenaSize(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_arenatype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_arenatype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_arenasize1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_arenasize1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_arenasize1 as text
%        str2double(get(hObject,'String')) returns contents of edit_arenasize1 as a double
size1 = str2double(get(hObject,'String'));
if isnan(size1) || size1 <= 0,
  warndlg('Arena size must be a positive number','Bad siz value');
  if strcmpi(handles.ArenaType,'Circle'),
    set(hObject,'String',num2str(handles.ArenaRadius));
  else
    set(hObject,'String',num2str(handles.ArenaWidth));
  end
  return;
end
if strcmpi(handles.ArenaType,'Circle'),
  handles.ArenaRadius = size1;
else
  handles.ArenaWidth = size1;
end
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_arenasize1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_arenasize1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_arenasize2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_arenasize2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_arenasize2 as text
%        str2double(get(hObject,'String')) returns contents of edit_arenasize2 as a double
size2 = str2double(get(hObject,'String'));
if isnan(size2) || size2 <= 0,
  warndlg('Arena size must be a positive number','Bad siz value');
  set(hObject,'String',num2str(handles.ArenaHeight));
  return;
end
handles.ArenaHeight = size2;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_arenasize2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_arenasize2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_ComputeArenaParameters.
function pushbutton_ComputeArenaParameters_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ComputeArenaParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.InputVideoFile),
  warndlg('Input video file has not yet been set, can''t show a sample frame');
  return;
end

SetBusy(handles,sprintf('Opening video file %s',handles.InputVideoFile));
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(handles.InputVideoFile);

% choose a frame
t = round(nframes/2);
im = readframe(t);
hfig = 1000;
figure(hfig);
clf;
hax = gca;
if size(im,3) > 1,
  image(im);
else
  imagesc(im,'Parent',hax);
  colormap gray;
end
axis(hax,'image');
hold(hax,'on');
if fid > 1,
  fclose(fid);
end

ClearBusy(handles);

hinstr = nan;

getpxpermm = strcmpi(handles.visible_edit_pxpermm,'on');

switch lower(handles.ArenaType),
  case 'circle',
    title('Label circular arena wall');
    hinstr = msgbox({'Click to enter points on the circular arena wall.'
      'Use normal button clicks to add points to the '
      'polyline.  A shift-, right-, or double-click adds '
      'a final point and ends the polyline selection.  '
      'Pressing RETURN or ENTER ends the polyline '
      'selection without adding a final point.  Pressing '
      'BACKSPACE or DELETE removes the previously '
      'selected point from the polyline.'},...
      'Label Circular Arena Wall');
        
    while true,
      [xc,yc,radius,h] = fitcircle_manual(hax);
      if isempty(xc),
        return;
      end
      axis(hax,'image');
      res = questdlg('Are you happy with the results?','','Yes','Redo','Cancel','Yes');
      switch lower(res),
        case 'redo',
          if any(ishandle(h)),
            delete(h(ishandle(h)))
          end
          continue;
        case 'cancel',
          if ishandle(hinstr), delete(hinstr); end
          return;
        case 'yes',
          break;
      end
    end

    handles.ArenaCenterX = xc;
    handles.ArenaCenterY = yc;
    handles.ArenaRadius = radius;
    guidata(hObject,handles);
    
    UpdateArenaParameters(handles);
    
    % get diameter
    if getpxpermm,
      while true,
        res = inputdlg('Diameter of arena in millimeters','',1);
        if isempty(res),
          if ishandle(hinstr), delete(hinstr); end
          return;
        end
        diameter_mm = str2double(res);
        if isnan(diameter_mm),
          warndlg('Please enter a number');
          continue;
        end
        break;
      end
      
      handles.pxpermm = 2*radius/diameter_mm;
      UpdateArenaParameters(handles);
    
      guidata(hObject,handles);
    end
    
  case 'rectangle',
    title('Label rectangular arena wall');
    hinstr = msgbox({'Outline the rectangular arena wall.'
      'Use the mouse to click and drag the desired rectangle.'},...
      'Label Rectangular Arena Wall');
    while true,

      try
        rect = getrect(hax);
      catch  %#ok<CTCH>
        return;
      end
      if isempty(rect),
        return;
      end
      xc = rect(1)+rect(3)/2;
      yc = rect(2)+rect(4)/2;
      width = rect(3);
      height = rect(4);

      % display the calculated center
      h = [];
      h(1) = plot(xc,yc,'mx','LineWidth',2);
      h(2) = text(xc,yc,sprintf('  (%.1f, %.1f)',xc,yc),'Color','m','FontWeight','bold');

      % plot the rectangle
      h(3) = plot(xc+width*.5*[-1,-1,1,1,-1],yc+height*.5*[-1,1,1,-1,-1],'m-','LineWidth',2);

      message = sprintf('Width is %.1f pixels, height is %.1f pixels',width,height);
      h(4) = text(15,15,message,'Color','m','FontWeight','bold','BackgroundColor','k');
      
      axis(hax,'image');
      res = questdlg('Are you happy with the results?','','Yes','Redo','Cancel','Yes');
      switch lower(res),
        case 'redo',
          if any(ishandle(h)),
            delete(h(ishandle(h)));
          end
          continue;
        case 'cancel',
          if ishandle(hinstr), delete(hinstr); end
          return;
        case 'yes',
          break;
      end
    end

    handles.ArenaCenterX = xc;
    handles.ArenaCenterY = yc;
    handles.ArenaWidth = width;
    handles.ArenaHeight = height;
    guidata(hObject,handles);
    
    UpdateArenaParameters(handles);
    
    if getpxpermm,
      % get width in millimeters
      while true,
        res = inputdlg('Width of arena in millimeters','',1);
        if isempty(res),
          if ishandle(hinstr), delete(hinstr); end
          return;
        end
        width_mm = str2double(res);
        if isnan(width_mm),
          warndlg('Please enter a number');
          continue;
        end
        break;
      end
      
      handles.pxpermm = width/width_mm;
      UpdateArenaParameters(handles);
      
      guidata(hObject,handles);
    end
  case 'none',
    
    title('Draw a line of known length');
    hinstr = msgbox({'Draw a line of known length, e.g. '
      'between two landmark points. The line shown is '
      'draggable and resizable. Double-click on line when '
      'done.'},...
      'Draw a Line of Known Length');
    hline = imdistline(hax);
    api = iptgetapi(hline);
    fcn = makeConstrainToRectFcn('imline',get(hax,'XLim'),get(hax,'YLim'));
    api.setPositionConstraintFcn(fcn);
    try
      pos = wait(hline);
    catch
      return;
    end
    if isempty(pos),
      return;
    end
    x = pos(:,1);
    y = pos(:,2);
    linelength = sqrt(diff(x).^2 + diff(y).^2);
    plot(x,y,'m-x','LineWidth',2)
    
    % get width in millimeters
    while true,
      res = inputdlg('Length of entered line in millimeters','',1);
      if isempty(res),
        if ishandle(hinstr), delete(hinstr); end
        return;
      end
      linelength_mm = str2double(res);
      if isnan(linelength_mm),
        warndlg('Please enter a number');
        continue;
      end
      break;
    end
    
    handles.pxpermm = linelength/linelength_mm;
    UpdateArenaParameters(handles);
    
end
  
if ishandle(hinstr), 
  delete(hinstr); 
end

% --- Executes on button press in pushbutton_Convert.
function pushbutton_Convert_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Convert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% save configuration
SaveConfiguration(handles,handles.rcfilename);

% check that all the files are ok
indata = get(handles.uitable_InputFiles,'data');
isok = [indata{:,3}];
if any(~isok),
  s = [{'The following input files are not set correctly: '}
    indata(~isok,1)];
  warndlg(s,'Problem with input files');
  return;
end

% check that file extensions make sense
[success,msg] = CheckOutputFileExtensions(handles,...
  'moviefilestr',handles.moviefilestr,...
  'trxfilestr',handles.trxfilestr);
if ~success,
  res = questdlg([msg,{'Convert anyways?'}]);
  if ~strcmpi(res,'Yes'),
    return;
  end
end

InputDataType = handles.InputDataTypes.(handles.InputDataType);
inputfiles = handles.InputFiles.(handles.InputDataType);
args = [{InputDataType.files.code}
  inputfiles];

SetBusy(handles,sprintf('Converting to output directory %s',handles.ExperimentDirectory));

try

[success,msg] = Convert2JAABAWrapper(handles.InputDataType,...
  'inmoviefile',handles.InputVideoFile,...
  args{:},...
  'expdir',handles.ExperimentDirectory,...
  'moviefilestr',handles.moviefilestr,...
  'trxfilestr',handles.trxfilestr,...
  'perframedirstr',handles.perframedirstr,...
  'overridefps',handles.OverRideFPS,...
  'overridearena',handles.OverRideArena,...
  'dosoftlink',handles.SoftLinkFiles,...
  'fliplr',handles.fliplr,...
  'flipud',handles.flipud,...
  'dotransposeimage',handles.dotransposeimage,...
  'fps',handles.fps,...
  'pxpermm',handles.pxpermm,...
  'arenatype',handles.ArenaType,...
  'arenacenterx',handles.ArenaCenterX,...
  'arenacentery',handles.ArenaCenterY,...
  'arenaradius',handles.ArenaRadius,...
  'arenawidth',handles.ArenaWidth,...
  'arenaheight',handles.ArenaHeight,...
  'frameinterval',[handles.CropFirstFrame,handles.CropEndFrame]);

catch ME,
  errordlg(getReport(ME),'Error converting');
  ClearBusy(handles);
  return;
end

ClearBusy(handles);

if ~success,
  uiwait(warndlg(msg,'Problem converting'));
else
  if isempty(msg),
    res = questdlg('Conversion succesful!','Conversion successful','Go to folder','Done','Go to folder');
  else
    res = questdlg(msg,'Conversion successful!','Go to folder','Done','Go to folder');
  end
  if strcmpi(res,'Go to folder'),
    if ispc,
      winopen(handles.ExperimentDirectory);
    else
      if isdeployed,
        myweb_nocheck(handles.ExperimentDirectory);
      else
        web(handles.ExperimentDirectory,'-browser');
      end
    end
  end
end

% --- Executes on button press in pushbutton_Save.
function pushbutton_Save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uiputfile('*.mat','Save Configuration');
if ~ischar(filename),
  return;
end
[success,msg] = SaveConfiguration(handles,fullfile(pathname,filename));
if ~success,
  warndlg(msg);
end

% --- Executes on button press in pushbutton_Cancel.
function pushbutton_Cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

SaveConfiguration(handles,handles.rcfilename);
delete(handles.figure1);

% --- Executes when selected object is changed in uipanel_InputDataType.
function uipanel_InputDataType_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_InputDataType 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

handles.InputDataTypeIndex = get(eventdata.NewValue,'UserData');
handles.InputDataType = handles.InputDataTypeNames{handles.InputDataTypeIndex};
handles = UpdateInputDataType(handles);
guidata(hObject,handles);


% --- Executes on button press in checkbox_ComputeArenaParameters.
function checkbox_ComputeArenaParameters_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_ComputeArenaParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_ComputeArenaParameters


% --- Executes when selected cell(s) is changed in uitable_InputFiles.
function uitable_InputFiles_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable_InputFiles (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if isempty(eventdata.Indices),
  return;
end
row = eventdata.Indices(1);
col = eventdata.Indices(2);
if col ~= 2,
  return;
end
if row == 1,
  exts = handles.VideoFileExts;
  %description = 'Input video file';
  name = 'Video file';
  inputfile = handles.InputVideoFile;
  inputfilestr = handles.inputmoviefilestr;
else
  i = row - 1;
  InputDataType = handles.InputDataTypes.(handles.InputDataType);
  exts = InputDataType.files(i).exts;
  %description = InputDataType.files(i).description;
  name = InputDataType.files(i).name;
  inputfile = handles.InputFiles.(handles.InputDataType){i};
  inputfilestr = handles.inputfilestrs.(handles.InputDataType){i};
end

if row > 1 && InputDataType.files(i).multiplefiles > 0,
  if isempty(inputfile),
    inputdir = handles.inputdir;
  else
    inputdir = myfileparts(inputfile{end});
  end
  args = {'FilterSpec',inputdir,'Prompt',sprintf('Choose input %s',lower(name)),'Output','cell'};
  if ~isempty(inputfile) && ~iscell(inputfile),
    inputfile = {inputfile};
  end
  if ~isempty(inputfile),
    args(end+1:end+2) = {'Append',inputfile};
  end
  if ~isempty(exts),
    exts = ReorderExts(exts,inputfile);
    args(end+1:end+2) = {'Type',exts};
  end  
  inputfile = uipickfiles(args{:});
  if isnumeric(inputfile),
    return;
  end
  [path,filename] = myfileparts(inputfile{end});
else
  if isempty(inputfile),
    inputfile = fullfile(handles.inputdir,inputfilestr);
  end
  if row > 1 && InputDataType.files(i).isdir,
    [inputparentdir,inputname] = fileparts(inputfile);
    inputfile = uigetdir2(inputparentdir,sprintf('Choose input %s',lower(name)),inputname);
    if ~ischar(inputfile),
      return;
    end
    [path,filename] = myfileparts(inputfile);
  else
    
    exts = ReorderExts(exts,inputfile);
    [filename,path,filteridx] = uigetfile(exts,sprintf('Choose input %s',lower(name)),inputfile);
    if ~ischar(filename),
      return;
    end
    inputfile = fullfile(path,filename);
  end
end

handles.inputdir = path;
if row > 1,
  handles.inputfilestrs.(handles.InputDataType){i} = filename;
  handles.InputFiles.(handles.InputDataType){i} = inputfile;
else
  handles.inputmoviefilestr = filename;
  handles.InputVideoFile = inputfile;
end
[isok,msg] = CheckInputFile(handles,row,inputfile);
data = get(handles.uitable_InputFiles,'data');
if row > 1 && InputDataType.files(i).multiplefiles > 0,
  s = sprintf('%s, ',inputfile{:});
  s = s(1:end-2);
  data{row,2} = s;
else
  data{row,2} = inputfile;
end
data{row,3} = isok;
set(handles.uitable_InputFiles,'data',data);
if ~isok,
  warndlg(msg,sprintf('Problem with %s',name));
  guidata(hObject,handles);
  return;
end

guidata(hObject,handles);

function exts = ReorderExts(exts,inputfile)

% reorder exts
if numel(exts) > 1 && ~isempty(inputfile),
  if iscell(inputfile),
    [~,~,ext] = fileparts(inputfile{end});
  else
    [~,~,ext] = fileparts(inputfile);
  end
  if ~isempty(ext),
    exts = [{['*',ext]};setdiff(exts,{['*',ext]})];
  end
end

function [isok,msg] = CheckInputFile(handles,row,inputfile)

isok = true;
msg = '';

if iscell(inputfile),
  for i = 1:numel(inputfile),
    isok = isok && exist(inputfile{i},'file') > 0;
    if ~isok,
      msg = sprintf('File %s does not exist.',inputfile{i});
      break;
    end
  end
else
  isok = exist(inputfile,'file') > 0;
end
if ~isok,
  msg = sprintf('File %s does not exist.',inputfile);
end

% 
% if row == 1,
%   try
%     [readframe,~,fid,headerinfo] = get_readframe_fcn(inputfile);
%     readframe(1);
%     isok = true;
%     if fid > 1,
%       fclose(fid);
%     end
%   catch ME,
%     isok = false;
%     msg = getReport(ME);
%   end
%   return;
% end
% 
% switch handles.InputDataType,
%   case 'Ctrax',
%     InputDataType = handles.InputDataTypes.Ctrax;
%     inputfiletype = InputDataType.files(row-1).name;
%     switch inputfiletype,
%       case 'Trx mat file',
%         try
%           [trx,matname,succeeded] = load_tracks(inputfile);
%           if ~succeeded,
%             isok = false;
%             msg = sprintf('Could not load trx from mat file %s',inputfile);
%           else
%             isok = true;
%           end
%         catch ME,
%           isok = false;
%           msg = getReport(ME);
%         end
%       case 'Ann file',
%         % WORKING HERE
%     end
% end

% --- Executes when entered data in editable cell(s) in uitable_InputFiles.
function uitable_InputFiles_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_InputFiles (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
fprintf('Bye\n');


% --- Executes on button press in checkbox_OverRideFPS.
function checkbox_OverRideFPS_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_OverRideFPS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_OverRideFPS
handles.OverRideFPS = get(hObject,'Value');
guidata(hObject,handles);

% --- Executes on button press in pushbutton_ReadArenaParameters.
function pushbutton_ReadArenaParameters_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ReadArenaParameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

InputDataType = handles.InputDataTypes.(handles.InputDataType);

if strcmpi(InputDataType.readarena,'no') && strcmpi(InputDataType.readpxpermm,'no'),
  warndlg(sprintf('Cannot read arena parameters for data type %s',InputDataType.name));
  return;
end

SetBusy(handles,'Reading arena parameters...');

inputfiles = handles.InputFiles.(handles.InputDataType);
args = [{InputDataType.files.code}
  inputfiles];

[success,msg,...
  handles.ArenaType,handles.ArenaCenterX,handles.ArenaCenterY,...
  handles.ArenaRadius,handles.ArenaWidth,handles.ArenaHeight,...
  handles.pxpermm] = ...
  ReadArenaParameters_Wrapper(handles.InputDataType,...
  'inmoviefile',handles.InputVideoFile,...
  args{:},...
  'pxpermm',handles.pxpermm,...
  'arenatype',handles.ArenaType,...
  'arenacenterx',handles.ArenaCenterX,...
  'arenacentery',handles.ArenaCenterY,...
  'arenaradius',handles.ArenaRadius,...
  'arenawidth',handles.ArenaWidth,...
  'arenaheight',handles.ArenaHeight);

if ~success,
  warndlg(msg,'Could not read arena parameters');
  ClearBusy(handles);
  return;
end

UpdateArenaParameters(handles);

ClearBusy(handles);

msgbox(msg,'Read arena parameters');
guidata(hObject,handles);

function UpdateArenaParameters(handles)

i = find(strcmpi(handles.ArenaTypes,handles.ArenaType),1);
if isempty(i),
  warndlg(sprintf('Illegal value %s for arenatype',handles.ArenaType),'Bad arenatype');
  return;
end
set(handles.popupmenu_arenatype,'Value',i);
set(handles.edit_centerx,'String',num2str(handles.ArenaCenterX));
set(handles.edit_centery,'String',num2str(handles.ArenaCenterY));
set(handles.edit_pxpermm,'String',num2str(handles.pxpermm));

UpdateArenaSize(handles);

% --- Executes on button press in pushbutton_ReadFPS.
function pushbutton_ReadFPS_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ReadFPS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

InputDataType = handles.InputDataTypes.(handles.InputDataType);

% try to read from movie
if isempty(handles.InputVideoFile) && InputDataType.videorequired,
  warndlg('Input video file has not been set yet');
  return;
end

SetBusy(handles,'Reading FPS...');

readfpsfrom = '';

if ~isempty(handles.InputVideoFile),
  try 
    [readframe,nframes,fid,headerinfo] = get_readframe_fcn(handles.InputVideoFile);
    if isfield(headerinfo,'timestamps'),
      handles.fps = 1/median(diff(headerinfo.timestamps));
      readfpsfrom = 'inmoviefile';
    elseif isfield(headerinfo,'FrameRate'),
      handles.fps = headerinfo.FrameRate;
      readfpsfrom = 'inmoviefile';
    elseif isfield(headerinfo,'m_fFps'),
      handles.fps = headerinfo.m_fFps;
      readfpsfrom = 'inmoviefile';
    else 
      
      f0 = ceil(nframes*.1);
      f1 = floor(nframes*.9);
      if f1 <= f0,
        f0 = 1; f1 = nframes;
      end
      [~,t0] = readframe(f0);
      [~,t1] = readframe(f1);
      if ~isnan(t0) && ~isnan(t1),
        handles.fps = (f1-f0)/(t1-t0);
        readfpsfrom = 'inmoviefile';
      end
      
    end
    if fid > 1,
      fclose(fid);
    end
  catch %#ok<CTCH>
    readfpsfrom = '';
  end
end

if ~isempty(readfpsfrom),
  msg = 'Read fps from moviefile';
else
  args = [{InputDataType.files.code}
    handles.InputFiles.(handles.InputDataType)];
  [success,msg,handles.fps] = ReadFPS_Wrapper(handles.InputDataType,args{:});
  if ~success,
    warndlg(msg,'Could not read FPS');
    ClearBusy(handles);
    return;
  end
end

set(handles.edit_fps,'String',num2str(handles.fps));

ClearBusy(handles);

guidata(hObject,handles);

msgbox(msg);


% --- Executes when selected cell(s) is changed in uitable_OutputFiles.
function uitable_OutputFiles_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable_OutputFiles (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

if isempty(eventdata.Indices),
  return;
end
row = eventdata.Indices(1);
col = eventdata.Indices(2);
if col ~= 2 || row ~= 1,
  return;
end

[parentdir,expname] = myfileparts(handles.ExperimentDirectory);
if ~exist(parentdir,'dir'),
  parentdir = pwd;
end
expdir = uigetdir2(parentdir,'Choose experiment directory',expname);
if ~ischar(expdir),
  return;
end
handles.ExperimentDirectory = expdir;
data = get(hObject,'Data');
data{1,2} = handles.ExperimentDirectory;
set(hObject,'Data',data);
guidata(hObject,handles);


% --- Executes when entered data in editable cell(s) in uitable_OutputFiles.
function uitable_OutputFiles_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable_OutputFiles (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

if isempty(eventdata.Indices),
  return;
end
row = eventdata.Indices(1);
col = eventdata.Indices(2);
if col ~= 2,
  return;
end
data = get(hObject,'Data');
if isempty(eventdata.NewData),
  switch row,
    case 1,
      data{1,2} = handles.ExperimentDirectory;
    case 2,
      data{2,2} = handles.moviefilestr;
    case 3,
      data{3,2} = handles.trxfilestr;
    case 4,
      data{4,2} = handles.perframedirstr;
  end
  set(hObject,'Data',data);
  return;
end

if row == 1,
  % make sure experiment directory is legal
  [parentdir,expname] = myfileparts(eventdata.NewData);
  if ~isempty(parentdir) && ~exist(parentdir,'dir'),
    res = questdlg(sprintf('Parent directory %s does not exist. Create?',parentdir));
    if ~strcmpi(res,'Yes'),
      data{1,2} = handles.ExperimentDirectory;
      set(hObject,'Data',data);
      return;
    end
    [success,msg] = mkdir(parentdir);
    if ~success,
      uiwait(warndlg(msg));
      data{1,2} = handles.ExperimentDirectory;
      set(hObject,'Data',data);
      return;
    end
  end
  if ~IsNiceFileName(expname),
    res = questdlg(sprintf('Experiment name %s is not a great file name. Are you sure you want to use it?',expname));
    if ~strcmpi(res,'Yes'),
      data{1,2} = handles.ExperimentDirectory;
      set(hObject,'Data',data);
      return;
    end
  end
  handles.ExperimentDirectory = eventdata.NewData;
  guidata(hObject,handles);
  return;
else
  
  name = data{row,2};
  
  % check the input name
  if ~IsNiceFileName(name),
    res = questdlg(sprintf('%s is not a great file name. Are you sure you want to use it?',name));
  else
    res = 'Yes';
  end
  if strcmpi(res,'Yes'),
    if row == 2,
      [isokext,msg] = CheckOutputFileExtensions(handles,...
        'moviefilestr',name);
    elseif row == 3,
      [isokext,msg] = CheckOutputFileExtensions(handles,...
        'trxfilestr',name);
    else
      isokext = true;
    end
    if ~isokext,
      res = questdlg([msg,{'Do you want to use it anyway?'}]);
    end
  end
  if ~strcmpi(res,'Yes'),
    switch row,
      case 2,
        data{2,2} = handles.moviefilestr;
      case 3,
        data{3,2} = handles.trxfilestr;
      case 4,
        data{4,2} = handles.perframedirstr;
    end
    set(hObject,'Data',data);
    return;
  end
  switch row,
    case 2,
      handles.moviefilestr = name;
    case 3,
      handles.trxfilestr = name;
    case 4,
      handles.perframedirstr = name;
  end
  guidata(hObject,handles);

end

function [success,msg] = SaveConfiguration(handles,filename)

SetBusy(handles,sprintf('Saving configuration to file %s',filename));

rcdata = struct;
for i = 1:numel(handles.config_fns),
  fn = handles.config_fns{i};
  if isfield(handles,fn),
    rcdata.(fn) = handles.(fn);
  end
end

try
  save(filename,'-struct','rcdata');
catch ME,
  msg = getReport(ME);
  success = false;
  ClearBusy(handles);
  return;
end

success = true;
msg = '';

ClearBusy(handles);

function [handles,success,msg] = LoadConfiguration(handles,filename,ignorefns)

if nargin < 3,
  ignorefns = {};
end

rc = struct;
if exist(filename,'file'),
  try
    rc = load(filename);
  catch ME,
    msg = getReport(ME);
    success = false;
    return;
  end
end

fns = setdiff(handles.config_fns,ignorefns);

% fill in read-in stuff
for i = 1:numel(fns),
  fn = fns{i};
  if isfield(rc,fn),
    handles.(fn) = rc.(fn);
  end
end

if ismember('inputfilestrs',fns),
  fns0 = fieldnames(handles.inputfilestrs);
  fns1 = fieldnames(handles.InputDataTypes);
  idxremove = ~ismember(fns0,fns1);
  if any(idxremove),
    handles.inputfilestrs = rmfield(handles.inputfilestrs,fns0(idxremove));
  end
  idxadd = find(~ismember(fns1,fns0));
  for i = idxadd(:)',
    handles.inputfilestrs.(fns1{i})= cell(1,numel(handles.InputDataTypes.(fns1{i}).files));  
  end
end

success = true;
msg = '';

function SetBusy(handles,msg)

if nargin < 2,
  msg = 'Thinking...';
end

oldpointer = get(handles.figure1,'Pointer');
if strcmpi(oldpointer,'watch'),
  oldpointer = 'arrow';
end
setappdata(handles.figure1,'OldPointer',oldpointer);
set(handles.figure1,'Pointer','watch');

set(handles.text_Status,'String',msg,'ForegroundColor','m');
drawnow;

function ClearBusy(handles)

oldpointer = getappdata(handles.figure1,'OldPointer');
set(handles.figure1,'Pointer',oldpointer);
set(handles.text_Status,'String','Ready.','ForegroundColor','g');
drawnow;


% --- Executes on button press in pushbutton_Load.
function pushbutton_Load_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[filename,pathname] = uigetfile('*.mat','Load Configuration');
if ~ischar(filename),
  return;
end

ignorefns = {'ExperimentDirectory'};

[handles,success,msg] = LoadConfiguration(handles,fullfile(pathname,filename),ignorefns);
if ~success,
  warndlg(msg);
  return;
end

handles = UpdateGUI(handles);
guidata(hObject,handles);


function handles = UpdateGUI(handles)

% set value for input data type
handles.InputDataTypeIndex = find(strcmpi(handles.InputDataTypeNames,handles.InputDataType),1);
if isempty(handles.InputDataTypeIndex),
  handles.InputDataTypeIndex = 1;
  handles.InputDataType = handles.InputDataTypeNames{1};
end
set(handles.radiobutton_InputDataTypes(handles.InputDataTypeIndex),'Value',1);

handles = UpdateInputDataType(handles);

handles = UpdateOutputFilesTable(handles);

% options
set(handles.checkbox_softlink,'Value',handles.SoftLinkFiles);
set(handles.checkbox_fliplr,'Value',handles.fliplr);
set(handles.checkbox_flipud,'Value',handles.flipud);
set(handles.checkbox_transpose,'Value',handles.dotransposeimage);
set(handles.edit_fps,'String',num2str(handles.fps));
set(handles.edit_pxpermm,'String',num2str(handles.pxpermm));
set(handles.checkbox_OverRideFPS,'Value',handles.OverRideFPS);
set(handles.checkbox_OverRideArena,'Value',handles.OverRideArena);
set(handles.edit_CropFirstFrame,'String',num2str(handles.CropFirstFrame));
set(handles.edit_CropEndFrame,'String',num2str(handles.CropEndFrame));

% arena parameters
UpdateArenaParameters(handles);


% --- Executes on button press in checkbox_OverRideArena.
function checkbox_OverRideArena_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_OverRideArena (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_OverRideArena
handles.OverRideArena = get(hObject,'Value');
guidata(hObject,handles);


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
SaveConfiguration(handles,handles.rcfilename);
delete(handles.figure1);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function [success,msg] = CheckOutputFileExtensions(handles,varargin)

success = true;
msg = {};

[moviefilestr,trxfilestr] = myparse(varargin,...
  'moviefilestr','','trxfilestr','');

if ~isempty(handles.InputVideoFile) && ~isempty(moviefilestr),
  [~,~,extin] = fileparts(handles.InputVideoFile);
  [~,~,extout] = fileparts(moviefilestr);
  if ~isempty(extin) && ~isempty(extout) && ~strcmpi(extin,extout),
    success = false;
    msg{end+1} = sprintf('Output video file extension %s does not match input video file extension %s',extout,extin);
  end
end
if ~isempty(trxfilestr),
  [~,~,extout] = fileparts(trxfilestr);
  if ~strcmpi(extout,'.mat'),
    success = false;
    msg{end+1} = sprintf('Output trx file extension is ''%s'' and should be ''.mat''',extout);
  end
end


% --- Executes on button press in checkbox_transpose.
function checkbox_transpose_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_transpose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_transpose
handles.dotransposeimage = get(hObject,'Value');
guidata(hObject,handles);

function edit_CropFirstFrame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CropFirstFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CropFirstFrame as text
%        str2double(get(hObject,'String')) returns contents of edit_CropFirstFrame as a double
v = str2double(get(hObject,'String'));
if isnan(v) || v > handles.CropEndFrame,
  warndlg('First frame must be a number between 1 and end frame');
  set(hObject,'String',num2str(handles.CropFirstFrame));
  return;
end
handles.CropFirstFrame = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_CropFirstFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CropFirstFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_CropEndFrame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_CropEndFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_CropEndFrame as text
%        str2double(get(hObject,'String')) returns contents of edit_CropEndFrame as a double

v = str2double(get(hObject,'String'));
if isnan(v) || v < handles.CropFirstFrame,
  warndlg('End frame must be a number greater than or equal to first frame. Enter "Inf" if the end of the trajectory should not be cropped.');
  set(hObject,'String',num2str(handles.CropFirstFrame));
  return;
end
handles.CropEndFrame = v;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_CropEndFrame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_CropEndFrame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
