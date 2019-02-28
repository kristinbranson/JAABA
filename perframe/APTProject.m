function varargout = APTProject(varargin)
% APTPROJECT MATLAB code for APTProject.fig
%      APTPROJECT, by itself, creates a new APTPROJECT or raises the existing
%      singleton*.
%
%      H = APTPROJECT returns the handle to a new APTPROJECT or the handle to
%      the existing singleton*.
%
%      APTPROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APTPROJECT.M with the given input arguments.
%
%      APTPROJECT('Property','Value',...) creates a new APTPROJECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before APTProject_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to APTProject_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help APTProject

% Last Modified by GUIDE v2.5 27-Feb-2019 16:21:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @APTProject_OpeningFcn, ...
                   'gui_OutputFcn',  @APTProject_OutputFcn, ...
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


% --- Executes just before APTProject is made visible.
function APTProject_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to APTProject (see VARARGIN)

% Choose default command line output for APTProject
lbl_file = myparse(varargin,'lbl_file','');
handles.lbl_file = lbl_file;

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Checks
if isempty(lbl_file)
  % close the gui if there is no label file.
  % the gui is closed in outputfcn.
  return
end

try
  apt = load(lbl_file,'-mat');
catch ME
  errmsg = sprintf('Could not load label file %s',ME.message);
  errordlg(errmsg, 'Unable to load label file');
  return
end
  
if isempty(apt.movieFilesAll)
  errordlg('APT Label file should have at least one movie','Empty APT project');
  return;
end

%
handles.apt = apt;
handles.n_pts = apt.cfg.NumLabelPoints;
handles.n_view = apt.cfg.NumViews;
handles.has_trx = apt.projectHasTrx;
if handles.has_trx
  handles.use_theta = apt.preProcParams.TargetCrop.AlignUsingTrxTheta;
  sz = apt.preProcParams.TargetCrop.Radius;
  handles.sz = [sz,sz];
elseif apt.cropIsCropMode
  % dont use theta if there is no trx file. 
  % might need updating later.
  handles.use_theta = False;
  assert(false, 'Update size when using crops')
  handles.sz = apt.movieFilesAllCropInfo;
else
  handles.use_theta = False;
  ht = apt.movieInfoAll{1}.info.max_height;
  wd = apt.movieInfoAll{1}.info.max_width;
  handles.sz = [ht,wd];
end

[featureList,xmlList] = getFeatureLexiconListsFromXML();
stFeatures = cellfun(@xmlParamsIsST,xmlList);
featureList = featureList(~stFeatures);
xmlList = xmlList(~stFeatures);
handles.trxList = {featureList, xmlList};

handles.global_ftrs = {'velmag','dvelmag','dphi'};
handles.body_ftrs = {'velmag','dvelmag','dphi',...
        'dist_center', 'ddist_center'};
handles.theta_ftrs = {'x','y','dx','dy','dtheta', ...
          'sin','cos'};
handles.pair_ftrs = {'pair_x','pair_y','pair_dist','pair_ddist',...
          'pair_cos','pair_sin','pair_dangle','pair_area','pair_dvelmag',...
          'pair_u1','pair_v1','pair_u2','pair_v2'};
handles.triad_ftrs = {'triad_cos','triad_sin','triad_dangle',...
          'triad_area','triad_darea','triad_dlen'};

handles = initialize(handles);
guidata(hObject,handles);

% UIWAIT makes APTProject wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = APTProject_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
if isempty(handles.lbl_file)
  figure1_CloseRequestFcn(hObject,eventdata,handles);
end

function handles = initialize(handles)
if handles.has_trx
  set(handles.trx_pop,'String',handles.trxList{1})
  set(handles.trx_pop,'Value',1);
else
  set(handles.trx_pop,'String',{''})
  set(handdles.trx_pop,'Value',0)
end
set(handles.trx_pop,'Enable',onoffif(handles.has_trx));

set(handles.body_checkbox,'Value',handles.has_trx);
set(handles.theta_checkbox,'Value',~handles.use_theta);
set(handles.theta_checkbox,'Enable',onoffif(handles.has_trx));
set(handles.global_checkbox,'Value',true);

% load a reference image and point locations.
apt = handles.apt;
has_ref_img = false;
for ndx = 1:size(apt.movieFilesAll,1)
  has_label = false;
  cur_pts = SparseLabelArray.full(apt.labeledpos{ndx});
  tndx_to_use = 1;
  locs = [];
  for tndx = 1:size(cur_pts,4)
    frm = find(~isnan(cur_pts(1,1,:,tndx)),1);
    if isempty(frm)
      continue;
    else
      has_label = true;
      tndx_to_use = tndx;
      locs = cur_pts(:,:,frm,tndx);
      break;
    end
  end
  if ~has_label
    continue;
  end
  has_ref_img = true;
  mov_file = unMacroise(apt.movieFilesAll{ndx,1},apt);
  rfn = get_readframe_fcn(mov_file);
  trx = load(getTrxFile(apt.trxFilesAll{ndx},mov_file),'-mat');
  trx = trx.trx(tndx_to_use);
  img = rfn(frm);
  off = trx.off;
  x = trx.x(frm+off);
  y = trx.y(frm+off);
  theta = trx(tndx_to_use).theta(frm+off);
  [img, aff_mat] = CropImAroundTrx(img,x,y,theta,handles.sz(2),handles.sz(1));
  ll = locs*aff_mat(1:2,1:2);
  ll(:,2) = ll(:,2)+aff_mat(3,2) + handles.sz(1);
  ll(:,1) = ll(:,1)+aff_mat(3,1) + handles.sz(2);
  locs = ll;
  break
end
if has_ref_img
  axes(handles.axes1);
  imshow(img);
  hold on;
  scatter(ll(:,1),ll(:,2),'.');
  hold off;
  handles.im = img;
  handles.locs = ll;
else
  handles.im = [];
  handles.locs = [];
end
axis(handles.axes1,'off');

function filename = unMacroise(filename,apt)
macros = fieldnames(apt.projMacros); 
for fndx = 1:numel(macros)
  f = macros{fndx};
  filename = strrep(filename,sprintf('$%s',f),apt.projMacros.(f));
end

function trxfile = getTrxFile(trxfile,movfile)
movdir = fileparts(movfile);
trxfile = strrep(trxfile,'$movdir',movdir);


function onoff = onoffif(in)
if in
  onoff = 'on';
else
  onoff = 'off';
end

% --- Executes on button press in trx_checkbox.
function trx_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to trx_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trx_checkbox


% --- Executes on button press in body_checkbox.
function body_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to body_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of body_checkbox


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5



function edit_pair_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pair as text
%        str2double(get(hObject,'String')) returns contents of edit_pair as a double
str = get(hObject,'String');
pair_str = strsplit(str,',');
pairs = {};
for ndx = 1:numel(pair_str)
  cur_str = strtrim(pair_str{ndx});
  cur_pair = strsplit(cur_str);
  if numel(cur_pair)~=2
    errordlg('Incorrect format for specifying pairs.');
    return;
  end
  pts = [];
  for pndx = 1:numel(cur_pair)
    cur_p = str2double(cur_pair{pndx});
    if round(cur_p) ~= cur_p || isnan(cur_p)
      errordlg('Please specify points as integer');
    end
    pts(end+1) = cur_p;
  end
  pairs{end+1} = pts;
end
handles.pairs = pairs;

img = handles.im;
ll = handles.locs;
axes(handles.axes1);
imshow(img);
hold on;
scatter(ll(:,1),ll(:,2),'.');
for pndx = 1:numel(handles.pairs)
  p = ll(handles.pairs{pndx},:);
  plot(p(:,1),p(:,2));
end
hold off;
axis(handles.axes1,'off');
  
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_pair_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_triad_Callback(hObject, eventdata, handles)
% hObject    handle to edit_triad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_triad as text
%        str2double(get(hObject,'String')) returns contents of edit_triad as a double
str = get(hObject,'String');
pair_str = strsplit(str,',');
triads = {};
for ndx = 1:numel(pair_str)
  cur_str = strtrim(pair_str{ndx});
  cur_pair = strsplit(cur_str);
  if numel(cur_pair)~=3
    errordlg('Incorrect format for specifying triads.');
    return;
  end
  pts = [];
  for pndx = 1:numel(cur_pair)
    cur_p = str2double(cur_pair{pndx});
    if round(cur_p) ~= cur_p || isnan(cur_p)
      errordlg('Please specify points as integer');
    end
    pts(end+1) = cur_p;
  end
  triads{end+1} = pts;
end
handles.triads = triads;
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function edit_triad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_triad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_custom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_custom as text
%        str2double(get(hObject,'String')) returns contents of edit_custom as a double
str = get(hObject,'String');
pair_str = strsplit(str,',');
custom_list = {};
for ndx = 1:numel(pair_str)
  custom_list = strtrim(pair_str{ndx});
end
handles.custom_list = custom_list;
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function edit_custom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trx_pop.
function trx_pop_Callback(hObject, eventdata, handles)
% hObject    handle to trx_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trx_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trx_pop


% --- Executes during object creation, after setting all properties.
function trx_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trx_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in global_checkbox.
function global_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to global_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of global_checkbox


% --- Executes on button press in theta_checkbox.
function theta_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to theta_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of theta_checkbox


% --- Executes on button press in done_pb.
function done_pb_Callback(hObject, eventdata, handles)
% hObject    handle to done_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in cancel_pb.
function cancel_pb_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);
