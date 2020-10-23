function varargout = spaceTimeOptions(varargin)
% SPACETIMEOPTIONS MATLAB code for spaceTimeOptions.fig
%      SPACETIMEOPTIONS, by itself, creates a new SPACETIMEOPTIONS or raises the existing
%      singleton*.
%
%      H = SPACETIMEOPTIONS returns the handle to a new SPACETIMEOPTIONS or the handle to
%      the existing singleton*.
%
%      SPACETIMEOPTIONS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPACETIMEOPTIONS.M with the given input arguments.
%
%      SPACETIMEOPTIONS('Property','Value',...) creates a new SPACETIMEOPTIONS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spaceTimeOptions_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spaceTimeOptions_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spaceTimeOptions

% Last Modified by GUIDE v2.5 22-Oct-2020 04:46:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spaceTimeOptions_OpeningFcn, ...
                   'gui_OutputFcn',  @spaceTimeOptions_OutputFcn, ...
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


% --- Executes just before spaceTimeOptions is made visible.
function spaceTimeOptions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spaceTimeOptions (see VARARGIN)

% Choose default command line output for spaceTimeOptions
handles.output = hObject;

[mov_file,trx_file,prev_st] = myparse(varargin,'mov_file','','trx_file','','prev_st',[]);
handles.mov_file = mov_file;
handles.trx_file = trx_file;

[rfn,nframes,fid,headerinfo] = get_readframe_fcn(mov_file);
handles.rfn = rfn;
handles.nframes = nframes;
ex_frame = rfn(1);
fr_sz = size(ex_frame);
handles.width = size(ex_frame,2);
handles.height = size(ex_frame,1);
handles.fid = fid;
handles.headerinfo = headerinfo;

if ~isempty(trx_file)
    trx = load_tracks(trx_file);
    handles.trx = trx;
    handles.n_trx = length(trx);
    handles.cur_fly = 1;
    handles.frame_num = randsample(trx(handles.cur_fly).nframes,1) - trx(handles.cur_fly).off;
else
    % if not trx then dummy tracks at the center of the frame.
    trx = struct('x',repmat(fr_sz(2)/2,[1,nframes]),...
                'y',repmat(fr_sz(1)/2,[1,nframes]),...
                'theta',repmat(90,[1,nframes]),...
                'firstframe',1,...
                'endframe',nframes,...
                'off',0,...
                'nframes',nframes);
                
    handles.trx = trx;
    handles.n_trx = 1;
    handles.cur_fly = 1;
    handles.frame_num = randsample(nframes,1);
end


if isempty(prev_st)
    prev_st = getSTParams();
end
if ~isfield(prev_st,'x_center')
  prev_st.x_center = 0;
  prev_st.y_center = 0;
  prev_st.imsz_display = prev_st.psize*max(prev_st.npatches_x,prev_st.npatches_y);
end
handles.params = prev_st;

updateGUI(handles);
updatePlots(handles);

handles.is_canceled = true;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spaceTimeOptions wait for user response (see UIRESUME)
uiwait(handles.figure1);


function updateGUI(handles)
params = handles.params;
set(handles.edit_patch,'String',sprintf('%d',params.psize));
set(handles.edit_horz,'String',sprintf('%d',params.npatches_x));
set(handles.edit_vert,'String',sprintf('%d',params.npatches_y));
set(handles.edit_frame,'String',sprintf('%d',handles.frame_num));
set(handles.edit_trx,'String',sprintf('%d',handles.cur_fly));
set(handles.edit_x_center,'String',sprintf('%d',params.x_center));
set(handles.edit_y_center,'String',sprintf('%d',params.y_center));
method_ndx = find(strcmp(handles.params.methods,handles.params.cur_method));
set(handles.menu_methods,'String',handles.params.method_names);
set(handles.menu_methods,'Value',method_ndx);

function handles = updateimsz(handles)
params = handles.params;
grid_sz_x = params.npatches_x*params.psize;
grid_sz_y = params.npatches_y*params.psize;
x_c = params.x_center; y_c = params.y_center;
max_sz_x = max(x_c+grid_sz_x/2,-x_c+grid_sz_x/2);
max_sz_y = max(y_c+grid_sz_y/2,-y_c+grid_sz_y/2);
max_sz = max(max_sz_x,max_sz_y)*2;
params.imsz_display = max(params.imsz_display,max_sz);
handles.params = params;
set(handles.edit_imsz,'String',sprintf('%d',params.imsz_display));


function handles = updatePlots(handles)

set(handles.figure1, 'pointer', 'watch')
handles = updateimsz(handles);
trx = handles.trx;
cur_fly = handles.cur_fly;
frame_num = handles.frame_num;
rfn = handles.rfn;
trx_fn = frame_num + trx(cur_fly).off;
cur_frame = rfn(frame_num);
x = trx(cur_fly).x(trx_fn);
y = trx(cur_fly).y(trx_fn);
theta = trx(cur_fly).theta(trx_fn);
psize = handles.params.psize;
nbins = handles.params.nbins;

ncells_x = handles.params.npatches_x;
ncells_y = handles.params.npatches_y;
x_c = handles.params.x_center;
y_c = handles.params.y_center;
winwidth = handles.params.imsz_display;
winheight = handles.params.imsz_display;
hsz = handles.params.imsz_display/2;
grid_hsz_x = ncells_x*psize/2;
grid_hsz_y = ncells_y*psize/2;
handles.cur_im = CropImAroundTrx(cur_frame,x,y,theta-pi/2,...
    (winwidth-1)/2,(winheight-1)/2);
cla(handles.axes_hog);
cla(handles.axes_hof);

ftrs = genSTFeatures(handles.rfn,handles.headerinfo,frame_num,frame_num+1,...
    handles.trx,handles.params.is_stationary,handles.params.cur_method,handles.params);

% hog ftrs
scale = 6;
im1 = handles.cur_im;

axes(handles.axes_hog);
him = imshow(imresize(im1,scale));
axis image;
colormap gray;
hold on;
axis off;

colors = hsv(nbins);
colors = colors([ (end/2+1):end 1:end/2],:);

[nr,nc,~] = size(im1);
% maxv2 = max(H(:));
maxv2 = 0.03;
wd = 0.5;
hogpatch = [wd wd -wd -wd wd;-psize psize psize -psize -psize]/2;
h = [];
H = ftrs.hogftrs{cur_fly}(:,:,:,1);
bincenters = linspace(0,pi,nbins+1);
bincenters = bincenters(1:nbins);
dt = mean(diff(bincenters));
binedges = [bincenters(1)-dt/2,(bincenters(1:end-1)+bincenters(2:end))/2,bincenters(end)+dt/2];

for xi = 1:ncells_x
  cx = (psize/2 + (xi-1)*psize - grid_hsz_x + x_c + hsz )*scale+ 1 ;
  if cx+psize/2 > nc*scale
    break;
  end
  for yi = 1:ncells_y
    cy = (psize/2 + (yi-1)*psize - grid_hsz_y + y_c + hsz)*scale+ 1 ;
    if cy+psize/2 > nr*scale
      break;
    end
    
    for bini = 1:nbins
      tmp = bincenters(bini);
      curpatch = [cos(tmp) -sin(tmp); sin(tmp) cos(tmp)]*hogpatch;
      xcurr = cx + curpatch(1,:)*scale;
      ycurr = cy + curpatch(2,:)*scale;
      h(yi,xi,bini) = patch(xcurr,ycurr,colors(bini,:),...
          'LineStyle','none',...
          'FaceAlpha',min(1,H(yi,xi,bini)/maxv2));
    end    
  end
end

for xi = 0:ncells_x
    plot([(xi*psize+x_c + hsz- grid_hsz_x)*scale,(xi*psize+x_c+hsz-grid_hsz_x)*scale],...
      [ (y_c+hsz-grid_hsz_y)*scale,(y_c+ncells_y*psize+hsz-grid_hsz_y)*scale],...
      'Color',[0.3,0.3,0.3,0.2]);
end
for yi = 0:ncells_y
    plot([(x_c+hsz-grid_hsz_x)*scale,(x_c+ncells_x*psize+hsz-grid_hsz_x)*scale],...
      [ (yi*psize+y_c+hsz-grid_hsz_y)*scale,(yi*psize+y_c+hsz-grid_hsz_y)*scale],...
      'Color',[0.3,0.3,0.3,0.2]);
end
scatter(hsz*scale,hsz*scale,'marker','*');
scatter((x_c+hsz)*scale,(y_c+hsz)*scale,'marker','o');


% Flow (hof) ftrs
if handles.params.is_stationary && frame_num < handles.nframes,
  ny = trx(cur_fly).y(trx_fn+1);
  nx = trx(cur_fly).x(trx_fn+1);
  ntheta = trx(cur_fly).theta(trx_fn+1);
  n_fr = rfn(frame_num+1);
else
  if frame_num == handles.nframes
      n_fr = cur_frame;
  else
      n_fr = rfn(frame_num+1);
  end
  ny = y; nx = x; ntheta = theta;
end

im2 = CropImAroundTrx(n_fr,nx,ny,ntheta-pi/2,(winwidth-1)/2,(winheight-1)/2);
axes(handles.axes_hof);
him = imshowpair(imresize(im1,scale),imresize(im2,scale),'Scaling','none');
axis image;
colormap gray;
hold on;

colors = hsv(nbins);
[nr,nc,~] = size(im1);
% maxv2 = max(F(:));
maxv2 = 1;

F = ftrs.flowftrs{cur_fly}(:,:,:,1);
h = [];
bincenters = linspace(0,pi,nbins+1);
bincenters = bincenters(1:nbins);

bincenters2 = bincenters*2;
dt = mean(diff(bincenters2));
binedges2 = [bincenters2(1)-dt/2,(bincenters2(1:end-1)+bincenters2(2:end))/2,bincenters2(end)+dt/2];

for xi = 1:ncells_x
  cx = (psize/2  + (xi-1)*psize  + x_c +hsz-grid_hsz_x)*scale+ 1;
  if cx+psize/2 > nc*scale
    break;
  end
  for yi = 1:ncells_y
    cy = (psize/2 + (yi-1)*psize + y_c+hsz-grid_hsz_y)*scale+ 1;
    if cy+psize/2 > nr*scale
      break;
    end
    
    for bini = 1:nbins
      tmp = linspace(binedges2(bini),binedges2(bini+1),20);
      xcurr = cx + [0,psize/2*cos(tmp),0]*scale;
      ycurr = cy + [0,psize/2*sin(tmp),0]*scale;
      h(yi,xi,bini) = patch(xcurr,ycurr,colors(bini,:),'LineStyle','none','FaceAlpha',min(1,F(yi,xi,bini)/maxv2));
    end
    
  end
end
for xi = 0:ncells_x
    plot([(xi*psize+x_c+hsz-grid_hsz_x)*scale,(xi*psize+x_c+hsz-grid_hsz_x)*scale],...
      [ (y_c+hsz-grid_hsz_y)*scale,(y_c+ncells_y*psize+hsz-grid_hsz_y)*scale],...
      'Color',[0.3,0.3,0.3,0.2]);
end
for yi = 0:ncells_y
    plot([(x_c+hsz-grid_hsz_x)*scale,(x_c+ncells_x*psize+hsz-grid_hsz_x)*scale],...
      [ (yi*psize+y_c+hsz-grid_hsz_y)*scale,(yi*psize+y_c+hsz-grid_hsz_y)*scale],...
      'Color',[0.3,0.3,0.3,0.2]);
end
scatter(hsz*scale,hsz*scale,'marker','*');
scatter((x_c+hsz)*scale,(y_c+hsz)*scale,'marker','o');

set(handles.figure1, 'pointer', 'arrow')


% --- Outputs from this function are returned to the command line.
function varargout = spaceTimeOptions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure

varargout{1} = ~handles.is_canceled;
varargout{2} = handles.params;
figure1_CloseRequestFcn(hObject,eventdata,handles);


function edit_patch_Callback(hObject, eventdata, handles)
% hObject    handle to edit_patch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_patch as text
%        str2double(get(hObject,'String')) returns contents of edit_patch as a double

patch_sz = str2double(get(hObject,'String'));
prev_patch_sz = handles.params.psize;
if isempty(patch_sz) || isnan(patch_sz) || round(patch_sz)~=patch_sz || patch_sz < 1
    warndlg('Cell size must be a positive integer');
    set(handles.edit_patch,'String',sprintf('%d',prev_patch_sz));
    return;
end
handles.params.psize = patch_sz;
handles = updatePlots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_patch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_patch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_horz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_horz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_horz as text
%        str2double(get(hObject,'String')) returns contents of edit_horz as a double
ncells = str2double(get(hObject,'String'));
prev_cells = handles.params.npatches_x;
if isempty(ncells) || isnan(ncells) || round(ncells)~=ncells || ncells < 1
    warndlg('Cell size must be a positive integer');
    set(handles.edit_patch,'String',sprintf('%d',prev_cells));
    return;
end
handles.params.npatches_x = ncells;
handles = updatePlots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_horz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_horz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frame as text
%        str2double(get(hObject,'String')) returns contents of edit_frame as a double

frame_num = str2double(get(hObject,'String'));
if isempty(frame_num) || isnan(frame_num) || round(frame_num)~=frame_num || frame_num < 1
    warndlg('Frame number must be a positive integer');
    set(handles.edit_frame,'String',sprintf('%d',handles.frame_num));
    return;
end
min_frame = handles.trx(handles.cur_fly).firstframe;
max_frame = handles.trx(handles.cur_fly).endframe;
if frame_num > max_frame || frame_num < min_frame
    warndlg(sprintf('Frame number must be between %d and %d for the current animal',min_frame,max_frame));
    set(handles.edit_frame,'String',sprintf('%d',handles.frame_num));
    return;
end    
handles.frame_num = frame_num;
handles = updatePlots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_trx_Callback(hObject, eventdata, handles)
% hObject    handle to edit_trx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_trx as text
%        str2double(get(hObject,'String')) returns contents of edit_trx as a double

fly_num = str2double(get(hObject,'String'));
if isempty(fly_num) || isnan(fly_num) || round(fly_num)~=fly_num || fly_num < 1
    warndlg('Animal number must be a positive integer');
    set(handles.edit_trx,'String',sprintf('%d',handles.cur_fly));
    return;
end
if fly_num > handles.n_trx
    warndlg(sprintf('Animal number must be less than %d',handles.n_trx));
    set(handles.edit_trx,'String',sprintf('%d',handles.cur_fly));
    return;
end    
handles.cur_fly = fly_num;
if handles.frame_num < handles.trx(fly_num).firstframe
    handles.frame_num = handles.trx(fly_num).firstframe;
end
if handles.frame_num > (handles.trx(fly_num).endframe-1)
    handles.frame_num = handles.trx(fly_num).endframe-1;
end
updateGUI(handles);
handles = updatePlots(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_trx_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_trx (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.is_canceled=true;
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.is_canceled=false;
guidata(hObject,handles);
uiresume(handles.figure1);


% --- Executes on button press in toggle_align.
function toggle_align_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_align
if get(hObject,'Value')
    handles.params.is_stationary= true;
else
    handles.params.is_stationary= false;
end
handles = updatePlots(handles);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function toggle_align_CreateFcn(hObject, eventdata, handles)
% hObject    handle to toggle_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
try
    close(handles.fid);
catch ME
   
end
delete(hObject);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit_vert_Callback(hObject, eventdata, handles)
% hObject    handle to edit_vert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_vert as text
%        str2double(get(hObject,'String')) returns contents of edit_vert as a double
ncells = str2double(get(hObject,'String'));
prev_cells = handles.params.npatches_y;
if isempty(ncells) || isnan(ncells) || round(ncells)~=ncells || ncells < 1
    warndlg('Cell size must be a positive integer');
    set(handles.edit_patch,'String',sprintf('%d',prev_cells));
    return;
end
handles.params.npatches_y = ncells;
handles = updatePlots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_vert_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_vert (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in menu_methods.
function menu_methods_Callback(hObject, eventdata, handles)
% hObject    handle to menu_methods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns menu_methods contents as cell array
%        contents{get(hObject,'Value')} returns selected item from menu_methods
method_ndx = get(hObject,'Value');
handles.params.cur_method = handles.params.methods{method_ndx};
handles = updatePlots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function menu_methods_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menu_methods (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_x_center_Callback(hObject, eventdata, handles)
% hObject    handle to edit_x_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_x_center as text
%        str2double(get(hObject,'String')) returns contents of edit_x_center as a double
ncells = str2double(get(hObject,'String'));
prev_cells = handles.params.x_center;
if isempty(ncells) || isnan(ncells) || round(ncells)~=ncells 
    warndlg('Cell size must be a integer');
    set(handles.edit_x_center,'String',sprintf('%d',prev_cells));
    return;
end
handles.params.x_center = ncells;
handles = updatePlots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_x_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_x_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_y_center_Callback(hObject, eventdata, handles)
% hObject    handle to edit_y_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_y_center as text
%        str2double(get(hObject,'String')) returns contents of edit_y_center as a double
ncells = str2double(get(hObject,'String'));
prev_cells = handles.params.y_center;
if isempty(ncells) || isnan(ncells) || round(ncells)~=ncells 
    warndlg('Cell size must be a integer');
    set(handles.edit_y_center,'String',sprintf('%d',prev_cells));
    return;
end
handles.params.y_center = ncells;
handles = updatePlots(handles);
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function edit_y_center_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_y_center (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_imsz_Callback(hObject, eventdata, handles)
% hObject    handle to edit_imsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_imsz as text
%        str2double(get(hObject,'String')) returns contents of edit_imsz as a double
ncells = str2double(get(hObject,'String'));
prev_cells = handles.params.imsz_display;
if isempty(ncells) || isnan(ncells) || round(ncells)~=ncells || ncells < 1
    warndlg('Cell size must be a positive integer');
    set(handles.edit_imsz,'String',sprintf('%d',prev_cells));
    return;
end
params = handles.params;
grid_sz_x = params.npatches_x*params.psize;
grid_sz_y = params.npatches_y*params.psize;
x_c = params.x_center; y_c = params.y_center;
max_sz_x = max(x_c+grid_sz_x/2,-x_c+grid_sz_x/2);
max_sz_y = max(y_c+grid_sz_y/2,-y_c+grid_sz_y/2);
max_sz = max(max_sz_x,max_sz_y)*2;
prev_sz = params.imsz_display;
if ncells < max_sz
  warndlg('Input size cannot include the whole grid. Select a size large enough to include the whole grid');
  set(handles.edit_imsz,'String',sprintf('%d',prev_sz));
  return;
end

handles.params.imsz_display = ncells;
handles = updatePlots(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_imsz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_imsz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
