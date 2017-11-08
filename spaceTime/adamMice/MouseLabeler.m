function varargout = MouseLabeler(varargin)
% MOUSELABELER MATLAB code for MouseLabeler.fig
%      MOUSELABELER, by itself, creates a new MOUSELABELER or raises the existing
%      singleton*.
%
%      H = MOUSELABELER returns the handle to a new MOUSELABELER or the handle to
%      the existing singleton*.
%
%      MOUSELABELER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MOUSELABELER.M with the given input arguments.
%
%      MOUSELABELER('Property','Value',...) creates a new MOUSELABELER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MouseLabeler_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MouseLabeler_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MouseLabeler

% Last Modified by GUIDE v2.5 25-Aug-2014 09:48:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MouseLabeler_OpeningFcn, ...
                   'gui_OutputFcn',  @MouseLabeler_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}) && exist(varargin{1}),
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before MouseLabeler is made visible.
function MouseLabeler_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MouseLabeler (see VARARGIN)

% Choose default command line output for MouseLabeler
global MOUSELABELERLASTPATH;
handles.output = hObject;

% parse inputs
[handles.expdirs,handles.moviefilestr,handles.frontside,handles.savefile] = myparse(varargin,'expdirs',{},...
  'moviefilestr','movie_comb.avi','frontside',true,'loadfile','');

didload = false;
if ~isempty(handles.savefile),

  savedata = load(handles.savefile);
  fns = fieldnames(savedata);
  for i = 1:numel(fns),
    fn = fns{i};
    handles.(fn) = savedata.(fn);
  end
  didload = true;
  
else
  
  % select video to label
  if isempty(handles.expdirs)
    handles.expdirs = SelectExpDirs(handles.expdirs);
    if ~iscell(handles.expdirs) || isempty(handles.expdirs),
      delete(hObject);
      return;
    end
  end
  
end

MOUSELABELERLASTPATH = fileparts(handles.expdirs{end});

handles.keypressmode = '';
handles.motionobj = [];
if ~didload,
  handles.expi = 1;
end

if handles.frontside,
  handles.npoints = 2;
  handles.template = nan(2,2);
  handles.templatecolors = [.7,0,.7;0,.7,.7];
  set(handles.togglebutton_labelfront,'Visible','on');
  set([handles.text_nframeslabeledfront,handles.text_firstframefront,handles.text_lastframefront],'Visible','on');
else
  handles.npoints = 1;
  handles.template = nan(1,2);
  handles.templatecolors = [.7,0,.7];
  set(handles.togglebutton_labelfront,'Visible','off');
  set([handles.text_nframeslabeledfront,handles.text_firstframefront,handles.text_lastframefront],'Visible','off');
end
handles.pointselected = false(1,handles.npoints);

handles = InitializeGUI(handles);

% initialize template -- two points

if ~didload,
  handles.labeledpos = nan([handles.npoints,2,handles.nframes]);
  handles.labeledpos_perexp = cell(1,numel(handles.expdirs));
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MouseLabeler wait for user response (see UIRESUME)
% uiwait(handles.figure);

function handles = SetTemplate(handles)
  
uiwait(msgbox('Click to create template points. First, click to create each point. Then you can drag points around. Hit escape when done.'));

handles.hpoly = [];
handles.keypressmode = 'settemplate';
handles.template = nan(0,2);

axes(handles.axes_curr);

while true,

  keydown = waitforbuttonpress;
  if get(0,'CurrentFigure') ~= handles.figure,
    continue;
  end
  if keydown == 0 && strcmpi(get(handles.figure,'SelectionType'),'normal'),
    tmp = get(handles.axes_curr,'CurrentPoint');
    x = tmp(1,1);
    y = tmp(1,2);
    handles.hpoly(end+1) = plot(handles.axes_curr,x,y,'w+','MarkerSize',12);%,...
      %'KeyPressFcn',handles.keypressfcn);
    handles.template(end+1,:) = [x,y];
  elseif keydown == 1 && double(get(handles.figure,'CurrentCharacter')) == 27,
    break;
  end
  
end

handles.npoints = numel(handles.hpoly);
handles.templatecolors = jet(handles.npoints)*.5+.5;
for i = 1:handles.npoints,

  set(handles.hpoly(i),'Color',handles.templatecolors(i,:),...
    'ButtonDownFcn',@(hObject,eventdata) PointButtonDownCallback(hObject,eventdata,handles.figure,i));
  %addNewPositionCallback(handles.hpoly(i),@(pos) UpdateLabels(pos,handles.figure,i));
  
end

function PointButtonDownCallback(hObject,eventdata,hfig,i)

handles = guidata(hfig);
handles = StopLabeling(handles);

handles.motionobj = i;
handles.didmovepoint = false;
guidata(hObject,handles);

function UpdateLabels(pos,hfig,i,handles)

if nargin < 4,
  handles = guidata(hfig);
end
handles.labeledpos(i,:,handles.f) = pos;
guidata(hfig,handles);

function handles = InitializeVideo(handles)

% open video
handles.moviefile = fullfile(handles.expdirs{handles.expi},handles.moviefilestr);
[handles.readframe,handles.nframes,handles.fid,handles.headerinfo] = ...
  get_readframe_fcn(handles.moviefile);

[handles.imcurr,~] = handles.readframe(handles.f);
handles.f_im = handles.f;
handles.imprev = handles.imcurr;
handles.fprev_im = handles.f;

handles.minv = max(handles.minv,0);
if isfield(handles.headerinfo,'bitdepth'),
  handles.maxv = min(handles.maxv,2^handles.headerinfo.bitdepth-1);
elseif isa(handles.imcurr,'uint16'),
  handles.maxv = min(2^16 - 1,handles.maxv);
elseif isa(handles.imcurr,'uint8'),
  handles.maxv = min(handles.maxv,2^8 - 1);
else
  handles.maxv = min(handles.maxv,2^(ceil(log2(max(handles.imcurr(:)))/8)*8));
end

handles.imsize = size(handles.imcurr);

set(handles.axes_curr,'CLim',[handles.minv,handles.maxv],...
  'XLim',[.5,size(handles.imcurr,2)+.5],...
  'YLim',[.5,size(handles.imcurr,1)+.5]);
set(handles.axes_prev,'CLim',[handles.minv,handles.maxv],...
  'XLim',[.5,size(handles.imcurr,2)+.5],...
  'YLim',[.5,size(handles.imcurr,1)+.5]);

sliderstep = [1/(handles.nframes-1),min(1,10/(handles.nframes-1))];
set(handles.slider_frame,'Value',0,'SliderStep',sliderstep);

if isfield(handles,'labeledpos_perexp') && numel(handles.labeledpos_perexp) >= handles.expi && ...
    ~isempty(handles.labeledpos_perexp{handles.expi}),
  handles.labeledpos = handles.labeledpos_perexp{handles.expi};
elseif isfield(handles,'npoints'),
  handles.labeledpos = nan([handles.npoints,2,handles.nframes]);
end

if handles.frontside,
  handles.template(:,2) = (handles.imsize(1)+1)/2;
  handles.template(1,1) = (handles.imsize(1)+1)/4;
  handles.template(2,1) = 3*(handles.imsize(1)+1)/4;
else
  handles.template(:,2) = (handles.imsize(1)+1)/2;
  handles.template(1,1) = (handles.imsize(1)+1)/2;
end
handles.pointselected(:) = false;

maxchars = 40;
if numel(handles.expdirs{handles.expi}) <= maxchars+3,
  s = [repmat(' ',maxchars+3-numel(handles.expdirs{handles.expi})),handles.expdirs{handles.expi}];
else
  s = ['...',handles.expdirs{handles.expi}(end-maxchars+1:end)];
end
set(handles.text_currexp,'String',s);

UpdateNLabeled(handles);

function UpdateNLabeled(handles)

if ~isfield(handles,'labeledpos') || isempty(handles.labeledpos),
  nlabels = 0;
  nlabelsfront = 0;
else
  nlabels = nnz(~isnan(handles.labeledpos(1,1,:)));
  nlabelsfront = nnz(~isnan(handles.labeledpos(2,1,:)));
end
if nlabels == 0,
  f0 = '--';
  f1 = '--';
else
  f0 = num2str(find(~isnan(handles.labeledpos(1,1,:)),1));
  f1 = num2str(find(~isnan(handles.labeledpos(1,1,:)),1,'last'));
end

set(handles.text_nframeslabeled,'String',num2str(nlabels));
set(handles.text_firstframe,'String',f0);
set(handles.text_lastframe,'String',f1);

if nlabelsfront == 0,
  f0 = '--';
  f1 = '--';
else
  f0 = num2str(find(~isnan(handles.labeledpos(2,1,:)),1));
  f1 = num2str(find(~isnan(handles.labeledpos(2,1,:)),1,'last'));
end

set(handles.text_nframeslabeledfront,'String',num2str(nlabelsfront));
set(handles.text_firstframefront,'String',f0);
set(handles.text_lastframefront,'String',f1);


function handles = InitializeGUI(handles)

handles.islabeling = 0;

set(handles.figure,'ToolBar','figure');

handles.image_curr = imagesc(0,'Parent',handles.axes_curr,...
  'ButtonDownFcn',get(handles.axes_curr,'ButtonDownFcn'));
axis(handles.axes_curr,'image','off');
hold(handles.axes_curr,'on');
colormap(handles.figure,gray(256));

handles.image_prev = imagesc(0,'Parent',handles.axes_prev);
axis(handles.axes_prev,'image','off');
hold(handles.axes_prev,'on');
handles.posprev = [];

linkaxes([handles.axes_prev,handles.axes_curr]);

fcn = get(handles.slider_frame,'Callback');
handles.hslider_listener = handle.listener(handles.slider_frame,...
  'ActionEvent',fcn);
set(handles.slider_frame,'Callback','');

handles.hpoly = [];

handles.f = 1;
handles.maxv = inf;
handles.minv = 0;

handles.labelbuttoncolor = get(handles.togglebutton_label,'BackgroundColor');

handles.hpoly = nan(1,handles.npoints);
for i = 1:handles.npoints,
  handles.hpoly(i) = plot(handles.axes_curr,nan,nan,'+','Color',handles.templatecolors(i,:),'MarkerSize',12,...
    'ButtonDownFcn',@(hObject,eventdata) PointButtonDownCallback(hObject,eventdata,handles.figure,i),'LineWidth',2);
end
handles.posprev = nan(1,handles.npoints);
for i = 1:handles.npoints,
  handles.posprev(i) = plot(handles.axes_prev,nan,nan,'+','Color',handles.templatecolors(i,:),'MarkerSize',8);
end

if handles.frontside,
  set([handles.text_side,handles.text_nframeslabeled,handles.text_firstframe,handles.text_lastframe],...
    'ForegroundColor',handles.templatecolors(1,:));
  set([handles.text_front,handles.text_nframeslabeledfront,handles.text_firstframefront,handles.text_lastframefront],...
    'ForegroundColor',handles.templatecolors(2,:));
else
  set([handles.text_front,handles.text_nframeslabeledfront,handles.text_firstframefront,handles.text_lastframefront],...
    'Visible','off');
  set(handles.togglebutton_labelfront,'Visible','off');
  set([handles.menu_go_nextmismatch,handles.menu_go_prevmismatch],'Visible','off');
end
  

handles = InitializeVideo(handles);

handles = UpdateFrame(handles);

hchil = findall(handles.figure,'-property','KeyPressFcn');
handles.keypressfcn = get(handles.figure,'KeyPressFcn');
set(hchil,'KeyPressFcn',handles.keypressfcn);

function handles = UpdateFrame(handles)

if handles.f > 1,
  if handles.fprev_im ~= handles.f-1,
    [handles.imprev,~] = handles.readframe(handles.f-1);
    handles.fprev_im = handles.f-1;
  end
end

if handles.f_im ~= handles.f,  
  [handles.imcurr,~] = handles.readframe(handles.f);
  handles.f_im = handles.f;
end

if handles.f > 1 && isfield(handles,'imprev'),
  set(handles.image_prev,'CData',handles.imprev);
else
  set(handles.image_prev,'CData',0);
end
if isfield(handles,'imcurr'),
  set(handles.image_curr,'CData',handles.imcurr);
end
set(handles.slider_frame,'Value',(handles.f-1)/(handles.nframes-1));
set(handles.edit_frame,'String',num2str(handles.f));

for i = 1:handles.npoints,
  set(handles.hpoly(i),'XData',handles.labeledpos(i,1,handles.f),...
    'YData',handles.labeledpos(i,2,handles.f));
end

if handles.f > 1,
  for i = 1:handles.npoints,
    set(handles.posprev(i),'XData',handles.labeledpos(i,1,handles.f-1),...
      'YData',handles.labeledpos(i,2,handles.f-1));
  end
else
  set(handles.posprev,'XData',nan,'YData',nan);
end
  


% --- Outputs from this function are returned to the command line.
function varargout = MouseLabeler_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider_frame_Callback(hObject, eventdata, handles)
% hObject    handle to slider_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
handles.f = round(1 + v * (handles.nframes - 1));
handles = UpdateFrame(handles);
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function slider_frame_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function menu_file_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_file_save_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = StopLabeling(handles);

handles = SaveState(handles);
guidata(hObject,handles);

function handles = SaveState(handles)

handles.labeledpos_perexp{handles.expi} = handles.labeledpos;

global MOUSELABELERSAVEFILE;

savedata = struct;
savedata.labeledpos_perexp = handles.labeledpos_perexp;
savedata.expdirs = handles.expdirs;
savedata.frontside = handles.frontside;
savedata.npoints = handles.npoints;
savedata.expi = handles.expi;
savedata.f = handles.f;
savedata.minv = handles.minv;
savedata.maxv = handles.maxv;

if ~isfield(handles,'savefile'),
  handles.savefile = '';
end

[f,p] = uiputfile('*.mat','Save labels to file',handles.savefile);
if ~ischar(f),
  return;
end
handles.savefile = fullfile(p,f);
MOUSELABELERSAVEFILE = handles.savefile;

save(handles.savefile,'-struct','savedata');

% --------------------------------------------------------------------
function menu_file_quit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_quit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CloseGUI(handles);


% --------------------------------------------------------------------
function menu_setup_Callback(hObject, eventdata, handles)
% hObject    handle to menu_setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_setup_settemplate_Callback(hObject, eventdata, handles)
% hObject    handle to menu_setup_settemplate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.template),
  
  res = questdlg('Changing template will result in all labels being cleared. Save before doing this?');
  if strcmpi(res,'Cancel'),
    return;
  elseif strcmpi(res,'Yes'),
    handles = SaveState(handles);
  end

  % delete the current template
  for i = 1:handles.npoints,
    try %#ok<TRYNC>
      delete(handles.hpoly(i));
    end
  end
  
  handles.template = [];
  handles.npoints = 0;
  
end

handles = SetTemplate(handles);

handles.labeledpos = nan([handles.npoints,2,handles.nframes]);
handles.labeledpos(:,:,handles.f) = handles.template;
handles.pointselected = false(1,handles.npoints);

delete(handles.posprev(ishandle(handles.posprev)));
handles.posprev = nan(1,handles.npoints);
for i = 1:handles.npoints,
  handles.posprev(i) = plot(handles.axes_prev,nan,nan,'+','Color',handles.templatecolors(i,:),'MarkerSize',8);%,...
    %'KeyPressFcn',handles.keypressfcn);
end

guidata(hObject,handles);

function CloseImContrast(hObject)

handles = guidata(hObject);
clim = get(handles.axes_curr,'CLim');
handles.minv = clim(1);
handles.maxv = clim(2);
set(handles.axes_prev,'CLim',[handles.minv,handles.maxv]);
set(handles.axes_curr,'CLim',[handles.minv,handles.maxv]);
delete(handles.adjustbrightness_listener);
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_setup_adjustbrightness_Callback(hObject, eventdata, handles)
% hObject    handle to menu_setup_adjustbrightness (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%set(handles.axes_curr,'CLim',[min(handles.imcurr(:)),max(handles.imcurr(:))]);
hcontrast = imcontrast_kb(handles.axes_curr);
handles.adjustbrightness_listener = addlistener(hcontrast,'ObjectBeingDestroyed',@(x,y) CloseImContrast(hObject));
guidata(hObject,handles);


% --- Executes on button press in togglebutton_label.
function togglebutton_label_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_label (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if get(hObject,'Value'),
  
  if handles.islabeling == 2,
    set(handles.togglebutton_labelfront,'Value',0,'BackgroundColor',handles.labelbuttoncolor,'String','Label front');
  end
  
  handles.islabeling = 1;
else
  handles.islabeling = 0;
end
  
if handles.islabeling == 1,
  set(hObject,'BackgroundColor',handles.templatecolors(1,:),'String','Labeling side');
else
  set(hObject,'BackgroundColor',handles.labelbuttoncolor,'String','Label side');
end
guidata(hObject,handles);
% 
% while true,
%   try
%     handles = guidata(hObject);
%   catch
%     return;
%   end
%   if handles.islabeling ~= 1,
%     break;
%   end
% 
%   LabelFrame(handles,hObject);
%   
%   if handles.f == handles.nframes,
%     handles.islabeling = 0;
%     break;
%   end
%   
% end
% 
% set(hObject,'Value',0,'BackgroundColor',[0,.6,0],'String','Label side');
% guidata(hObject,handles);

% function LabelFrame(handles,hObject)
% 
% if handles.islabeling == 0,
%   return;
% end
% 
% % handles.labeledpos(handles.islabeling,:,handles.f) = nan;
% % set(handles.hpoly(handles.islabeling),'XData',nan,'YData',nan);
% 
% i = handles.islabeling;
% 
% while true,
%       
%   try
%     keydown = waitforbuttonpress;
%     if ~ishandle(handles.figure),
%       return;
%     end
%     objclicked = gco;
%     eventdata = struct;
%     eventdata.Character = get(handles.figure,'CurrentCharacter');
%     eventdata.Modifier = get(handles.figure,'CurrentModifier');
%     eventdata.Key = get(handles.figure,'CurrentKey');
%     tmp = get(handles.axes_curr,'CurrentPoint');
%     x = tmp(1,1);
%     y = tmp(1,2);
%     handles = guidata(hObject);
%   catch ME
%     if strcmpi(ME.identifier,'MATLAB:UI:CancelWaitForKeyOrButtonPress'),
%       return;
%     else
%       warning(getReport(ME));
%       return;
%     end
%   end
%   if get(0,'CurrentFigure') ~= handles.figure,
%     continue;
%   end
%   if keydown == 0 && strcmpi(get(handles.figure,'SelectionType'),'normal'),
%     if objclicked ~= handles.image_curr,
%       j = find(objclicked == handles.hpoly);
%       if ~isempty(j),
%         PointButtonDownCallback(objclicked,[],handles.figure,j);
%         return;
%       end
%     else
%       xlim = get(handles.axes_curr,'XLim');
%       ylim = get(handles.axes_curr,'YLim');
%       if x < xlim(1) || x > xlim(2) || y < ylim(1) || y > ylim(2),
%         continue;
%       end
%       set(handles.hpoly(i),'XData',x,'YData',y);
%       handles.labeledpos(i,:,handles.f) = [x,y];
%       UpdateNLabeled(handles);
%       break;
%     end
%   elseif keydown == 1,
%     figure_KeyPressFcn(hObject, eventdata, handles);
%   end
%   
% end
% 
% if handles.islabeling==0 || handles.f == handles.nframes,
%   
%   if handles.islabeling == 1,
%     handles.islabeling = 0;
%     set(handles.togglebutton_label,'Value',0);
%     togglebutton_label_Callback(handles.togglebutton_label,[],handles);
%   elseif handles.islabeling == 2,
%     set(handles.togglebutton_labelfront,'Value',0);
%     togglebutton_labelfront_Callback(handles.togglebutton_labelfront,[],handles);
%   end
%   guidata(hObject,handles);
%   return;
%   
% end
% 
% handles.f = min(handles.f + 10,handles.nframes);
% handles = UpdateFrame(handles);
% guidata(hObject,handles);

function handles = StopLabeling(handles)

if ~isfield(handles,'islabeling') || handles.islabeling == 0,
  return;
end

handles.labeledpos_perexp{handles.expi} = handles.labeledpos;

if handles.islabeling == 1,
  set(handles.togglebutton_label,'Value',0);
  handles.islabeling = 0;
  togglebutton_label_Callback(handles.togglebutton_label, [], handles);
elseif handles.islabeling == 2,
  set(handles.togglebutton_labelfront,'Value',0);
  handles.islabeling = 0;
  togglebutton_labelfront_Callback(handles.togglebutton_labelfront, [], handles);
end
  
if nargout > 0,
  handles = guidata(handles.togglebutton_label);
end

% Hint: get(hObject,'Value') returns toggle state of togglebutton_label


% --- Executes when user attempts to close figure.
function figure_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
CloseGUI(handles);


function CloseGUI(handles)

handles = StopLabeling(handles);

res = questdlg('Save before closing?');
if strcmpi(res,'Cancel'),
  return;
elseif strcmpi(res,'Yes'),
  SaveState(handles);
end

delete(handles.figure);


% --------------------------------------------------------------------
function menu_file_load_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = StopLabeling(handles);

res = questdlg('Save before closing current video?');
if strcmpi(res,'Cancel'),
  return;
elseif strcmpi(res,'Yes'),
  handles = SaveState(handles);
end

global MOUSELABELERSAVEFILE;
if isempty(MOUSELABELERSAVEFILE),
  defaultfile = '';
else
  defaultfile = MOUSELABELERSAVEFILE;
end
[f,p] = uigetfile('*.mat','Load file...',defaultfile);
if ~ischar(f),
  return;
end
handles.savefile = fullfile(p,f);
if ~exist(handles.savefile,'file'),
  warndlg(sprintf('File %s does not exist',handles.savefile),'File does not exist','modal');
  return;
end

MOUSELABELERSAVEFILE = handles.savefile;

savedata = load(handles.savefile);
fns = fieldnames(savedata);
olddata = struct;
for i = 1:numel(fns),
  fn = fns{i};
  olddata.(fn) = handles.(fn);
  handles.(fn) = savedata.(fn);
end

if isfield(handles,'fid') && isnumeric(handles.fid) && handles.fid > 0,
  fclose(handles.fid);
end
if handles.frontside,
  handles.npoints = 2;
  handles.template = nan(2,2);
  handles.templatecolors = [.7,0,.7;0,.7,.7];
  set(handles.togglebutton_labelfront,'Visible','on');
  set([handles.text_nframeslabeledfront,handles.text_firstframefront,handles.text_lastframefront],'Visible','on');
else
  handles.npoints = 1;
  handles.template = nan(1,2);
  handles.templatecolors = [.7,0,.7];
  set(handles.togglebutton_labelfront,'Visible','off');
  set([handles.text_nframeslabeledfront,handles.text_firstframefront,handles.text_lastframefront],'Visible','off');
end
handles.pointselected = false(1,handles.npoints);
handles = InitializeVideo(handles);
handles.f_im = nan;
handles.fprev_im = nan;

handles = UpdateFrame(handles);

guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_file_switchexp_Callback(hObject, eventdata, handles, expi)
% hObject    handle to menu_file_switchexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = StopLabeling(handles);

handles.labeledpos_perexp{handles.expi} = handles.labeledpos;

if nargin < 4,
  expi = SelectExp(handles);  
end
if isempty(expi),
  return;
end

if handles.fid > 1,
  try
    fclose(handles.fid);
  end
end
handles.expi = expi;

handles.labeledpos = [];
handles.f = 1;
handles.minv = 0;
handles.maxv = inf;

handles = InitializeVideo(handles);
handles.f_im = nan;
handles.fprev_im = nan;

handles = UpdateFrame(handles);

guidata(hObject,handles);

function expi = SelectExp(handles)

maxchars = 60;
ss = cell(1,numel(handles.expdirs));
for i = 1:numel(handles.expdirs),
  if numel(handles.expdirs{i}) <= maxchars+3,
    ss{i} = [repmat(' ',maxchars+3-numel(handles.expdirs{i})),handles.expdirs{i}];
  else
    ss{i} = ['...',handles.expdirs{i}(end-maxchars+1:end)];
  end
  if isempty(handles.labeledpos_perexp{i}),
    nlabels = 0;
  else
    nlabels = nnz(all(~isnan(handles.labeledpos_perexp{i}(:,1,:)),1));
  end
  ss{i} = sprintf('%s - %d labels',ss{i},nlabels);
end
[expi,v] = listdlg('ListString',ss,'SelectionMode','single','InitialValue',handles.expi,...
  'Name','Select experiment','PromptString','Switch to experiment:',...
  'ListSize',[600,300]);
if v == 0,
  expi = [];
  return;
end

function expdirs = SelectExpDirs(expdirs)

global MOUSELABELERLASTPATH;

if isempty(MOUSELABELERLASTPATH),
  MOUSELABELERLASTPATH = '';
end

expdirs = uipickfiles('FilterSpec',MOUSELABELERLASTPATH,...
  'DirsOnly',true,'Append',expdirs);



function moviefile = SelectVideo()

global MOUSELABELERLASTPATH;

if isempty(MOUSELABELERLASTPATH),
  MOUSELABELERLASTPATH = '';
end

[f,p] = uigetfile('*.*','Select video to label',MOUSELABELERLASTPATH);
if ~ischar(f),
  return;
end
MOUSELABELERLASTPATH = p;
moviefile = fullfile(p,f);



function edit_frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frame as text
%        str2double(get(hObject,'String')) returns contents of edit_frame as a double
f = str2double(get(hObject,'String'));
if isnan(f),
  set(hObject,'String',num2str(handles.f));
  return;
end
f = min(max(1,round(f)),handles.nframes);
set(hObject,'String',num2str(f));
if f ~= handles.f,
  handles.f = f;
  handles = UpdateFrame(handles);
end
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

function handles = GoToNextLabeledFrame(handles)
        
df = find(any(~isnan(handles.labeledpos(:,1,handles.f+1:end)),1),1);
if isempty(df),
  
  if handles.expi < numel(handles.expdirs),
    res = questdlg('Go to next video?');
    if strcmpi(res,'Yes'),
      menu_file_switchexp_Callback(handles.figure, [], handles, handles.expi+1)
    end
  end
  return;
  
end

f = min(handles.f+df,handles.nframes);
handles = GoToFrame(handles,f);

function handles = GoToPreviousLabeledFrame(handles)

f = find(any(~isnan(handles.labeledpos(:,1,1:handles.f-1)),1),1,'last');
if isempty(f),
  
  if handles.expi > 1,
    res = questdlg('Go to previous video?');
    if strcmpi(res,'Yes'),
      menu_file_switchexp_Callback(handles.figure, [], handles, handles.expi-1)
    end
  end
  return;
else
  df = handles.f-f;
end
f = max(handles.f-df,1);
handles = GoToFrame(handles,f);
    
function handles = GoToNextMismatchFrame(handles)
        
if ~handles.frontside,
  return;
end

df = find(isnan(handles.labeledpos(1,1,handles.f+1:end))~=isnan(handles.labeledpos(2,1,handles.f+1:end)),1);
if isempty(df),
  
  if handles.expi < numel(handles.expdirs),
    res = questdlg('Go to next video?');
    if strcmpi(res,'Yes'),
      menu_file_switchexp_Callback(handles.figure, [], handles, handles.expi+1)
    end
  end
  return;
  
end

f = min(handles.f+df,handles.nframes);
handles = GoToFrame(handles,f);

function handles = GoToPreviousMismatchFrame(handles)
        
if ~handles.frontside,
  return;
end

f = find(isnan(handles.labeledpos(1,1,1:handles.f-1))~=isnan(handles.labeledpos(2,1,1:handles.f-1)),1,'last');
if isempty(f),
  
  if handles.expi > 1,
    res = questdlg('Go to previous video?');
    if strcmpi(res,'Yes'),
      menu_file_switchexp_Callback(handles.figure, [], handles, handles.expi-1)
    end
  end
  return;
else
  df = handles.f-f;
end
f = max(handles.f-df,1);
handles = GoToFrame(handles,f);

function handles = GoToFrame(handles,f)

if f ~= handles.f,
  handles.f = f;
  handles = UpdateFrame(handles);
  guidata(handles.figure,handles);
end


% --- Executes on key press with focus on figure and none of its controls.
function figure_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)  

%fprintf('Pressed %s, > %s<\n',eventdata.Key,sprintf('%s ',eventdata.Modifier{:}));
  
  
switch eventdata.Key,
  case 'rightarrow',
    if any(handles.pointselected),
      xlim = get(handles.axes_curr,'XLim');
      dx = diff(xlim);
      if ismember('shift',eventdata.Modifier),
        dx = dx / 50;
      else
        dx = dx / 500;
      end
      for i = find(handles.pointselected),
        x = get(handles.hpoly(i),'XData');
        handles.labeledpos(i,1,handles.f) = x + dx;
        set(handles.hpoly(i),'XData',handles.labeledpos(i,1,handles.f));
        guidata(hObject,handles);
      end
    else
      if ismember('control',eventdata.Modifier),
        GoToNextLabeledFrame(handles);
      elseif ismember('alt',eventdata.Modifier),
        GoToNextMismatchFrame(handles);
      else
        if ismember('shift',eventdata.Modifier),
          df = 10;
        else
          df = 1;
        end
        f = min(handles.f+df,handles.nframes);
        GoToFrame(handles,f);
      end
    end
  case 'equal',
    if ismember('control',eventdata.Modifier),
      GoToNextLabeledFrame(handles);
    elseif ismember('alt',eventdata.Modifier),
      GoToNextMismatchFrame(handles);
    else
      if ismember('shift',eventdata.Modifier),
        df = 10;
      else
        df = 1;
      end
      f = min(handles.f+df,handles.nframes);
      GoToFrame(handles,f);
    end
  case 'hyphen',
    if ismember('control',eventdata.Modifier),
      GoToPreviousLabeledFrame(handles);
    elseif ismember('alt',eventdata.Modifier),
      GoToPreviousMismatchFrame(handles);
    else
      if ismember('shift',eventdata.Modifier),
        df = 10;
      else
        df = 1;
      end
      f = max(handles.f-df,1);
      GoToFrame(handles,f);
    end
  case 'leftarrow',
    if any(handles.pointselected),
      xlim = get(handles.axes_curr,'XLim');
      dx = diff(xlim);
      if ismember('shift',eventdata.Modifier),
        dx = dx / 50;
      else
        dx = dx / 500;
      end
      for i = find(handles.pointselected),
        x = get(handles.hpoly(i),'XData');
        handles.labeledpos(i,1,handles.f) = x - dx;
        set(handles.hpoly(i),'XData',handles.labeledpos(i,1,handles.f));
        guidata(hObject,handles);
      end
    else
      if ismember('control',eventdata.Modifier),
        GoToPreviousLabeledFrame(handles);
      elseif ismember('alt',eventdata.Modifier),
        GoToPreviousMismatchFrame(handles);
      else
        if ismember('shift',eventdata.Modifier),
          df = 10;
        else
          df = 1;
        end
        f = max(handles.f-df,1);
        GoToFrame(handles,f);
      end
    end
  case 'uparrow',
    if any(handles.pointselected),
      ylim = get(handles.axes_curr,'YLim');
      dy = diff(ylim);
      if ismember('shift',eventdata.Modifier),
        dy = dy / 50;
      else
        dy = dy / 500;
      end
      for i = find(handles.pointselected),
        y = get(handles.hpoly(i),'YData');
        handles.labeledpos(i,2,handles.f) = y - dy;
        set(handles.hpoly(i),'YData',handles.labeledpos(i,2,handles.f));
        guidata(hObject,handles);
      end
    end
  case 'downarrow',
    if any(handles.pointselected),
      ylim = get(handles.axes_curr,'YLim');
      dy = diff(ylim);
      if ismember('shift',eventdata.Modifier),
        dy = dy / 50;
      else
        dy = dy / 500;
      end
      for i = find(handles.pointselected),
        y = get(handles.hpoly(i),'YData');
        handles.labeledpos(i,2,handles.f) = y + dy;
        set(handles.hpoly(i),'YData',handles.labeledpos(i,2,handles.f));
        guidata(hObject,handles);
      end
    end
  case 's'
    set(handles.togglebutton_label,'Value',~get(handles.togglebutton_label,'Value'));
    togglebutton_label_Callback(handles.togglebutton_label,[],handles);
  case 'f'
    if handles.frontside,
      set(handles.togglebutton_labelfront,'Value',~get(handles.togglebutton_labelfront,'Value'));
      togglebutton_labelfront_Callback(handles.togglebutton_labelfront,[],handles);
    end
  case 'escape',
    for i = find(handles.pointselected),
      handles.pointselected(i) = false;
      UpdatePointSelected(handles,i);
    end
    if handles.islabeling == 1,
      set(handles.togglebutton_label,'Value',0);
      togglebutton_label_Callback(handles.togglebutton_label,[],handles);
    elseif handles.islabeling == 2,
      set(handles.togglebutton_labelfront,'Value',0);
      togglebutton_labelfront_Callback(handles.togglebutton_labelfront,[],handles);
    end
      
    guidata(hObject,handles);
end


% --- Executes on mouse motion over figure - except title and menu.
function figure_WindowButtonMotionFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'motionobj') || isempty(handles.motionobj), return; end

if isnumeric(handles.motionobj),
  handles.didmovepoint = true;
  tmp = get(handles.axes_curr,'CurrentPoint');
  pos = tmp(1,1:2);
  set(handles.hpoly(handles.motionobj),'XData',pos(1),'YData',pos(2));
  UpdateLabels(pos,hObject,handles.motionobj,handles);
end


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure_WindowButtonUpFcn(hObject, eventdata, handles)
% hObject    handle to figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ~isempty(handles.motionobj) && isnumeric(handles.motionobj),
  if ~handles.didmovepoint,
    handles.pointselected(handles.motionobj) = ~handles.pointselected(handles.motionobj);
    UpdatePointSelected(handles,handles.motionobj);    
  end
  handles.motionobj = [];
  guidata(hObject,handles);
end

function UpdatePointSelected(handles,i)

if handles.pointselected(i),
  set(handles.hpoly(i),'LineWidth',3,'Marker','o');
else
  set(handles.hpoly(i),'LineWidth',2,'Marker','+');
end


% --------------------------------------------------------------------
function menu_help_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_help_keyboardshortcuts_Callback(hObject, eventdata, handles)
% hObject    handle to menu_help_keyboardshortcuts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = {};
s{end+1} = '* When point is selected, LEFT, RIGHT, UP, and DOWN move the';
s{end+1} = '  selected point a small amount.';
s{end+1} = '* When point is selected, SHIFT+LEFT, RIGHT, UP, and DOWN move';
s{end+1} = '  the selected point a large amount.';
s{end+1} = '* When no template point is selected, LEFT and RIGHT decrement';
s{end+1} = '  and increment the frame shown.';
s{end+1} = '* MINUS (-) and EQUAL (=) always decrement and increment the';
s{end+1} = '  frame shown.';
s{end+1} = '* When no template point is selected, SHIFT+LEFT and';
s{end+1} = '  SHIFT+RIGHT decrease and increase the frame shown by 10.';
s{end+1} = '* SHIFT+MINUS and SHIFT+EQUAL always decrease and increase the';
s{end+1} = '  frame shown by 10, resp.';
s{end+1} = '* CTRL+MINUS and CTRL+EQUAL go to the previous and next labeled';
s{end+1} = '  frame, resp.';
s{end+1} = '* ALT+MINUS and ALT+EQUAL go to the previous and next';
s{end+1} = '  "mismatched" frame, resp. A "mismatched frame is one for';
s{end+1} = '  which the front view has been labeled but not the side view, ';
s{end+1} = '  or vice-versa.';
s{end+1} = '* S toggles whether you are labeling the side view for the ';
s{end+1} = '  current frame.';
s{end+1} = '* F toggles whether you are labeling the front view for the ';
s{end+1} = '  current frame.';

msgbox(s,'Keyboard shortcuts','help','modal');


% --------------------------------------------------------------------
function menu_file_addexp_Callback(hObject, eventdata, handles)
% hObject    handle to menu_file_addexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = StopLabeling(handles);

% select video to label
oldexpdirs = handles.expdirs;
newexpdirs = SelectExpDirs(handles.expdirs);
if ~iscell(handles.expdirs),
  return;
end

% add the new exps
[ism1] = ismember(newexpdirs,oldexpdirs);
[ism2] = ismember(oldexpdirs,newexpdirs);
handles.expdirs = [handles.expdirs,newexpdirs(~ism1)];
handles.labeledpos_perexp = [handles.labeledpos_perexp,cell(1,nnz(~ism1))];

% remove old exps
oldexpdir = handles.expdirs{handles.expi};
handles.expdirs(~ism2) = [];
handles.labeledpos_perexp(~ism2) = [];
i = find(strcmp(oldexpdir,handles.expdirs));
if isempty(i),
  handles.expi = 1;
  
  if handles.fid > 1,
    try
      fclose(handles.fid);
    end
  end
  
  handles.labeledpos = [];
  handles.f = 1;
  handles.minv = 0;
  handles.maxv = inf;
  
  handles = InitializeVideo(handles);
  handles.f_im = nan;
  handles.fprev_im = nan;
  
  handles = UpdateFrame(handles);
  
else
  handles.expi = i;
end

guidata(hObject,handles);


% --- Executes on button press in togglebutton_labelfront.
function togglebutton_labelfront_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton_labelfront (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~handles.frontside,
  return;
end

% Hint: get(hObject,'Value') returns toggle state of togglebutton_labelfront
if get(hObject,'Value'),
  if handles.islabeling == 1,
    set(handles.togglebutton_label,'Value',0,'BackgroundColor',handles.labelbuttoncolor,'String','Label side');
  end
  handles.islabeling = 2;
else
  handles.islabeling = 0;
end
  
if handles.islabeling == 2,
  set(hObject,'BackgroundColor',handles.templatecolors(2,:),'String','Labeling front');
else
  set(hObject,'BackgroundColor',handles.labelbuttoncolor,'String','Label front');
end
guidata(hObject,handles);
% 
% while true,
%   try
%     handles = guidata(hObject);
%   catch
%     return;
%   end
%   if handles.islabeling ~= 2,
%     break;
%   end
% 
%   LabelFrame(handles,hObject);
%   
%   if handles.f == handles.nframes,
%     handles.islabeling = 0;
%     break;
%   end
%   
% end
% 
% set(hObject,'Value',0,'BackgroundColor',[0,.6,0],'String','Label front');
% guidata(hObject,handles);


% --- Executes on mouse press over axes background.
function axes_curr_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes_curr (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.islabeling ~= 0,
  tmp = get(handles.axes_curr,'CurrentPoint');
  x = tmp(1,1);
  y = tmp(1,2);
  i = handles.islabeling;
  set(handles.hpoly(i),'XData',x,'YData',y);
  handles.labeledpos(i,:,handles.f) = [x,y];
  UpdateNLabeled(handles);
  
  if handles.f < handles.nframes,
    handles.f = min(handles.f + 10,handles.nframes);
    handles = UpdateFrame(handles);
  end
  
end
guidata(hObject,handles);
  


% --------------------------------------------------------------------
function menu_edit_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_go_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_go_nextlabeled_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_nextlabeled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

GoToNextLabeledFrame(handles);

% --------------------------------------------------------------------
function menu_go_previouslabeled_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_previouslabeled (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

GoToPreviousLabeledFrame(handles);

% --------------------------------------------------------------------
function menu_go_nextmismatch_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_nextmismatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

GoToNextMismatchFrame(handles);


% --------------------------------------------------------------------
function menu_go_prevmismatch_Callback(hObject, eventdata, handles)
% hObject    handle to menu_go_prevmismatch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

GoToPreviousMismatchFrame(handles);


% --------------------------------------------------------------------
function menu_edit_clearvideo_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_clearvideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'labeledpos') || isempty(handles.labeledpos),
  return;
end

res = questdlg('Delete all labels from the current video?');
if strcmpi(res,'Yes'),
  handles.labeledpos(:) = nan;
  UpdateNLabeled(handles);
end


% --------------------------------------------------------------------
function menu_edit_clearall_Callback(hObject, eventdata, handles)
% hObject    handle to menu_edit_clearall (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

res = questdlg('Delete all labels from the all video?');
if strcmpi(res,'Yes'),
  if isfield(handles,'labeledpos'),
    handles.labeledpos(:) = nan;
  end
  if ~isfield(handles,'labeledpos_perexp'),
    for i = 1:numel(handles.labeledpos_perexp),
      handles.labeledpos_perexp{i}(:) = nan;
    end
  end
  UpdateNLabeled(handles);
end
