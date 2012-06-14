function varargout = PlayJAABAResults(varargin)
% PLAYJAABARESULTS MATLAB code for PlayJAABAResults.fig
%      PLAYJAABARESULTS, by itself, creates a new PLAYJAABARESULTS or raises the existing
%      singleton*.
%
%      H = PLAYJAABARESULTS returns the handle to a new PLAYJAABARESULTS or the handle to
%      the existing singleton*.
%
%      PLAYJAABARESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PLAYJAABARESULTS.M with the given input arguments.
%
%      PLAYJAABARESULTS('Property','Value',...) creates a new PLAYJAABARESULTS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PlayJAABAResults_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PlayJAABAResults_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PlayJAABAResults

% Last Modified by GUIDE v2.5 10-Jun-2012 21:47:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PlayJAABAResults_OpeningFcn, ...
                   'gui_OutputFcn',  @PlayJAABAResults_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1}) && exist(varargin{1}) %#ok<EXIST>
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PlayJAABAResults is made visible.
function PlayJAABAResults_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PlayJAABAResults (see VARARGIN)

% Choose default command line output for PlayJAABAResults
handles.output = hObject;

persistent expdir;
persistent classifierparamsfile;

[handles.expdir,...
  handles.classifierparamsfile] = ...
  myparse(varargin,'expdir','',...
  'classifierparamsfile','');

if isempty(handles.expdir),
  handles.expdir = uigetdir(expdir,'Experiment directory');
  if ~ischar(handles.expdir),
    delete(handles.figure_PJR);
    return;
  end
  expdir = handles.expdir;
end
if isempty(handles.classifierparamsfile),
  [filename,path] = uigetfile(classifierparamsfile,'Classifier parameters file');
  if ~ischar(filename),
    delete(handles.figure_PJR);
    return;
  end
  handles.classifierparamsfile = fullfile(path,filename);
  classifierparamsfile = handles.classifierparamsfile;
end

handles = InitializeState(handles);

handles = InitializePlots(handles);

handles = SetCurrentFly(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PlayJAABAResults wait for user response (see UIRESUME)
% uiwait(handles.figure_PJR);

function handles = InitializeState(handles)

score_norm_prctile = 80;
  
%% read config parameters
handles.classifierparams = ReadClassifierParamsFile(handles.classifierparamsfile);
handles.nbehaviors = numel(handles.classifierparams);
handles.behaviors = cell(1,handles.nbehaviors);
handles.scorefns = cell(1,handles.nbehaviors);
for i = 1:handles.nbehaviors,
  if iscell(handles.classifierparams(i).behaviors.names),
    handles.behaviors{i} = sprintf('%s_',handles.classifierparams(i).behaviors.names{:});
    handles.behaviors{i} = handles.behaviors{i}(1:end-1);
  else
    handles.behaviors{i} = handles.classifierparams(i).behaviors.names;
  end
  handles.behaviors{i} = regexprep(handles.behaviors{i},'(_|^)([a-z])','${upper($2)}');
  scorefn = handles.classifierparams(i).file.scorefilename;
  handles.scorefns{i} = regexprep(scorefn,'\.mat$','');
end

%% create trx structure

global STORED_TRX;

if isempty(STORED_TRX),
handles.trx = Trx('trxfilestr',handles.classifierparams(1).file.trxfilename,...
  'moviefilestr',handles.classifierparams(1).file.moviefilename,...
  'perframedir',handles.classifierparams(1).file.perframedir);
if isfield(handles.classifierparams,'perframe') && isfield(handles.classifierparams(1).perframe,'params'),
  handles.trx.SetPerFrameParams(handles.classifierparams(1).perframe.params);
end  
if isfield(handles.classifierparams,'perframe') && isfield(handles.classifierparams(1).perframe,'landmark_params'),
  handles.trx.SetLandmarkParams(handles.classifierparams(1).perframe.landmark_params);
end  

handles.trx.AddExpDir(handles.expdir,'openmovie',false);
STORED_TRX = handles.trx;
else
  handles.trx = STORED_TRX;
end

%% for reading video

moviename = fullfile(handles.expdir,handles.trx.moviefilestr);
[handles.readframe,handles.nframes,handles.fid,handles.headerinfo] = get_readframe_fcn(moviename);

%% current fly, frame, etc. 

handles.target = 1;
handles.t = handles.trx(handles.target).firstframe;
handles.nflies = handles.trx.nflies;

%% playback fps

handles.fps = 3;
handles.dozoom = true;

%% choose weights for each behavior

handles.behaviorweight = nan(handles.nbehaviors,1);
for behaviori = 1:handles.nbehaviors,
  tmp = handles.trx.(handles.scorefns{behaviori});
  if iscell(tmp),
    tmp = [tmp{:}];
  end
  handles.behaviorweight(behaviori) = nnz(tmp>0);
end
handles.behaviorweight = handles.behaviorweight / sum(handles.behaviorweight);
handles.behaviorweight = handles.behaviorweight(:);

%% choose colors for each behavior

tmp = jet(256)*.7;
handles.behavior_colors = tmp(round(linspace(1,256,handles.nbehaviors)),:);

%% initial zoom

a = handles.trx.a;
if iscell(a),
  a = [a{:}];
end
meana = nanmedian(a)*4;
awidthradius = 3;
aheightradius = 3;
aborder = 1;
aback = .25;
handles.pxwidthradius = awidthradius*meana;
handles.pxheightradius = aheightradius*meana;
handles.pxborder = aborder*meana;
handles.pxback = meana*aback;
handles.meana = meana;

%% get the score normalizations

handles.score_norms = nan(handles.nbehaviors,1);
for i = 1:handles.nbehaviors,
  tmp = handles.trx.(handles.scorefns{i});
  if iscell(tmp),
    tmp = [tmp{:}];
  end
  handles.score_norms(i) = prctile(abs(tmp),score_norm_prctile);
end

%% not playing

handles.isplaying = false;

function handles = InitializePlots(handles)

trx_linewidth = 1;
target_linewidth = 2;
trxcolor = 'k';
textcolor = 'k';
fontsize = 30;
bar_width = .8;
barbackcolor = [.5,.5,.5];

%% main axis
% frame
im = handles.readframe(handles.t);
ncolors = size(im,3);
if ncolors == 1,
  handles.him = imagesc(im2double(im),'Parent',handles.haxmain,[0,1]);
  colormap(handles.haxmain,'gray');
else
  handles.him = image(im2double(im),'Parent',handles.haxmain);
end
hold(handles.haxmain,'on');

% plot current fly
handles.htrx = plot(handles.haxmain,nan,nan,'.-','LineWidth',trx_linewidth,'color',trxcolor,'HitTest','off');
handles.htarget = plot(handles.haxmain,nan,nan,'-','LineWidth',target_linewidth,'color','k','HitTest','off');
set(handles.haxmain,'XTick',[],'YTick',[]);

% plot other flies
handles.hothers = nan(1,handles.nflies);
for i = 1:handles.nflies,
  handles.hothers(i) = plot(handles.haxmain,nan,nan,'mo','markerfacecolor','m');
  set(handles.hothers(i),'ButtonDownFcn',@(hObject,eventdata) PlayJAABAResults('fly_ButtonDownFcn',hObject,eventdata,guidata(hObject),i));
end

% zoom
set(handles.haxmain,'Units','pixels');
handles.mainpos = get(handles.haxmain,'Position');
set(handles.haxmain,'Units','normalized');
pxheightradius1 = min((handles.pxwidthradius)*handles.mainpos(4) / handles.mainpos(3),(handles.headerinfo.nr-1)/2);
pxwidthradius1 = min((handles.pxheightradius)*handles.mainpos(3) / handles.mainpos(4),(handles.headerinfo.nc-1)/2);
handles.pxheightradius = max(pxheightradius1,handles.pxheightradius);
handles.pxwidthradius = max(pxwidthradius1,handles.pxwidthradius);
axis(handles.haxmain,'image');
%axis(handles.haxmain,handles.ax);

% target info
handles.hinfo = text(0,0,'','Parent',handles.haxmain,'Color',textcolor,...
  'FontUnits','pixels','FontSize',fontsize,...
  'HorizontalAlignment','left',...
  'VerticalAlignment','bottom',...
  'HitTest','off');

% current behaviors
handles.htext = text(nan,nan,'','Parent',handles.haxmain,'Color',textcolor,...
  'FontUnits','pixels','FontSize',fontsize);

%% bar plots
handles.hbar = nan(1,handles.nbehaviors);
handles.hbarback = nan(1,handles.nbehaviors);
for i = 1:handles.nbehaviors,
  handles.hbarback(i) = patch(i+bar_width/2*[-1,1,1,-1,-1],zeros(1,5)+.1,barbackcolor,'EdgeColor','none',...
    'Parent',handles.haxbottom);
  if i == 1,
    hold(handles.haxbottom,'on');
  end
  handles.hbar(i) = patch(i+bar_width/2*[-1,1,1,-1,-1],zeros(1,5)+.1,handles.behavior_colors(i,:),'EdgeColor','none','Parent',handles.haxbottom);
end

% zero line
plot(handles.haxbottom,[0,handles.nbehaviors+1],[0,0],'--','color',[.5,.5,.5]);
set(handles.haxbottom,'XTick',1:handles.nbehaviors,'XTickLabel',handles.behaviors,'YTick',[],...
  'XLim',[0,handles.nbehaviors+1],'YLim',[-1,1],'XColor','w','Color','k','TickLength',[0,0]);

% slider
set(handles.slider_frame,'Min',1,'Max',handles.nframes,...
  'Value',1,...
  'SliderStep',[1/(handles.nframes-1),100/(handles.nframes-1)]);

RecursiveSetKeyPressFcn(handles.figure_PJR);


%% set current fly

function handles = SetCurrentFly(handles)

handles.x = handles.trx(handles.target).x;
handles.y = handles.trx(handles.target).y;
handles.a = handles.trx(handles.target).a;
handles.b = handles.trx(handles.target).b;
handles.theta = handles.trx(handles.target).theta;
fil = [0.000263865082737   0.106450771973592   0.786570725887342   0.106450771973592   0.000263865082737];
fil = fil / sum(fil(:));
handles.smooththeta = imfilter(handles.theta,fil,'same','replicate');
handles.nframes_curr = handles.trx(handles.target).nframes;
handles.firstframe = handles.trx(handles.target).firstframe;
handles.endframe = handles.trx(handles.target).endframe;
handles.scores = nan(handles.nbehaviors,handles.nframes_curr);
handles.sex = handles.trx(handles.target).sex;
for behaviori = 1:handles.nbehaviors,
  handles.scores(behaviori,:) = min(1,max(-1,handles.trx(handles.target).(handles.scorefns{behaviori})/handles.score_norms(behaviori)));
end

if handles.t < handles.firstframe,
  handles.t = handles.firstframe;
end
if handles.t > handles.endframe,
  handles.t = handles.endframe;
end

SetCurrentFrame(handles);

%% set current frame

function SetCurrentFrame(handles,hObject)

if nargin < 2,
  hObject = nan;
end
if hObject ~= handles.slider_frame,
  set(handles.slider_frame,'Value',handles.t);
end
if hObject ~= handles.edit_frame,
  set(handles.edit_frame,'String',handles.t);
end
UpdatePlots(handles);

%% update plots

function UpdatePlots(handles,resetax)

persistent ax;

if nargin > 1 && resetax,
  ax = [];
end

trx_rad = 50;
printnone = true;
  
im = handles.readframe(handles.t);
set(handles.him,'CData',im2double(im));
i = handles.t+handles.trx(handles.target).off;
    
i0 = max(1,min(handles.nframes_curr,i-trx_rad(1)));
i1 = max(1,min(handles.nframes_curr,i+trx_rad(1)));
    
set(handles.htrx,'XData',handles.x(i0:i1),...
  'YData',handles.y(i0:i1));
    
updatefly(handles.htarget,handles.x(i),handles.y(i),...
  handles.theta(i),handles.a(i),handles.b(i));
    
for behaviori = 1:handles.nbehaviors,
  set(handles.hbar(behaviori),'YData',[0,0,1,1,0]*handles.scores(behaviori,i));
  set(handles.hbarback(behaviori),'YData',[0,0,1,1,0]*sign(handles.scores(behaviori,i)));
end
w = max(0,handles.scores(:,i)).*handles.behaviorweight;
z = sum(w);
if z == 0,
  colortmp = [0,0,0];
else
  colortmp = sum(bsxfun(@times,w,handles.behavior_colors),1)/z*max(handles.scores(:,i));
end
set(handles.htarget,'Color',colortmp);
    
tailpos = [handles.x(i)-(handles.pxback+2*handles.a(i))*cos(handles.theta(i)),...
  handles.y(i)-(handles.pxback+2*handles.a(i))*sin(handles.theta(i))];
sectori = floor(modrange((handles.smooththeta(i)-pi/8),0,2*pi)/pi*4)+1;
    
if ismember(sectori,[1,7,8]),
  halign = 'right';
elseif ismember(sectori,[2,6]),
  halign = 'center';
elseif ismember(sectori,[3,4,5]),
  halign = 'left';
else
  error('Sanity check: sectori > 8 or < 1');
end
if ismember(sectori,[1,2,3]),
  valign = 'bottom';
elseif ismember(sectori,[4,8]),
  valign = 'middle';
elseif ismember(sectori,[5,6,7]),
  valign = 'top';
else
  error('Sanity check: sectori > 8 or < 1');
end

labelis = find(handles.scores(:,i) > 0);
s = cell(1,numel(labelis));
for iii = 1:numel(labelis),
  s{iii} = sprintf('\\color[rgb]{%f %f %f}%s',handles.behavior_colors(labelis(iii),:),handles.behaviors{labelis(iii)});
end
if printnone && isempty(s),
  s = {'None'};
end
    
set(handles.htext,'String',s,...
  'Position',tailpos,...
  'HorizontalAlignment',halign,...
  'VerticalAlignment',valign);
    
set(handles.hinfo,'String',sprintf('Target %d, sex = %s',handles.target,handles.sex{i}));

for fly = 1:handles.nflies,
  if fly == handles.target || ...
      handles.t > handles.trx(fly).endframe || ...
      handles.t < handles.trx(fly).firstframe,
    set(handles.hothers(fly),'XData',nan,'YData',nan);
  else
    i = handles.t + handles.trx(fly).off;
    set(handles.hothers(fly),'XData',handles.trx(fly).x(i),...
      'yData',handles.trx(fly).y(i));
  end
end

outofbounds = handles.dozoom && (isempty(ax) || ...
  handles.x(i)-handles.pxborder < ax(1) || ...
  handles.x(i)+handles.pxborder > ax(2) || ...
  handles.y(i)-handles.pxborder < ax(3) || ...
  handles.y(i)+handles.pxborder > ax(4));
if outofbounds,
  ax = ResetAxis(handles.x,handles.y,i,handles.headerinfo.nr,handles.headerinfo.nc,handles.pxwidthradius,handles.pxheightradius,handles.pxborder);
  axis(handles.haxmain,ax);
  set(handles.hinfo,'Position',ax([1,4]));
end

%%

% --- Outputs from this function are returned to the command line.
function varargout = PlayJAABAResults_OutputFcn(hObject, eventdata, handles) 
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
v = round(get(hObject,'Value'));
if v > handles.endframe,
  v = handles.endframe;
  set(hObject,'Value',v);
elseif v < handles.firstframe,
  v = handles.firstframe;
  set(hObject,'Value',v);
end
handles.t = v;
SetCurrentFrame(handles,hObject);
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



function edit_frame_Callback(hObject, eventdata, handles)
% hObject    handle to edit_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_frame as text
%        str2double(get(hObject,'String')) returns contents of edit_frame as a double
v = str2double(get(hObject,'String'));
if isnan(v),
  set(hObject,'Value',num2str(handles.t));
  return;
end
v = max(handles.firstframe,v);
v = min(handles.endframe,v);
v = round(v);
handles.t = v;
SetCurrentFrame(handles,hObject);
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

function handles = stop(handles)

if handles.isplaying,
  pushbutton_play_Callback(handles.pushbutton_play,[],handles);
  handles = guidata(handles.pushbutton_play);
end

function playstop(handles)

pushbutton_play_Callback(handles.pushbutton_play,[],handles);


% --- Executes on button press in pushbutton_play.
function pushbutton_play_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.isplaying,
  handles.isplaying = false;
  guidata(hObject,handles);
else
  handles.isplaying = true;
  guidata(hObject,handles);
  play(handles);
end

function play(handles)

set(handles.pushbutton_play,'BackgroundColor',[.5,0,0],'String','Stop');

tic;
dt_sec = 1/handles.fps;

while true,

  if ~handles.isplaying,
    break;
  end
  
  if handles.t == handles.endframe,
    break;
  end
  
  handles = guidata(handles.pushbutton_play);
  dt = max(1,round(handles.fps*dt_sec));
  %fprintf('dt_sec = %f, dt = %d\n',dt_sec,dt);
  handles.t = handles.t+dt;
  guidata(handles.pushbutton_play,handles);
  SetCurrentFrame(handles);
  drawnow;
  dt_sec = toc;
  tic;
  if dt == 1,
    pause_time = 1/handles.fps - dt_sec;
    if pause_time > 0,
      %fprintf('pausing for %f s\n',pause_time);
      pause(pause_time);
    end
  end
  
end

set(handles.pushbutton_play,'BackgroundColor',[0,.5,0],'String','Play');
handles = guidata(handles.pushbutton_play);
handles.isplaying = false;
guidata(handles.pushbutton_play,handles);

function ax = ResetAxis(x,y,i,nr,nc,pxwidthradius,pxheightradius,pxborder)

n = numel(x);

minx = x(i);
maxx = x(i);
miny = y(i);
maxy = y(i);
w = pxwidthradius-2*pxborder;
h = pxheightradius-2*pxborder;
for j = i+1:n,
  minx = min(minx,x(j));
  maxx = max(maxx,x(j));
  dx = ceil(maxx-minx)+1;
  if dx > w,
    break;
  end
  miny = min(miny,y(j));
  maxy = max(maxy,y(j));
  dy = ceil(maxy-miny)+1;
  if dy > h,
    break;
  end
end
mux = (minx+maxx)/2;
muy = (miny+maxy)/2;
ax = [mux+pxwidthradius*[-1,1],muy+pxheightradius*[-1,1]];

if ax(2)>nc && ax(1)<1,
  ax(1:2) = (nc+1)/2+pxwidthradius*[-1,1];
elseif ax(2)>nc,
  ax(2) = nc;
  ax(1) = nc-(2*pxwidthradius+1)+1;
elseif ax(1)<1,
  ax(1) = 1;
  ax(2) = 2*pxwidthradius+1;
end
if ax(2)>nc || ax(1)<1,
  ax(1:2) = (nc+1)/2+pxwidthradius*[-1,1];
end
if ax(4)>nr && ax(3)<1,
  ax(3:4) = (nr+1)/2+pxwidthradius*[-1,1];
elseif ax(4)>nr,
  ax(4) = nr;
  ax(3) = nr-(2*pxheightradius+1)+1;
elseif ax(3)<1,
  ax(3) = 1;
  ax(4) = 2*pxheightradius+1;
end
if ax(4)>nr || ax(3)<1,
  ax(3:4) = (nr+1)/2+pxheightradius*[-1,1];
end


% --- Executes when figure_PJR is resized.
function figure_PJR_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure_PJR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'pxwidthradius'),
  return;
end

awidthradius = 3;
aheightradius = 3;
aborder = 1;
aback = .25;
meana = handles.meana;
handles.pxwidthradius = awidthradius*meana;
handles.pxheightradius = aheightradius*meana;
handles.pxborder = aborder*meana;
handles.pxback = meana*aback;

set(handles.haxmain,'Units','pixels');
handles.mainpos = get(handles.haxmain,'Position');
set(handles.haxmain,'Units','normalized');
pxheightradius1 = min((handles.pxwidthradius)*handles.mainpos(4) / handles.mainpos(3),(handles.headerinfo.nr-1)/2);
pxwidthradius1 = min((handles.pxheightradius)*handles.mainpos(3) / handles.mainpos(4),(handles.headerinfo.nc-1)/2);
handles.pxheightradius = max(pxheightradius1,handles.pxheightradius);
handles.pxwidthradius = max(pxwidthradius1,handles.pxwidthradius);
UpdatePlots(handles,true);
guidata(hObject,handles);

function fly_ButtonDownFcn(hObject,eventdata,handles,fly)

if fly == handles.target,
  return;
end

if handles.isplaying,
  handles = stop(handles);
end
handles.target = fly;
handles = SetCurrentFly(handles);
guidata(hObject,handles);


% --------------------------------------------------------------------
function menu_options_Callback(hObject, eventdata, handles)
% hObject    handle to menu_options (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to menu_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function menu_zoom_showwholevideo_Callback(hObject, eventdata, handles)
% hObject    handle to menu_zoom_showwholevideo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.dozoom = strcmpi(get(hObject,'Checked'),'on');
guidata(hObject,handles);
UpdateZoom(handles);

function UpdateZoom(handles)

if handles.dozoom,
  set(handles.menu_zoom_showwholevideo,'Checked','off');
  set(handles.menu_zoom_keepflyinaxes,'Checked','on');
  UpdatePlots(handles,true);
else
  set(handles.menu_zoom_showwholevideo,'Checked','on');
  set(handles.menu_zoom_keepflyinaxes,'Checked','off');
  axis(handles.haxmain,'image');
  set(handles.hinfo,'Position',[1,handles.headerinfo.nr]);
end


% --------------------------------------------------------------------
function menu_zoom_keepflyinaxes_Callback(hObject, eventdata, handles)
% hObject    handle to menu_zoom_keepflyinaxes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.dozoom = strcmpi(get(hObject,'Checked'),'off');
guidata(hObject,handles);
UpdateZoom(handles);


% --- Executes during object deletion, before destroying properties.
function figure_PJR_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure_PJR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

CleanUp(handles);

function CleanUp(handles)

if isfield(handles,'fid'),
  fclose(handles.fid);
end


% --- Executes when user attempts to close figure_PJR.
function figure_PJR_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_PJR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);


% --------------------------------------------------------------------
function menu_playbackspeed_Callback(hObject, eventdata, handles)
% hObject    handle to menu_playbackspeed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

while true,
  
  v = inputdlg('Playback speed (fps):','Playback speed',1,{num2str(handles.fps)});
  v = v{1};
  if ~ischar(v),
    return;
  end
  v = str2double(v);
  if isnan(v) || v <= 0,
    continue;
  end
  break;
end
handles.fps = v;
guidata(hObject,handles);


% --- Executes on key press with focus on figure_PJR and none of its controls.
function figure_PJR_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure_PJR (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Key,
  
  case 'leftarrow',
    if ~handles.isplaying,
      handles.t = max(handles.t-1,handles.firstframe);
      guidata(hObject,handles);
      SetCurrentFrame(handles);
    end
  case 'rightarrow',
    if ~handles.isplaying,
      handles.t = min(handles.t+1,handles.endframe);
      guidata(hObject,handles);
      SetCurrentFrame(handles);
    end  
  case 'uparrow',
    if ~handles.isplaying,
      handles.t = max(handles.t-30,handles.firstframe);
      guidata(hObject,handles);
      SetCurrentFrame(handles);
    end  
  case 'downarrow',
    if ~handles.isplaying,
      handles.t = min(handles.t+30,handles.endframe);
      guidata(hObject,handles);
      SetCurrentFrame(handles);
    end
  case {'space','p'},
    playstop(handles);
end


function RecursiveSetKeyPressFcn(hfig)

hchil = findall(hfig,'-property','KeyPressFcn');
goodidx = true(1,numel(hchil));
for i = 1:numel(hchil),
  if strcmpi(get(hchil(i),'Type'),'uicontrol') && strcmpi(get(hchil(i),'Style'),'edit'),
    goodidx(i) = false;
  end
end
set(hchil(goodidx),'KeyPressFcn',get(hfig,'KeyPressFcn'));


% --------------------------------------------------------------------
function menu_switchtarget_Callback(hObject, eventdata, handles)
% hObject    handle to menu_switchtarget (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

s = cellstr(num2str((1:handles.nflies)'));
[sel,ok] = listdlg('ListString',s,'SelectionMode','single',...
  'InitialValue',handles.target,'Name','Switch to target',...
  'PromptString','Switch to target: ');
if ~ok,
  return;
end
handles.target = sel;
handles = SetCurrentFly(handles);
guidata(hObject,handles);
