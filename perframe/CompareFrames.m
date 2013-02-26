function varargout = CompareFrames(varargin)
% COMPAREFRAMES MATLAB code for CompareFrames.fig
%      COMPAREFRAMES, by itself, creates a new COMPAREFRAMES or raises the existing
%      singleton*.
%
%      H = COMPAREFRAMES returns the handle to a new COMPAREFRAMES or the handle to
%      the existing singleton*.
%
%      COMPAREFRAMES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COMPAREFRAMES.M with the given input arguments.
%
%      COMPAREFRAMES('Property','Value',...) creates a new COMPAREFRAMES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CompareFrames_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CompareFrames_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CompareFrames

% Last Modified by GUIDE v2.5 26-Feb-2013 13:40:05

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CompareFrames_OpeningFcn, ...
                   'gui_OutputFcn',  @CompareFrames_OutputFcn, ...
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


% --- Executes just before CompareFrames is made visible.
function CompareFrames_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CompareFrames (see VARARGIN)

% Choose default command line output for CompareFrames
handles.output = hObject;

[JLabelH,expnum,fly,t] = myparse(varargin,...
  'JLabelH',[],...
  'expnum',0,...
  'fly',0,...
  't',0);

handles.JLabelH = JLabelH;
handles.data = handles.JLabelH.guidata.data;

handles.expnum = expnum;
handles.fly = fly;
handles.t = t;

handles = CacheFrames(handles);
handles = initialize(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CompareFrames wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CompareFrames_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function handles = CacheFrames(handles)
movie_filename = handles.JLabelH.guidata.movie_filename;
readframe_fcn = get_readframe_fcn(movie_filename);

firstFrame = handles.data.GetTrxFirstFrame(handles.expnum,handles.fly);
lastFrame = handles.data.GetTrxEndFrame(handles.expnum,handles.fly);

numcache = 30;

ts = max(firstFrame,handles.t-numcache):min(lastFrame,handles.t+numcache);
handles.ts = ts;

handles.centralframe = find(ts == handles.t);

img = readframe_fcn(ts(1));
handles.imgcache = zeros(handles.JLabelH.guidata.movie_height,...
  handles.JLabelH.guidata.movie_width,size(img,3),numel(ts));

for ndx = 1:numel(ts)
  img = readframe_fcn(ts(ndx));
  handles.imgcache(:,:,:,ndx) = img;
end

handles.trxcache = handles.data.GetTrxValues('Trx',handles.expnum,handles.fly,ts);


handles.labelsCache = handles.data.GetLabelIdx(handles.expnum,handles.fly,ts(1),ts(end));
handles.predCache = handles.data.GetPredictedIdx(handles.expnum,handles.fly,ts(1),ts(end));
handles.scoresCache = handles.data.NormalizeScores(handles.predCache.scoresidx);

function handles = initialize(handles)

set(handles.popupmenu_jump,'String',{'Current Fly','Current Experiment','All experiments'},'Value',1);

handles.curFrame = handles.centralframe;
handles.himage_preview = imagesc(0,'Parent',handles.axes_preview,[0,255]);
set(handles.himage_preview,'HitTest','off');
axis(handles.axes_preview,'equal');
axis(handles.axes_preview,'off');
hold(handles.axes_preview,'on');
colormap(handles.axes_preview,'gray');

handles.theta_jlabel = handles.data.GetTrxValues('Theta1',...
  handles.JLabelH.guidata.expi,handles.JLabelH.guidata.flies,...
  handles.JLabelH.guidata.ts);

handles.trx_plot = plot(handles.axes_preview,nan,nan,'-o','MarkerSize',3);
handles.fly_plot = plot(handles.axes_preview,nan,nan);
handles.fly_plot_extra = plot(handles.axes_preview,nan,nan);

handles.timeline = imagesc(0,'Parent',handles.axes_timeline,[0,255]);
set(handles.timeline,'HitTest','off');
hold(handles.axes_timeline,'on');
handles.centerline = plot(handles.axes_timeline,[nan nan],[0.5 3.5],'g');
handles.curline = plot(handles.axes_timeline,[nan nan],[0.5 3.5],'y');

UpdatePlots(handles,handles.centralframe);
handles.maxFrame = numel(handles.ts);
handles.minFrame = 1;
UpdateTimeline(handles);




function ChangeFly(handles)



function Jump(handles,jtype)
sel = get(handles.popupmenu_jump,'Value');
curT = handles.JLabelH.guidata.ts;
curFly = handles.data.flies;
curExp = handles.data.expi;

if sel == 1
  
  [nextT distT] = handles.data.NextClosestBagFly(jtype,curT,curExp,curFly);
  if isempty(nextT), return; end
  JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
    nextT,handles.JLabelH.figure_JLabel);

elseif sel ==2

  [nextT distT,newfly] = handles.data.NextClosestBagExp(jtype,curT,curExp);
  if isempty(nextT), return; end
  JLabel('SetCurrentFlies',handles.JLabelH, newfly);
  JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
    nextT,handles.JLabelH.figure_JLabel);

else
  
end



function align()


function matchZoom(handles)
sz = get(handles.JLabelH.axes_preview,'Position');
xlim = get(handles.JLabelH.axes_preview,'XLim');
scale = sz(3)/(xlim(2)-xlim(1));
this_sz = get(handles.axes_preview,'Position');
xsz = this_sz(3)/scale;
ysz = this_sz(4)/scale;

isz = size(handles.imgcache(:,:,:,1));
x = isz(2)/2;
y = isz(1)/2;
xlim = [x-xsz/2,x+xsz/2];
ylim = [y-ysz/2,y+ysz/2];

set(handles.axes_preview,'XLim',xlim,'YLim',ylim);


function UpdatePlots(handles,t)

img = handles.imgcache(:,:,:,t);
theta = handles.trxcache.theta(handles.centralframe);
x = handles.trxcache.x(handles.centralframe);
y = handles.trxcache.y(handles.centralframe);

theta_jlabel = handles.theta_jlabel;

xtr = size(img,2)/2 - x;
ytr = size(img,1)/2 - y;
tform = maketform('affine',[1 0 0; 0 1 0; xtr ytr 1]);
timg = imtransform(img,tform,'XData',[1 size(img,2)],'YData',[1 size(img,1)]);

if handles.align,
  dtheta = theta-theta_jlabel;
  timg = imrotate(timg,dtheta*180/pi,'bilinear','crop');
  rotmat = [cos(-dtheta) sin(-dtheta); -sin(-dtheta) cos(-dtheta)];
  temp = [(handles.trxcache.x-x)' (handles.trxcache.y-y)']*rotmat;
  rotatedtrxx = temp(:,1)+size(img,2)/2; 
  rotatedtrxy = temp(:,2)+size(img,1)/2;
  set(handles.trx_plot,'XData',rotatedtrxx,'Ydata',rotatedtrxy);
  pos.x = rotatedtrxx(t);
  pos.y = rotatedtrxy(t);
  pos.theta = handles.trxcache.theta(t)-dtheta;
  pos.a = handles.trxcache.a(t);
  pos.b = handles.trxcache.b(t);
else
  set(handles.trx_plot,'XData',handles.trxcache.x-x+size(img,2)/2,...
                       'Ydata',handles.trxcache.y-y+size(img,1)/2);  
  pos.x = handles.trxcache.x(t);
  pos.y = handles.trxcache.y(t);
  pos.theta = handles.trxcache.theta(t);
  pos.a = handles.trxcache.a(t);
  pos.b = handles.trxcache.b(t);
end

set(handles.himage_preview,'CData',uint8(timg));

set(handles.curline,'XData',[t t]);
set(handles.centerline,'XData',[handles.centralframe handles.centralframe]);

UpdateTargetPosition(handles.data.targettype,handles.fly_plot,...
  handles.fly_plot_extra,pos);
matchZoom(handles);


function UpdateTimeline(handles)
timg = zeros(3,numel(handles.ts),3);

for behaviori = 1:handles.data.nbehaviors
  curcolor = handles.JLabelH.guidata.labelcolors(behaviori,:);
  idx = handles.labelsCache.vals == behaviori;
  
  idxPred = handles.predCache.predictedidx == behaviori;
  
  for channel = 1:3
    timg(1,idx,channel) = curcolor(channel);
    timg(3,idxPred,channel) = curcolor(channel);
    scoreNdx = ceil(handles.scoresCache(idxPred)*31)+32;
    timg(2,idxPred,channel) = handles.JLabelH.guidata.scorecolor(scoreNdx,channel,1);
  end
  
  
end

set(handles.timeline,'CData',timg);
set(handles.axes_timeline,'XLim',[1 numel(handles.ts)]);
set(handles.axes_timeline,'YLim',[0.5 3.5]);
set(handles.axes_timeline,'YTick',[1 2 3],'YTickLabel',{'Manual','Scores','Predictions'});

xticks = 0:10:handles.centralframe;
xticks = [-fliplr(xticks(2:end)) xticks] + handles.centralframe;
xticks(xticks<0) = [];

xticks(xticks>numel(handles.ts))= [];
xtickLabel = {};
for ndx = 1:numel(xticks)
  xtickLabel{end+1} = sprintf('%d',xticks(ndx) - handles.centralframe + handles.t);
end
set(handles.axes_timeline,'XTick',xticks,'XTickLabel',xtickLabel);


% --- Executes on button press in pushbutton_close.
function pushbutton_close_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in popupmenu_jump.
function popupmenu_jump_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_jump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_jump contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_jump


% --- Executes during object creation, after setting all properties.
function popupmenu_jump_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_jump (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_curexp.
function pushbutton_curexp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_curexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in pushbutton_play.
function pushbutton_play_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_play (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pushbutton_play

if get(hObject,'Value')==1
  handles.play = true;
else
  handles.play = false;
end
handles.play_FPS = handles.JLabelH.guidata.play_FPS;
guidata(hObject,handles);
play(hObject);


function play(hObject)

handles = guidata(hObject);
curFrame = handles.curFrame;
jlabelframe = handles.JLabelH.guidata.ts;
handles.theta_jlabel = handles.data.GetTrxValues('Theta1',...
  handles.JLabelH.guidata.expi,handles.JLabelH.guidata.flies,...
  handles.JLabelH.guidata.ts+handles.centralframe-curFrame);
guidata(hObject,handles);

lastFrame = curFrame;
while true
  handles = guidata(hObject);
  if handles.play 
    ticker = tic;
    UpdatePlots(handles,curFrame);
    JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
      jlabelframe+curFrame-handles.curFrame,handles.JLabelH.figure_JLabel);
    lastFrame = curFrame;
    dt_sec = toc(ticker);
    dt = dt_sec*handles.play_FPS;
    curFrame = curFrame + ceil(dt);
    if curFrame > handles.maxFrame
      curFrame = handles.minFrame;
    end
    pause_time = (ceil(dt)-dt)/handles.play_FPS - dt_sec;
    if pause_time < 0,
      drawnow;
    else
      pause(pause_time);
    end
  else
    break;
  end
end
handles.curFrame = lastFrame;
guidata(hObject,handles);
  
  
% --- Executes on selection change in popupmenu_align.
function popupmenu_align_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_align contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_align


% --- Executes during object creation, after setting all properties.
function popupmenu_align_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_curfly.
function pushbutton_curfly_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_curfly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton_allexp.
function pushbutton_allexp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_allexp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in align.
function align_Callback(hObject, eventdata, handles)
% hObject    handle to align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of align


% --- Executes on button press in pushbutton_setcurrent.
function pushbutton_setcurrent_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

switch eventdata.Key
  
  case 'leftarrow'
    handles.curFrame = handles.curFrame - 1;
    if handles.curFrame < 1
      return;
    end
    jlabelframe = handles.JLabelH.guidata.ts-1;
    UpdatePlots(handles,handles.curFrame);
    JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
      jlabelframe,handles.JLabelH.figure_JLabel);
    
  case 'rightarrow'
    handles.curFrame = handles.curFrame + 1;
    if handles.curFrame > numel(handles.ts)
      return;
    end
    jlabelframe = handles.JLabelH.guidata.ts+1;
    UpdatePlots(handles,handles.curFrame);
    JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
      jlabelframe,handles.JLabelH.figure_JLabel);
    
  case 'uparrow'
    Jump(handles,'next');
  case 'downarrow'
    Jump(handles,'prev');
  
end
guidata(hObject,handles);


% --- Executes on button press in pushbutton_Bag.
function pushbutton_Bag_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Bag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.data.DoFastBagging();
handles.data.SetCurrentFlyForBag(handles.expnum,handles.fly,handles.t)