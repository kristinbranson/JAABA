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

% Last Modified by GUIDE v2.5 08-Mar-2013 10:40:15

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

handles.cropSz = 201;

[JLabelH,expnum,fly,t] = myparse(varargin,...
  'JLabelH',[],...
  'expnum',0,...
  'fly',0,...
  't',0);

handles.JLabelH = JLabelH;
handles.data = handles.JLabelH.guidata.data;
handles.JLabelH.guidata.NJObj.SetCompareFramesHandle(hObject);

handles.expnum = expnum;
handles.fly = fly;
handles.t = t;
handles.theta = handles.data.GetTrxValues('Theta1',expnum,fly,t);


set(handles.figure1,'Pointer','watch');

set(handles.edit_ignore,'String','0'); 
set(handles.radiobutton_behavior,'String',handles.data.labelnames{1});
set(handles.popupmenu_jump,'String',{'Current Fly','Current Experiment','All experiments','Training Data'},'Value',1);

handles.data.SetCurrentFlyForBag(handles.expnum,handles.fly,handles.t);
handles = CacheFrames(handles);
handles = initialize(handles);

set(handles.figure1,'Pointer','arrow');

fgColor = get(handles.text_status,'ForegroundColor');
bgColor = get(handles.text_status,'BackgroundColor');
set([handles.uipanel_jump,handles.radiobutton_all,...
    handles.radiobutton_behavior,handles.radiobutton_none,...
    handles.radiobutton_shortcut],....
  'BackgroundColor',bgColor,'ForegroundColor',fgColor);
set(handles.radiobutton_all,'Value',1);
handles.jump_restrict = 'all';

% Make sure the colors are OK on Mac
adjustColorsIfMac(hObject);

% Update handles structure
guidata(hObject, handles);
JLabel('UpdatePrediction',handles.JLabelH);

% UIWAIT makes CompareFrames wait for user response (see UIRESUME)
% uiwait(handles.figure1);
return


% --- Outputs from this function are returned to the command line.
function varargout = CompareFrames_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
return


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
handles.pfList = handles.data.allperframefns;

handles.pfcache = {};
for ndx = 1:numel(handles.data.perframedata)
  handles.pfcache{ndx} = handles.data.perframedata{ndx}(handles.ts-firstFrame+1);
end
set(handles.popupmenu_pf,'String',handles.pfList);
set(handles.text_info,'String',sprintf('Animal:%d, Frame:%d, Exp-Number:%d Exp-Name:%s',handles.fly,handles.t,handles.expnum,handles.data.expnames{handles.expnum}));
return


function handles = initialize(handles)

handles.align = get(handles.radiobutton_align,'Value');

handles.curFrame = handles.centralframe;
handles.himage_preview = imagesc(0,'Parent',handles.axes_preview,[0,255]);
set(handles.himage_preview,'HitTest','off');
axis(handles.axes_preview,'equal');
axis(handles.axes_preview,'off');
hold(handles.axes_preview,'on');
colormap(handles.axes_preview,'gray');

handles.trx_plot = plot(handles.axes_preview,nan,nan,'-o','MarkerSize',3);
handles.fly_plot = plot(handles.axes_preview,nan,nan,'r','LineWidth',2);
handles.fly_plot_extra = plot(handles.axes_preview,nan,nan);

handles.timeline = imagesc(0,'Parent',handles.axes_timeline,[0,255]);
set(handles.timeline,'HitTest','off');
fgColor = get(handles.text_status,'ForegroundColor');
set(handles.axes_timeline,'XColor',fgColor,'YColor',fgColor);
hold(handles.axes_timeline,'on');
handles.centerline = plot(handles.axes_timeline,[nan nan],[0.5 3.5],'g');
handles.curline = plot(handles.axes_timeline,[nan nan],[0.5 3.5],'y');

handles.pf_plot = plot(handles.axes_pf,nan,nan,'-o','MarkerSize',3);
set(handles.axes_pf,'Color',get(handles.text_status,'BackgroundColor'),'XColor',fgColor,'YColor',fgColor);
set(handles.pf_plot,'Color',fgColor);

handles.jumpList = [];
handles.jumpNum = 0;
handles.maxFrame = numel(handles.ts);
handles.minFrame = 1;

UpdatePlots(handles,handles.centralframe);
UpdateTimeline(handles);
popupmenu_pf_Callback(handles.popupmenu_pf,[],handles);




function handles = Jump(handles,jtype)
sel = get(handles.popupmenu_jump,'Value');
if handles.jumpNum<1,
  curT = handles.t;
  curExp = handles.expnum;
  curFly = handles.fly;
else
  curT = handles.jumpList(handles.jumpNum).t;
  curFly = handles.jumpList(handles.jumpNum).fly;
  curExp = handles.jumpList(handles.jumpNum).exp;
end

ignore = str2double(get(handles.edit_ignore,'String'));
if strcmp(jtype,'prev'),
  if (handles.jumpNum ==0),
    return;
  end
  handles.jumpNum = handles.jumpNum - 1;
  if handles.jumpNum == 0
    pushbutton_reset_Callback(handles.figure1, [], handles)
    return;
  end
  JLabel('SetCurrentMovie',handles.JLabelH, handles.jumpList(handles.jumpNum).exp);
  JLabel('SetCurrentFlies',handles.JLabelH, handles.jumpList(handles.jumpNum).fly);
  JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
    handles.jumpList(handles.jumpNum).t,handles.JLabelH.figure_JLabel);
  
else
  
  nextT = curT; newFly = curFly; newExp = curExp;
  
  jumpList = handles.jumpList(1:handles.jumpNum);
  jumpList(end+1).exp = handles.expnum;
  jumpList(end).fly = handles.fly;
  jumpList(end).t = handles.t;
  jumpList(end).dist = 0;
  jumpList(end).theta = handles.theta;
  if sel == 1
    
    [nextT distT] = handles.data.NextClosestBagFly(...
      jtype,curT,curExp,curFly,[],ignore,jumpList,handles.jump_restrict);
    if isempty(nextT), return; end
    JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
      nextT,handles.JLabelH.figure_JLabel);
    
  elseif sel ==2
    
    [nextT distT,newFly] = handles.data.NextClosestBagExp(...
      jtype,curT,curExp,curFly,ignore,jumpList,handles.jump_restrict);
    if isempty(nextT), return; end
    JLabel('SetCurrentFlies',handles.JLabelH, newFly);
    JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
      nextT,handles.JLabelH.figure_JLabel);
    
  elseif sel == 3
    
  elseif sel == 4
    [nextT distT,newFly,newExp] = handles.data.NextClosestBagTraining(...
      jtype,curT,curExp,curFly,ignore,jumpList,handles.jump_restrict);
    if isempty(nextT), return; end
    JLabel('SetCurrentMovie',handles.JLabelH, newExp);
    JLabel('SetCurrentFlies',handles.JLabelH, newFly);
    JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
      nextT,handles.JLabelH.figure_JLabel);
    
    
  end
  handles.jumpNum = handles.jumpNum + 1;
  handles.jumpList(handles.jumpNum).exp = newExp;
  handles.jumpList(handles.jumpNum).fly = newFly;
  handles.jumpList(handles.jumpNum).t = nextT;
  handles.jumpList(handles.jumpNum).dist = distT;
  handles.jumpList(handles.jumpNum).theta = handles.data.GetTrxValues('Theta1',newExp,newFly,nextT);

end

handles.curFrame = handles.centralframe;
UpdatePlots(handles,handles.centralframe);


function matchZoom(handles)
sz = get(handles.JLabelH.axes_preview,'Position');
xlim = get(handles.JLabelH.axes_preview,'XLim');
ylim = get(handles.JLabelH.axes_preview,'YLim');
scalex = sz(3)/(xlim(2)-xlim(1));
scaley = sz(4)/(ylim(2)-ylim(1));

scale = min(scalex,scaley);
this_sz = get(handles.axes_preview,'Position');
xsz = this_sz(3)/scale;
ysz = this_sz(4)/scale;

x = handles.cropSz/2;
y = handles.cropSz/2;
xlim = [x-xsz/2,x+xsz/2];
ylim = [y-ysz/2,y+ysz/2];

set(handles.axes_preview,'XLim',xlim,'YLim',ylim);


function UpdatePlots(handles,t)

img = handles.imgcache(:,:,:,t);
theta = handles.trxcache.theta(handles.centralframe);
x = handles.trxcache.x(handles.centralframe);
y = handles.trxcache.y(handles.centralframe);

if handles.jumpNum < 1,
  theta_jlabel = handles.theta;
else
  theta_jlabel = handles.jumpList(handles.jumpNum).theta;
end

smsz = (handles.cropSz-1)/2;
[h,w,d] = size(img);
padImg = zeros(h+2*smsz,w+2*smsz,d);
padImg(smsz + (1:h),smsz+(1:w),:) = img;
timg = padImg(round(y) + smsz+ (-smsz:smsz),round(x)+ smsz + (-smsz:smsz),:);

% tform = maketform('affine',[1 0 0; 0 1 0; xtr ytr 1]);
% timg = imtransform(img,tform,'XData',[1 size(img,2)],'YData',[1 size(img,1)]);

if handles.align,
%   dtheta = theta-theta_jlabel;
%   timg = imrotate(timg,dtheta*180/pi,'bilinear','crop');
%   rotmat = [cos(-dtheta) sin(-dtheta); -sin(-dtheta) cos(-dtheta)];
%   temp = [(handles.trxcache.x-x)' (handles.trxcache.y-y)']*rotmat;
%   rotatedtrxx = temp(:,1)+size(img,2)/2; 
%   rotatedtrxy = temp(:,2)+size(img,1)/2;
%   set(handles.trx_plot,'XData',rotatedtrxx,'Ydata',rotatedtrxy);
%   pos.x = rotatedtrxx(t);
%   pos.y = rotatedtrxy(t);
%   pos.theta = handles.trxcache.theta(t)-dtheta;
%   pos.a = handles.trxcache.a(t);
%   pos.b = handles.trxcache.b(t);

  dtheta = theta-theta_jlabel;
  timg = imrotate(timg,dtheta*180/pi,'bilinear','crop');
  rotmat = [cos(-dtheta) sin(-dtheta); -sin(-dtheta) cos(-dtheta)];
  temp = [(handles.trxcache.x-x)' (handles.trxcache.y-y)']*rotmat;
  rotatedtrxx = temp(:,1)+smsz+0.5; 
  rotatedtrxy = temp(:,2)+smsz+0.5;
  set(handles.trx_plot,'XData',rotatedtrxx,'Ydata',rotatedtrxy);
  pos.x = rotatedtrxx(t);
  pos.y = rotatedtrxy(t);
  pos.theta = handles.trxcache.theta(t)-dtheta;
  pos.a = handles.trxcache.a(t);
  pos.b = handles.trxcache.b(t);
else
  set(handles.trx_plot,'XData',handles.trxcache.x-x+smsz+0.5,...
                       'Ydata',handles.trxcache.y-y+smsz+0.5);  
  pos.x = handles.trxcache.x(t) - handles.trxcache.x(handles.centralframe)+smsz;
  pos.y = handles.trxcache.y(t) - handles.trxcache.y(handles.centralframe)+smsz;
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
if handles.jumpNum > 0
  status_str = sprintf('Jump Order:%d \n Distance:%.2f', ...
    handles.jumpNum,handles.jumpList(handles.jumpNum).dist);
else
  status_str = sprintf('Jump Order:0 \n Distance:0');
end

set(handles.text_status,'String',status_str);





function UpdateTimeline(handles)
timg = zeros(3,numel(handles.ts),3);

for behaviori = 1:handles.data.nbehaviors
  curcolor = handles.JLabelH.guidata.labelcolors(behaviori,:);
  idx = handles.labelsCache.vals == behaviori;
  
  idxPred = handles.predCache.predictedidx == behaviori;
  
  for channel = 1:3
    timg(1,idx,channel) = curcolor(channel);
    timg(2,idxPred,channel) = curcolor(channel);
    scoreNdx = ceil(handles.scoresCache(idxPred)*31)+32;
    timg(3,idxPred,channel) = handles.JLabelH.guidata.scorecolor(scoreNdx,channel,1);
  end
  
  
end

set(handles.timeline,'CData',timg);
set(handles.axes_timeline,'XLim',[1 numel(handles.ts)]);
set(handles.axes_pf,'XLim',[1 numel(handles.ts)]);
set(handles.axes_timeline,'YLim',[0.5 3.5]);
set(handles.axes_timeline,'YTick',[1 2 3],'YTickLabel',{'Manual','Predictions','Scores'});

xticks = 0:10:handles.centralframe;
xticks = [-fliplr(xticks(2:end)) xticks] + handles.centralframe;
xticks(xticks<0) = [];

xticks(xticks>numel(handles.ts))= [];
xtickLabel = {};
for ndx = 1:numel(xticks)
  xtickLabel{end+1} = sprintf('%d',xticks(ndx) - handles.centralframe + handles.t);
end
set(handles.axes_timeline,'XTick',xticks,'XTickLabel',xtickLabel);
set(handles.axes_pf,'Xtick',xticks,'XTickLabel',{});


function check = CheckBag(handles)
check = false;
  if isempty(handles.data.bagModels),
    uiwait(warndlg('You need to bag!'));
    return;
  end
  
  if handles.data.fastPredictBag.ts < handles.data.classifierTS
    uiwait(warndlg('Bagging was done before the current classifier was trained. Things maybe out of sync'));
    return;
  end
  
check = true;


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
pushbutton_reset_Callback(hObject, eventdata, handles);


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
if handles.jumpNum < 1,
  jlabelframe = handles.t;
else
  jlabelframe = handles.jumpList(handles.jumpNum).t;
end
guidata(hObject,handles);

lastFrame = curFrame;
while true
  handles = guidata(hObject);
  if handles.play 
    ticker = tic;
    UpdatePlots(handles,curFrame);
    if get(handles.radiobutton_sync,'Value');
      
      JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
        jlabelframe+curFrame-handles.centralframe,handles.JLabelH.figure_JLabel);
    end
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


% --- Executes on button press in pushbutton_setcurrent.
function pushbutton_setcurrent_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setcurrent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'Pointer','watch');
drawnow();
expnum = handles.JLabelH.guidata.expi;
fly = handles.JLabelH.guidata.flies;
t = handles.JLabelH.guidata.ts;

handles.expnum = expnum;
handles.fly = fly;
handles.t = t;
handles.theta = handles.data.GetTrxValues('Theta1',expnum,fly,t);

handles.data.SetCurrentFlyForBag(handles.expnum,handles.fly,handles.t)

handles = CacheFrames(handles);
handles = initialize(handles);

guidata(hObject,handles);
pushbutton_reset_Callback(hObject, eventdata, handles)
set(handles.figure1,'Pointer','arrow');



% --- Executes on key press with focus on figure1 or any of its controls.
function figure1_WindowKeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'Pointer','watch');
drawnow();

switch eventdata.Key
  
  case 'leftarrow'
    handles.curFrame = handles.curFrame - 1;
    if handles.curFrame < 1
      set(handles.figure1,'Pointer','arrow');
      return;
    end
    if handles.jumpNum < 1
      jlabelframe = handles.t + handles.curFrame - handles.centralframe;
    else
      jlabelframe = handles.jumpList(handles.jumpNum).t + handles.curFrame - handles.centralframe;;
    end
    UpdatePlots(handles,handles.curFrame);
    if get(handles.radiobutton_sync,'Value');
      JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
      jlabelframe,handles.JLabelH.figure_JLabel);
    end
    
  case 'rightarrow'
    handles.curFrame = handles.curFrame + 1;
    if handles.curFrame > numel(handles.ts)
      set(handles.figure1,'Pointer','arrow');
      return;
    end
    UpdatePlots(handles,handles.curFrame);
    if get(handles.radiobutton_sync,'Value');
      if handles.jumpNum < 1
        jlabelframe = handles.t + handles.curFrame - handles.centralframe;
      else
        jlabelframe = handles.jumpList(handles.jumpNum).t + handles.curFrame - handles.centralframe;;
      end
      JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
        jlabelframe,handles.JLabelH.figure_JLabel);
    end    
  case 'uparrow'
    if CheckBag(handles)
      handles = Jump(handles,'next');
      JLabel('UpdatePrediction',handles.JLabelH);
    end
  case 'downarrow'
    if CheckBag(handles)
      handles = Jump(handles,'prev');
      JLabel('UpdatePrediction',handles.JLabelH);
    end
end
guidata(hObject,handles);
set(handles.figure1,'Pointer','arrow');


% --- Executes on button press in pushbutton_Bag.
function pushbutton_Bag_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_Bag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'Pointer','watch');
drawnow();
handles.data.DoFastBagging();
handles.data.SetCurrentFlyForBag(handles.expnum,handles.fly,handles.t)
set(handles.figure1,'Pointer','arrow');
pushbutton_reset_Callback(hObject, eventdata, handles)


function edit_ignore_Callback(hObject, eventdata, handles)
% hObject    handle to edit_ignore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_ignore as text
%        str2double(get(hObject,'String')) returns contents of edit_ignore as a double


val = str2double(get(hObject,'String'));
if isempty(val) || (round(val)-val)~=0 || val < 0,
  uiwait(warndlg('Frames to ignore should be 0 or positive integer'));
  set(hObject,'String','0');
end
pushbutton_reset_Callback(hObject,eventdata,handles);

% --- Executes during object creation, after setting all properties.
function edit_ignore_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_ignore (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_reset.
function pushbutton_reset_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
JLabel('SetCurrentMovie',handles.JLabelH, handles.expnum);
JLabel('SetCurrentFlies',handles.JLabelH, handles.fly);
JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
  handles.t,handles.JLabelH.figure_JLabel);
handles.jumpList = [];
handles.jumpNum = 0;
handles.curFrame = handles.centralframe;

guidata(hObject,handles);
UpdatePlots(handles,handles.centralframe);
JLabel('UpdatePrediction',handles.JLabelH);


% --- Executes on button press in radiobutton_sync.
function radiobutton_sync_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_sync (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_sync
if handles.jumpNum > 0,
  jlabelframe = handles.jumpList(handles.jumpNum).t;
else
  jlabelframe = handles.t;
end

if ~get(hObject,'Value');
  JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
    jlabelframe,handles.JLabelH.figure_JLabel);
else
  JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
    jlabelframe + handles.curFrame-handles.centralframe,handles.JLabelH.figure_JLabel);
  
end
guidata(hObject,handles);

% --- Executes on button press in radiobutton_align.
function radiobutton_align_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_align (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_align
handles.align = get(hObject,'Value');
guidata(hObject,handles);
UpdatePlots(handles,handles.curFrame);


% --- Executes on button press in pushbutton_central.
function pushbutton_central_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_central (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.curFrame = handles.centralframe;
UpdatePlots(handles,handles.curFrame);
if get(handles.radiobutton_sync,'Value');
  if handles.jumpNum < 1
    jlabelframe = handles.t;
  else
    jlabelframe = handles.jumpList(handles.jumpNum).t;
  end
  JLabel('SetCurrentFrame',handles.JLabelH, 1, ...
    jlabelframe,handles.JLabelH.figure_JLabel);
end
guidata(hObject,handles);

% --- Executes on selection change in popupmenu_pf.
function popupmenu_pf_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_pf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_pf contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_pf
ndx = get(hObject,'Value');
ylim(1) = min(handles.pfcache{ndx});
ylim(2) = max(handles.pfcache{ndx});
set(handles.axes_pf,'YLim',ylim);
set(handles.pf_plot,'YData',handles.pfcache{ndx},'XData',1:numel(handles.ts));

% --- Executes during object creation, after setting all properties.
function popupmenu_pf_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_pf (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel_jump.
function uipanel_jump_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_jump 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

set(handles.figure1,'Pointer','watch');
drawnow();

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'radiobutton_all'
      handles.jump_restrict = 'all';
    case 'radiobutton_behavior'
      handles.jump_restrict = 'behavior';
    case 'radiobutton_none'
      handles.jump_restrict = 'none';
end

pushbutton_reset_Callback(hObject,eventdata,handles);
guidata(hObject,handles);
set(handles.figure1,'Pointer','arrow');


% --- Executes on button press in radiobutton_shortcut.
function radiobutton_shortcut_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton_shortcut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton_shortcut
if get(hObject,'Value')
  handles.JLabelH.guidata.NJObj.SetCurrentType('Jump To Similar Frames');
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.UnsetCurrentFlyForBag();
JLabel('UpdatePrediction',handles.JLabelH);

% Hint: delete(hObject) closes the figure
delete(hObject);


% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
