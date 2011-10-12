function varargout = showSimilarFrames(varargin)
% SHOWSIMILARFRAMES MATLAB code for showSimilarFrames.fig
%      SHOWSIMILARFRAMES, by itself, creates a new SHOWSIMILARFRAMES or raises the existing
%      singleton*.
%
%      H = SHOWSIMILARFRAMES returns the handle to a new SHOWSIMILARFRAMES or the handle to
%      the existing singleton*.
%
%      SHOWSIMILARFRAMES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOWSIMILARFRAMES.M with the given input arguments.
%
%      SHOWSIMILARFRAMES('Property','Value',...) creates a new SHOWSIMILARFRAMES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before showSimilarFrames_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to showSimilarFrames_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help showSimilarFrames

% Last Modified by GUIDE v2.5 05-Oct-2011 14:07:42

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @showSimilarFrames_OpeningFcn, ...
                   'gui_OutputFcn',  @showSimilarFrames_OutputFcn, ...
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


% --- Executes just before showSimilarFrames is made visible.
function showSimilarFrames_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to showSimilarFrames (see VARARGIN)

% Choose default command line output for showSimilarFrames
handles.output = hObject;

handles.maxFrames = 20;
handles.numSimilar = 4;
handles.halfSize = 30;
handles.imgX = 1024;
handles.imgY = 1024;
handles.frameNo = handles.maxFrames+1;
handles.isPlaying = false;
handles.fwdLimit = 5;
handles.revLimit = 5;
handles.fps = 3;
handles.maxfps = 10;

set(handles.frameSlider,'SliderStep',[1/(2*handles.maxFrames) 3/(2*handles.maxFrames)]);
set(handles.beforeSlider,'SliderStep',[1/handles.maxFrames 3/handles.maxFrames]);
set(handles.afterSlider,'SliderStep',[1/handles.maxFrames 3/handles.maxFrames]);
set(handles.afterSlider,'Value',1-handles.fwdLimit/handles.maxFrames);
set(handles.beforeSlider,'Value',handles.revLimit/handles.maxFrames);
set(handles.spfSlider,'Value',(handles.fps-1)/(handles.maxfps-1));
set(handles.FrameText,'String','frame:0');
set(handles.beforeTxt,'String',num2str(handles.revLimit));
set(handles.afterTxt,'String',num2str(handles.fwdLimit));

sz = handles.halfSize;

for ndx = 1:handles.numSimilar
  tName = sprintf('top%d',ndx);
  mName = sprintf('middle%d',ndx);
  bName = sprintf('bottom%d',ndx);
  handles.hAxes(1,ndx) = initAxes(handles.(tName),sz);
  handles.hAxes(2,ndx) = initAxes(handles.(mName),sz);
  handles.hAxes(3,ndx) = initAxes(handles.(bName),sz);
  
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes showSimilarFrames wait for user response (see UIRESUME)
% uiwait(handles.figure_showSimilarFrames);


% --- Outputs from this function are returned to the command line.
function varargout = showSimilarFrames_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function loadmovies(hObject, movieNames)
% Loads the movies.

handles = guidata(hObject);
pointer = []; trx = {}; labels = {};
for movieNum = 1:length(movieNames)
  fprintf('Loading..');
[pointer(movieNum).readframe,pointer(movieNum).nframes,pointer(movieNum).movie_fid,pointer(movieNum).movieheaderinfo] = ...
  get_readframe_fcn(sprintf('%s/movie.ufmf',movieNames{movieNum}),'interruptible',false);
  trx{movieNum} = load_tracks(sprintf('%s/registered_trx.mat',movieNames{movieNum}));
  fprintf('Done %s\n',movieNames{movieNum});
  labels{movieNum} = load(sprintf('%s/labeledsharpturns.mat',movieNames{movieNum}));
end
handles.movieNames = movieNames;
handles.pointer = pointer;
handles.origTrx = trx;
handles.labels = labels;

guidata(hObject,handles);

function setFrames(hObject,frameData)
% Sets the frames to be shown.

handles = guidata(hObject);

handles.frames{1} = frameData.posFrames;
handles.frames{2} = frameData.curFrame;
handles.frames{3} = frameData.negFrames;
guidata(hObject,handles);

readFrames(hObject);
handles = guidata(hObject);
updatePlots(hObject,handles,handles.frameNo);

function readFrames(hObject)
  handles = guidata(hObject);
  handles.isPlaying = false;
  guidata(hObject,handles);
  curImg = getRotatedFrame(handles,handles.frames{2});
  mtrx = getTracks(handles,handles.frames{2});
  [ptrx ntrx] = getLabels(handles,handles.frames{2});
  for ndx = 1:handles.numSimilar
    handles.im{2,ndx} = curImg;
    handles.trx{2,ndx} = mtrx;
    x = handles.trx{2,ndx}.X;
    y = handles.trx{2,ndx}.Y;
    set(handles.hAxes(2,ndx).labelPos,'XData',x.*ptrx,'Ydata',y.*ptrx);
    set(handles.hAxes(2,ndx).labelNeg,'XData',x.*ntrx,'Ydata',y.*ntrx);
    set(handles.hAxes(2,ndx).trax,'XData',x,'YData',y);
  end
  
  for row = [1 3]
    for ndx = 1:handles.numSimilar
      handles.im{row,ndx} = getRotatedFrame(handles,handles.frames{row}(ndx));
      handles.trx{row,ndx} = getTracks(handles,handles.frames{row}(ndx));
      
      x = handles.trx{row,ndx}.X;
      y = handles.trx{row,ndx}.Y;
      ptrx = handles.trx{row,ndx}.ptrx;
      ntrx = handles.trx{row,ndx}.ntrx;
      set(handles.hAxes(row,ndx).labelPos,'XData',x.*ptrx,'YData',y.*ptrx);    
      set(handles.hAxes(row,ndx).labelNeg,'XData',x.*ntrx,'YData',y.*ntrx);    
      set(handles.hAxes(row,ndx).trax,'XData',x,'YData',y);
      fprintf('.');
    end
    fprintf('\n');
  end
  
  guidata(hObject,handles);
  
function updatePlots(hObject,handles,frameNo)

  for row = 1:3
    for ndx = 1:handles.numSimilar
      set(handles.hAxes(row,ndx).image,'CData',uint8(handles.im{row,ndx}(:,:,:,frameNo)));
      updatefly(handles.hAxes(row,ndx).fly,handles.trx{row,ndx},frameNo);
      set(handles.FrameText,'String',sprintf('frame:%d',frameNo-handles.maxFrames-1));
    end
  end
  
%   guidata(hObject,handles);
  
function play(hObject)
  handles = guidata(hObject);
  frameNo = handles.frameNo;
  ticLoop = tic; % for loop
  ticUpdate = tic; % for updates;
  while(handles.isPlaying)
    if(frameNo-handles.maxFrames>handles.fwdLimit)
      frameNo = handles.maxFrames+1-handles.revLimit;
    else
      frameNo = frameNo+1;   
    end
    guidata(hObject,handles);
    tmp = toc(ticLoop);
    if tmp < 1/handles.fps
      pause(1/handles.fps-tmp);
    end
    ticLoop = tic;
    updatePlots(hObject,handles,frameNo);
    set(handles.frameSlider,'Value', frameNo/(2*handles.maxFrames+1));
    
    if(toc(ticUpdate)>0.5)
      handles = guidata(hObject);
    end
  end
  handles.frameNo = frameNo;
  guidata(hObject,handles);
 
function [ptrx ntrx]= getLabels(handles,curFrame)
  curExp = curFrame.expNum;
  curFly = curFrame.flyNum;
  curTime = curFrame.curTime;
  sz = handles.maxFrames;
  
  tSlice = curTime-sz:curTime+sz;
  ltrx = nan(1,length(tSlice));
  ptrx = nan(1,length(tSlice));
  ntrx = nan(1,length(tSlice));
  
  curLabel = handles.labels{curExp};
  if ~ismember(curFly,curLabel.flies),return; end
  tflyNum = find(curLabel.flies==curFly);
  curT0 = curLabel.t0s{tflyNum};
  curT1 = curLabel.t1s{tflyNum};
  curNames = curLabel.names{tflyNum};
  
  for bnum = 1:length(curNames)
    if strcmpi(curNames{bnum},'None')
      curVal = -1;
    else
      curVal = 1;
    end
    curBoutSlice = curT0(bnum):curT1(bnum);
    overlap = ismember(tSlice,curBoutSlice);
    ltrx(overlap)=curVal;
    
  end
  ptrx(ltrx>0) = 1;
  ntrx(ltrx<0) = 1;
  
  
function trax = getTracks(handles,curFrame)
  curExp = curFrame.expNum;
  curFly = curFrame.flyNum;
  curTime = curFrame.curTime;
  sz = handles.maxFrames;
  
  tSlice = curTime-sz:curTime+sz;
  trax.X = handles.origTrx{curExp}(curFly).x(tSlice);
  trax.Y = handles.origTrx{curExp}(curFly).y(tSlice);
  trax.theta = handles.origTrx{curExp}(curFly).theta(tSlice);
  trax.maj = handles.origTrx{curExp}(curFly).a(tSlice);
  trax.min = handles.origTrx{curExp}(curFly).b(tSlice);

  % Rotate
  curA = handles.origTrx{curExp}(curFly).theta(curTime)-pi/2;
  trax.theta = trax.theta-curA;

  rotMat = [cos(curA) sin(curA) ; -sin(curA) cos(curA) ];
  trax.X = trax.X - trax.X(sz+1);
  trax.Y = trax.Y - trax.Y(sz+1);
  R = rotMat*[trax.X; trax.Y];
  trax.X = R(1,:); trax.Y = R(2,:);
  trax.X = trax.X + handles.halfSize+1;
  trax.Y = trax.Y + handles.halfSize+1;
  [trax.ptrx trax.ntrx] = getLabels(handles,curFrame);
  
function im = getRotatedFrame(handles,curFrame)
% Reads frames around curFrame and rotates them so that the fly is vertical
  sz = handles.halfSize;
  curExp = curFrame.expNum;
  curFly = curFrame.flyNum;
  curTime = curFrame.curTime;
  curX = floor(handles.origTrx{curExp}(curFly).x(curTime));
  curY = floor(handles.origTrx{curExp}(curFly).y(curTime));
  curA = handles.origTrx{curExp}(curFly).theta(curTime);
  tt = zeros(handles.imgY+2*sz,handles.imgX+2*sz,1);
  im = zeros(2*sz+1,2*sz+1,1,2*handles.maxFrames+1);
  bBoxX = (curX-2*sz:curX+2*sz)+sz;
  bBoxY = (curY-2*sz:curY+2*sz)+sz;
  
  for offset = -handles.maxFrames:handles.maxFrames
    tt(sz+1:end-sz,sz+1:end-sz,1) = handles.pointer(curExp).readframe(curTime+offset);
    timg = tt(bBoxY,bBoxX,:);
    rotI = imrotate(timg,curA*180/pi-90,'bilinear','crop');
    im(:,:,:,offset+handles.maxFrames+1) = rotI(sz+1:end-sz,sz+1:end-sz,:);
  end

  
function axesH = initAxes(ax,sz)
  axesH.image = []; 
  hold(ax,'on');
  set(ax,'XLimMode','manual');  xlim(ax,[1 2*sz+1]);
  set(ax,'YLimMode','manual');  ylim(ax,[1 2*sz+1]);
  axesH.image = imagesc(zeros(2*sz+1),'Parent',ax,[0,255]);
  axesH.labelNeg = plot(ax,nan,nan,'Linestyle','-','Color',[0.5 0.1 0.1],'Linewidth',3);
  axesH.labelPos = plot(ax,nan,nan,'Linestyle','-','Color',[0.1 0.5 0.1],'Linewidth',3);
  axesH.trax = plot(ax,nan,nan,'Linestyle','-','Marker','.',...
    'Color',[0.1 0.1 0.1],'MarkerSize',4,'Linewidth',0.1);
  axesH.fly = plot(ax,nan,nan,'Linestyle','-','Color',[0.7 0.2 0.2]);
  colormap(ax,'gray');
  axis(ax,'image','off');

function updatefly(h,trx,t)
% Coped from JCtrax.

% draw an isosceles triangle with center (x,y)
% with rotation theta
% with height maj*4
% with base min*4

x = trx.X(t);
y = trx.Y(t);
theta = trx.theta(t);
maj = trx.maj(t);
min = trx.min(t);

% isosceles triangle not yet rotated or centered
pts = [-maj*2,-min*2,
       -maj*2,min*2,
       maj*2,0];

% rotate
costheta = cos(theta);
sintheta = sin(theta);
R = [costheta,sintheta;-sintheta,costheta];
pts = pts*R;

% translate
pts(:,1) = pts(:,1) + x;
pts(:,2) = pts(:,2) + y;

% plot
set(h,'xdata',pts([1:3,1],1),'ydata',pts([1:3,1],2));
colr = [0.2 0.2 0.2];
if ~isnan(trx.ptrx(t)); colr = [0.1 0.5 0.1]; end
if ~isnan(trx.ntrx(t)); colr = [0.5 0.1 0.1]; end
set(h,'Color',colr);
  
% --- Executes on slider movement.
function frameSlider_Callback(hObject, eventdata, handles)
% hObject    handle to frameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
handles.frameNo = round(v*(2*handles.maxFrames))+1;
guidata(hObject,handles);
updatePlots(hObject,handles,handles.frameNo);



% --- Executes during object creation, after setting all properties.
function frameSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function beforeSlider_Callback(hObject, eventdata, handles)
% hObject    handle to beforeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
handles.revLimit = round(v*handles.maxFrames);
set(handles.beforeTxt,'String',sprintf('%d',handles.revLimit));
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function beforeSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to beforeSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function afterSlider_Callback(hObject, eventdata, handles)
% hObject    handle to afterSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
handles.fwdLimit = round((1-v)*handles.maxFrames);
set(handles.afterTxt,'String',sprintf('%d',handles.fwdLimit));
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function afterSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to afterSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in playButton.
function playButton_Callback(hObject, eventdata, handles)
% hObject    handle to playButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of playButton
if(get(hObject,'Value'));
  handles.isPlaying = true;
  guidata(hObject,handles);
  play(hObject);
else
  handles.isPlaying = false;
  guidata(hObject,handles);  
end


% --- Executes on slider movement.
function spfSlider_Callback(hObject, eventdata, handles)
% hObject    handle to spfSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
v = get(hObject,'Value');
handles.fps = round(v*(handles.maxfps-1))+1;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function spfSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spfSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
