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
handles.tsize = 500;
handles.cache = containers.Map('keyType','char','valueType','any');

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

function SetJLabelData(hObject,obj)
handles = guidata(hObject);
handles.JLDobj = obj;
guidata(hObject,handles);

%{
function loadmovies(hObject)
% Loads the movies.

handles = guidata(hObject);
pointer = []; 
movieNames = handles.JLDobj.expdirs;
for movieNum = 1:length(movieNames)
[pointer(movieNum).readframe,pointer(movieNum).nframes,pointer(movieNum).movie_fid,pointer(movieNum).movieheaderinfo] = ...
  get_readframe_fcn(sprintf('%s/movie.ufmf',movieNames{movieNum}),'interruptible',false);
end
handles.movieNames = movieNames;
handles.pointer = pointer;

guidata(hObject,handles);
%}

function CacheTracksLabeled(hObject,exps)
% Load the tracks.

handles = guidata(hObject);
blk = handles.tsize;
off = handles.maxFrames;
fprintf('Loading ');
if nargin<2
  exps = 1:handles.JLDobj.nexps;
end
for expNum = exps
  fprintf('.');
  trxfile = handles.JLDobj.GetFile('trx',expNum);
  trx = load_tracks(trxfile);
  
  handles.firstframe{expNum} = [trx.firstframe];
  handles.endframe{expNum} = [trx.endframe];

  % Get labels. Store only those parts of tracks that have been labeled.
  for curFly = 1:length(trx)
    labels = handles.JLDobj.GetLabels(expNum,curFly);
    for curNdx = 1:length(labels.t0s)
      trxOffset = handles.firstframe{expNum}(curFly);
      curblockStart = floor( (labels.t0s(curNdx)-off-trxOffset+1)/blk)+1;
      curblockEnd = ceil( (labels.t1s(curNdx)+off-trxOffset+1)/blk);
      for curBlk = curblockStart:curblockEnd
        idStr = sprintf('%d_%d_%d',expNum,curFly,curBlk);
        if ~isKey(handles.cache,idStr),
          tSlice = ((curBlk-1)*blk+1):(curBlk*blk);
          trax = [];
          trax.X = trx(curFly).x(tSlice);
          trax.Y = trx(curFly).y(tSlice);
          trax.theta = trx(curFly).theta(tSlice);
          trax.maj = trx(curFly).a(tSlice);
          trax.min = trx(curFly).b(tSlice);
          handles.cache(idStr) = trax;
        end
      end
    end
    
  end
end
fprintf(' Done loading tracks\n');

guidata(hObject,handles);


function CacheTracks(hObject,expNum,curFly,t0,t1)
handles = guidata(hObject);
blk = handles.tsize;
off = handles.maxFrames;
trxfile = handles.JLDobj.GetFile('trx',expNum);
trx = load_tracks(trxfile);

handles.firstframe{expNum} = [trx.firstframe];
handles.endframe{expNum} = [trx.endframe];
trxOffset = handles.firstframe{expNum}(curFly);
curblockStart = floor( (t0-off-trxOffset+1)/blk)+1;
curblockEnd = ceil( (t1+off-trxOffset+1)/blk);

for curBlk = curblockStart:curblockEnd
  idStr = sprintf('%d_%d_%d',expNum,curFly,curBlk);
  if ~isKey(handles.cache,idStr),
    tSlice = ((curBlk-1)*blk+1):(curBlk*blk);
    tSliceValid = tSlice >= 1 & tSlice<=(trx(curFly).endframe-trxOffset+1);
    trax = [];
    trax.X = nan(1,length(tSlice));
    trax.Y = nan(1,length(tSlice));
    trax.theta = nan(1,length(tSlice));
    trax.maj = nan(1,length(tSlice));
    trax.min = nan(1,length(tSlice));
    
    validNdx = tSlice(tSliceValid);
    trax.X(tSliceValid) = trx(curFly).x(validNdx);
    trax.Y(tSliceValid) = trx(curFly).y(validNdx);
    trax.theta(tSliceValid) = trx(curFly).theta(validNdx);
    trax.maj(tSliceValid) = trx(curFly).a(validNdx);
    trax.min(tSliceValid) = trx(curFly).b(validNdx);
    handles.cache(idStr) = trax;
  end
end
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
  relTrx = getRelTrx(handles,handles.frames{2});
  [labels predictions] = getLabels(handles,handles.frames{2});
  ptrx = labels.ptrx;
  ntrx = labels.ntrx;
  for ndx = 1:handles.numSimilar
    handles.im{2,ndx} = curImg;
    handles.trx{2,ndx} = relTrx;
    set(handles.hAxes(2,ndx).labelPos,'XData',relTrx.X.*ptrx,'Ydata',relTrx.Y.*ptrx);
    set(handles.hAxes(2,ndx).labelNeg,'XData',relTrx.X.*ntrx,'Ydata',relTrx.Y.*ntrx);
    set(handles.hAxes(2,ndx).trax,'XData',relTrx.X,'YData',relTrx.Y);
  end
  
  for row = [1 3]
    for ndx = 1:handles.numSimilar
      relTrx = getRelTrx(handles,handles.frames{row}(ndx));
      handles.im{row,ndx} = getRotatedFrame(handles,handles.frames{row}(ndx));
      
      handles.trx{row,ndx} = relTrx;
      [labels predictions] = getLabels(handles,handles.frames{row}(ndx));
      ptrx = labels.ptrx;
      ntrx = labels.ntrx;
      set(handles.hAxes(row,ndx).labelPos,'XData',relTrx.X.*ptrx,'YData',relTrx.Y.*ptrx);    
      set(handles.hAxes(row,ndx).labelNeg,'XData',relTrx.X.*ntrx,'YData',relTrx.Y.*ntrx);    
      set(handles.hAxes(row,ndx).trax,'XData',relTrx.X,'YData',relTrx.Y);
      fprintf('.');
    end
    fprintf('\n');
  end
  
  guidata(hObject,handles);
  
function updatePlots(hObject,handles,frameNo)

  for row = 1:3
    for ndx = 1:handles.numSimilar
      set(handles.hAxes(row,ndx).image,'CData',uint8(handles.im{row,ndx}(:,:,:,frameNo)));
      updatefly(handles.hAxes(row,ndx).fly,handles.hAxes(row,ndx).flyPred,handles.trx{row,ndx},frameNo);
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
 
  
function [labels predictions]= getLabels(handles,curFrame)
% Reads labels and predictions.

  curExp = curFrame.expNum;
  curFly = curFrame.flyNum;
  curTime = curFrame.curTime;
  sz = handles.maxFrames;
  
  tSlice = curTime-sz:curTime+sz;
  ltrx = nan(1,length(tSlice));
  labels.ptrx = nan(1,length(tSlice));
  labels.ntrx = nan(1,length(tSlice));
  
  curLabel = handles.JLDobj.GetLabels(curExp,curFly);
%   windowNdx = find( (handles.JLDobj.windowdata.exp == curExp) & ...
%       (handles.JLDobj.windowdata.flies == curFly) & ...
%       (handles.JLDobj.windowdata.t == curTime) ,1);
    
  idxcurr = handles.JLDobj.windowdata.exp == curExp & ...
    all(bsxfun(@eq,handles.JLDobj.windowdata.flies,curFly),2) & ...
    handles.JLDobj.windowdata.t >= (curTime-sz) & ...
    handles.JLDobj.windowdata.t <= (curTime+sz) & ...
    handles.JLDobj.windowdata.isvalidprediction;

  tNdx = handles.JLDobj.windowdata.t(idxcurr);  
  predictions = nan(1,length(tSlice));
  predictions(  (tNdx-curTime+sz+1) ) = handles.JLDobj.windowdata.predicted(idxcurr);
  curT0 = curLabel.t0s;
  curT1 = curLabel.t1s;
  curNames = curLabel.names;
  
  for bnum = 1:length(curNames)
    if strcmpi(curNames{bnum},'None')
      curVal = -1;
    else
      curVal = 1;
    end
    curBoutSlice = curT0(bnum):(curT1(bnum)-1);
    overlap = ismember(tSlice,curBoutSlice);
    ltrx(overlap)=curVal;
    
  end
  labels.ptrx(ltrx>0) = 1;
  labels.ntrx(ltrx<0) = 1;

  
function relTrax = getRelTrx(handles,curFrame)
  curExp = curFrame.expNum;
  curFly = curFrame.flyNum;
  curTime = curFrame.curTime;
  sz = handles.maxFrames;
  
  tSlice = curTime-sz:curTime+sz;
  trax = readCache(handles,curExp,curFly,curTime);
  
  % Rotate
  curA = trax.theta(handles.maxFrames+1)-pi/2;
  trax.theta = trax.theta-curA;

  rotMat = [cos(curA) sin(curA) ; -sin(curA) cos(curA) ];
  trax.X = trax.X - trax.X(sz+1);
  trax.Y = trax.Y - trax.Y(sz+1);
  R = rotMat*[trax.X; trax.Y];
  trax.X = R(1,:); trax.Y = R(2,:);
  trax.X = trax.X + handles.halfSize+1;
  trax.Y = trax.Y + handles.halfSize+1;
  [trax.labels trax.predictions] = getLabels(handles,curFrame);
  relTrax = trax;

  
function im = getRotatedFrame(handles,curFrame)
% Reads frames around curFrame and rotates them so that the fly is vertical
  sz = handles.halfSize;
  curExp = curFrame.expNum;
  curFly = curFrame.flyNum;
  curTime = curFrame.curTime;
  curTrx = readCache(handles,curExp,curFly,curTime);
  curX = round(curTrx.X(handles.maxFrames+1));
  curY = round(curTrx.Y(handles.maxFrames+1));
  curA = curTrx.theta(handles.maxFrames+1);
  tt = zeros(handles.imgY+2*sz,handles.imgX+2*sz,1);
  im = zeros(2*sz+1,2*sz+1,1,2*handles.maxFrames+1);
  bBoxX = (curX-2*sz:curX+2*sz)+sz;
  bBoxY = (curY-2*sz:curY+2*sz)+sz;
  
  pointer = [];
  curMovie = handles.JLDobj.GetFile('movie',curExp);
  [pointer.readframe, pointer.nframes, pointer.movie_fid, pointer.movieheaderinfo] = ...
    get_readframe_fcn(curMovie,'interruptible',false);

  
  for offset = -handles.maxFrames:handles.maxFrames
    if curTime+offset<handles.firstframe{curExp}(curFly) || ...
        curTime+offset>handles.endframe{curExp}(curFly)
      tt(sz+1:end-sz,sz+1:end-sz,1) = ...
        zeros(pointer.movieheaderinfo.nr,pointer.movieheaderinfo.nc);
    else
      tt(sz+1:end-sz,sz+1:end-sz,1) = ...
        pointer.readframe(curTime+offset);
    end
    timg = tt(bBoxY,bBoxX,:);
    rotI = imrotate(timg,curA*180/pi-90,'bilinear','crop');
    im(:,:,:,offset+handles.maxFrames+1) = rotI(sz+1:end-sz,sz+1:end-sz,:);
  end
  
  fclose(pointer.movie_fid);

  
function trx = readCache(handles,expNum,flyNum,curT)

  blk = handles.tsize;
  sz = handles.maxFrames;
  trxOffset = handles.firstframe{expNum}(flyNum);
  blkNumStart = floor( (curT-sz-trxOffset+1)/blk)+1;
  blkNumEnd = ceil( (curT+sz-trxOffset+1)/blk);
  off = curT-trxOffset+1 - (blkNumStart-1)*blk - sz;
  trx.X = []; trx.Y = []; trx.theta = []; trx.maj = []; trx.min = [];
  for curBlk = blkNumStart:blkNumEnd
    idStr = sprintf('%d_%d_%d',expNum,flyNum,curBlk);
    if ~isKey(handles.cache,idStr)
      CacheTracks(handles.output,expNum,flyNum,curT,curT);
      handles = guidata(handles.output);
    end
    curTrx = handles.cache(idStr);
    trx.X = [trx.X curTrx.X];
    trx.Y = [trx.Y curTrx.Y];
    trx.theta = [trx.theta curTrx.theta];
    trx.maj = [trx.maj curTrx.maj];
    trx.min = [trx.min curTrx.min];
  end
 
  trx.X = trx.X(off:(off+2*sz));
  trx.Y = trx.Y(off:(off+2*sz));
  trx.theta = trx.theta(off:(off+2*sz));
  trx.maj = trx.maj(off:(off+2*sz));
  trx.min = trx.min(off:(off+2*sz));

  
function axesH = initAxes(ax,sz)
  axesH.image = []; 
  hold(ax,'on');
  set(ax,'XLimMode','manual');  xlim(ax,[1 2*sz+1]);
  set(ax,'YLimMode','manual');  ylim(ax,[1 2*sz+1]);
  axesH.image = imagesc(zeros(2*sz+1),'Parent',ax,[0,255]);
  axesH.labelNeg = plot(ax,nan,nan,'Linestyle','-','Color',[0.0 0.0 0.7],'Linewidth',3);
  axesH.labelPos = plot(ax,nan,nan,'Linestyle','-','Color',[0.7 0.0 0.0],'Linewidth',3);
  axesH.trax = plot(ax,nan,nan,'Linestyle','-','Marker','.',...
    'Color',[0.1 0.1 0.1],'MarkerSize',4,'Linewidth',0.1);
  axesH.fly = plot(ax,nan,nan,'Linestyle','-','Color',[0.7 0.2 0.2]);
  axesH.flyPred = plot(ax,nan,nan,'Linestyle','-','Color',[0.7 0.2 0.2],...
    'Linewidth',2);
  colormap(ax,'gray');
  axis(ax,'image','off');

  
function updatefly(h1,h2,trx,t)
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

% plot for labels.
set(h1,'xdata',pts([1 2],1),'ydata',pts([1 2],2));
colr = [0.2 0.2 0.2];
if ~isnan(trx.labels.ptrx(t)); colr = [0.7 0.0 0.0]; end
if ~isnan(trx.labels.ntrx(t)); colr = [0.0 0.0 0.7]; end
set(h1,'Color',colr);

% plot for labels
set(h2,'xdata',pts([2 3 1], 1),'ydata',pts([2 3 1],2));

colr = [0.2 0.2 0.2];
if isnan(trx.predictions(t)), colr = [0.2 0.2 0.2];
elseif trx.predictions(t)<1.5, colr = [0.7 0.0 0.0]; 
else colr = [0.0 0.0 0.7]; 
end

set(h2,'Color',colr);


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
