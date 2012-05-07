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

% Last Modified by GUIDE v2.5 15-Mar-2012 16:36:51

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
handles.halfSize = 30*5;
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
handles.display = 'manual';

set(handles.frameSlider,'SliderStep',[1/(2*handles.maxFrames) 3/(2*handles.maxFrames)]);

set(handles.frame_range_edit, 'string', '5');
set(handles.fps_edit, 'string', '5');
set(handles.FrameText,'String','frame:0');

sz = handles.halfSize;
tlsz = handles.maxFrames;

linkaxesGroup = [];

for ndx = 1:handles.numSimilar
  tName = sprintf('top%d',ndx);
  mName = sprintf('middle%d',ndx);
  bName = sprintf('bottom%d',ndx);
  handles.hAxes(1,ndx) = initAxes(handles,handles.(tName),sz,1,ndx);
  handles.hAxes(2,ndx) = initAxes(handles,handles.(mName),sz,2,ndx);
  handles.hAxes(3,ndx) = initAxes(handles,handles.(bName),sz,3,ndx);
  handles.hTimeline(1,ndx) = initTimeline(handles.(['tl_',tName]),tlsz,handles.frameNo);
  handles.hTimeline(2,ndx) = initTimeline(handles.(['tl_',mName]),tlsz,handles.frameNo);
  handles.hTimeline(3,ndx) = initTimeline(handles.(['tl_',bName]),tlsz,handles.frameNo);
  linkaxesGroup = [linkaxesGroup, handles.(tName)];
  linkaxesGroup = [linkaxesGroup, handles.(mName)];
  linkaxesGroup = [linkaxesGroup, handles.(bName)];
end

linkaxes(linkaxesGroup,'xy');
%plot(handles.axes_slider_zero,[0,1],[0,0],'y-','LineWidth',5);
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

function SetJLabelData(hObject,obj,HJLabel)
handles = guidata(hObject);
handles.JLDobj = obj;
handles.JLabelHandles = HJLabel;
handles.imgX = HJLabel.movie_width;
handles.imgY = HJLabel.movie_height;
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
          tStart = max((curBlk-1)*blk+1,1);
          tEnd = min(curBlk*blk,trx(curFly).endframe-trxOffset+1);
          tSlice = tStart:tEnd;
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

function add_prep_list(hObject,exps)
handles = guidata(hObject);
if ~isempty(handles.JLDobj.allperframefns)
  for i = 1:numel(handles.JLDobj.allperframefns),
    handles.timeline_prop_options{i} = handles.JLDobj.allperframefns{i};
  end
  %handles.d = handles.data.allperframefns(1);
  handles.perframepropis = 1;
  set(handles.timeline_properties,'String',handles.timeline_prop_options,'Value',1);
  %handles.timeline_data_ylims = nan(2,numel(handles.data.allperframefns));
end
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

%zoom to the center  
limVal = [handles.halfSize*0.8, handles.halfSize*1.2];
set(handles.middle1,'xLim',limVal,'yLim',limVal);


function readFrames(hObject)
handles = guidata(hObject);
handles.isPlaying = false;
guidata(hObject,handles);
curImg = getRotatedFrame(handles,handles.frames{2});
relTrx = getRelTrx(handles,handles.frames{2});
[labels predictions] = getLabels(handles,handles.frames{2});
ptrx = labels.ptrx;
ntrx = labels.ntrx;

plotPredictionsInAutoMode(handles);
poltLabelsInManualMode(handles);

%for row 2
for t = 1:length(ptrx)
    % tranfer data to image for manual timeline
    timeline{2}(1,t,:) = [0.2 0.2 0.2];
    if ~isnan(ptrx(t)); timeline{2}(1,t,:) = [0.7 0.0 0.0]; end
    if ~isnan(ntrx(t)); timeline{2}(1,t,:) = [0.0 0.0 0.7]; end
    
    % tranfer data to image for automatic timeline
    if isnan(predictions(t)), timeline{2}(2,t,:) = [0.2 0.2 0.2];
    elseif predictions(t)<1.5, timeline{2}(2,t,:) = [0.7 0.0 0.0];
    else timeline{2}(2,t,:) = [0.0 0.0 0.7];
    end
end

for ndx = 1:handles.numSimilar
    handles.im{2,ndx} = curImg;
    handles.trx{2,ndx} = relTrx;
    set(handles.hAxes(2,ndx).labelPos,'visible', 'on');
    set(handles.hAxes(2,ndx).labelNeg,'visible', 'on');
    set(handles.hAxes(2,ndx).predPos,'visible', 'off');
    set(handles.hAxes(2,ndx).predNeg,'visible', 'off');
    set(handles.hAxes(2,ndx).trax,'XData',relTrx.X,'YData',relTrx.Y);
    set(handles.hTimeline(2,ndx).image, 'CData', timeline{2});
    set(handles.hTimeline(2,ndx).hcurr,'XData',[handles.frameNo,handles.frameNo]);
end

%create property plot
%set(handles.axes_timeline_properties,'YLimMode','auto');
hold(handles.axes_timeline_properties,'on');
handles.htimeline_propdata{2} = plot(handles.axes_timeline_properties,nan,nan,'.-');

%choose the first property on the popup list 
s = handles.timeline_prop_options{1};
prop = find(strcmpi(s,handles.JLDobj.allperframefns),1);
handles.perframepropis = prop;
handles.perframeprops = s;

[perframedata,T0,T1] = handles.JLDobj.GetPerFrameData(handles.frames{2}.expNum,handles.frames{2}.flyNum,prop);
TStart = handles.frames{2}.curTime - handles.maxFrames;
IStart = TStart - T0 + 1;
TEnd = handles.frames{2}.curTime + handles.maxFrames;
IEnd = TEnd - T0 + 1;
perframedata = perframedata(IStart:IEnd);
set(handles.htimeline_propdata{2},'XData',-handles.maxFrames:handles.maxFrames,'YData',perframedata,'Color','k', 'LineWidth', 2);

%update text_timeline_prop
handles.cur_rb_row = 2;
handles.cur_rb_col = 1;
perframedata = handles.JLDobj.GetPerFrameData1(handles.frames{2}.expNum,handles.frames{2}.flyNum,prop,handles.frames{2}.curTime);
s = sprintf('%.3f',perframedata);
set(handles.text_timeline_prop,'String',s);

%for row 1 and 3
for row = [1 3]
    for ndx = 1:handles.numSimilar
        relTrx = getRelTrx(handles,handles.frames{row}(ndx));
        handles.im{row,ndx} = getRotatedFrame(handles,handles.frames{row}(ndx));
        
        handles.trx{row,ndx} = relTrx;
        [labels predictions] = getLabels(handles,handles.frames{row}(ndx));
        ptrx = labels.ptrx;
        ntrx = labels.ntrx;
        
        % plot for manual timeline
        for t = 1:length(ptrx)
            timeline{row,ndx}(1,t,:) = [0.2 0.2 0.2];
            if ~isnan(ptrx(t)); timeline{row,ndx}(1,t,:) = [0.7 0.0 0.0]; end
            if ~isnan(ntrx(t)); timeline{row,ndx}(1,t,:) = [0.0 0.0 0.7]; end
            
            % plot for automatic timeline
            if isnan(predictions(t)), timeline{row,ndx}(2,t,:) = [0.2 0.2 0.2];
            elseif predictions(t)<1.5, timeline{row,ndx}(2,t,:) = [0.7 0.0 0.0];
            else timeline{row,ndx}(2,t,:) = [0.0 0.0 0.7];
            end
        end
        
        set(handles.hAxes(row,ndx).labelPos,'visible', 'on');
        set(handles.hAxes(row,ndx).labelNeg,'visible', 'on');
        set(handles.hAxes(row,ndx).predPos,'visible', 'off');
        set(handles.hAxes(row,ndx).predNeg,'visible', 'off');        
        
        set(handles.hAxes(row,ndx).trax,'XData',relTrx.X,'YData',relTrx.Y);
        set(handles.hTimeline(row,ndx).image, 'CData', timeline{row,ndx});
        fprintf('.');
        
        
        %create property plots
        if row == 1
            propColor = 'b';
        else
            propColor = 'r';
        end
        
        
        handles.htimeline_propdata{row}(ndx) = plot(handles.axes_timeline_properties,nan,nan,'.-','color', propColor);
        
        [perframedata,T0,T1] = handles.JLDobj.GetPerFrameData(handles.frames{row}(ndx).expNum,handles.frames{row}(ndx).flyNum,prop);
        TStart = handles.frames{row}(ndx).curTime - handles.maxFrames;
        IStart = TStart - T0 + 1;
        TEnd = handles.frames{row}(ndx).curTime + handles.maxFrames;
        IEnd = TEnd - T0 + 1;
        perframedata = perframedata(IStart:IEnd);
        set(handles.htimeline_propdata{row}(ndx),'XData',-handles.maxFrames:handles.maxFrames,'YData',perframedata);
    end
    fprintf('\n');
end

%plot hcurr for axes_timeline_properties
yLimVal=ylim(handles.axes_timeline_properties);
frameNo = handles.frameNo - handles.maxFrames - 1;
handles.hPropOrig= plot(handles.axes_timeline_properties,[frameNo,frameNo],[yLimVal(1),yLimVal(2)],'g--');
handles.hPropCurr= plot(handles.axes_timeline_properties,[frameNo,frameNo],[yLimVal(1),yLimVal(2)],'k-');
guidata(hObject,handles);
  
function updatePlots(hObject,handles,frameNo)

  for row = 1:3
    for ndx = 1:handles.numSimilar
      set(handles.hAxes(row,ndx).image,'CData',uint8(handles.im{row,ndx}(:,:,:,frameNo)));
      set(handles.hTimeline(row,ndx).hcurr,'XData', [frameNo,frameNo]);
      set(handles.hPropCurr,'XData',[frameNo- handles.maxFrames - 1,frameNo- handles.maxFrames - 1]);
      updatefly(handles.hAxes(row,ndx).fly,handles.hAxes(row,ndx).flyPred,handles.trx{row,ndx},frameNo, handles.display);
    end
    set(handles.FrameText,'String',sprintf('frame:%d',frameNo-handles.maxFrames-1));
  end
  
  %update the text_timeline_prop value
  curExperiment = handles.frames{handles.cur_rb_row}(handles.cur_rb_col).expNum;
  curFly = handles.frames{ handles.cur_rb_row}(handles.cur_rb_col).flyNum;
  curTime = handles.frames{ handles.cur_rb_row}(handles.cur_rb_col).curTime+handles.frameNo-handles.maxFrames-1;
  curProp= handles.perframepropis;
  perframedata = handles.JLDobj.GetPerFrameData1(curExperiment,curFly,curProp,curTime);
  s = sprintf('%.3f',perframedata);
  set(handles.text_timeline_prop,'String',s);

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
    handles.frameNo = frameNo;
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
  curA = double(curTrx.theta(handles.maxFrames+1));
  tt = zeros(handles.imgY+4*sz,handles.imgX+4*sz,1);
  im = zeros(2*sz+1,2*sz+1,1,2*handles.maxFrames+1);
  bBoxX = (curX-2*sz:curX+2*sz)+2*sz;
  bBoxY = (curY-2*sz:curY+2*sz)+2*sz;
  
  pointer = [];
  curMovie = handles.JLDobj.GetFile('movie',curExp);
  [pointer.readframe, pointer.nframes, pointer.movie_fid, pointer.movieheaderinfo] = ...
    get_readframe_fcn(curMovie,'interruptible',false);

  
  for offset = -handles.maxFrames:handles.maxFrames
    if curTime+offset<handles.firstframe{curExp}(curFly) || ...
        curTime+offset>handles.endframe{curExp}(curFly)
      tt(2*sz+1:end-2*sz,2*sz+1:end-2*sz,1) = ...
        zeros(pointer.movieheaderinfo.nr,pointer.movieheaderinfo.nc);
    else
      tt(2*sz+1:end-2*sz,2*sz+1:end-2*sz,1) = ...
        pointer.readframe(curTime+offset);
    end
    timg = tt(bBoxY,bBoxX,:);
    rotI = imrotate(timg,curA*180/pi-90,'bilinear','crop');
    im(:,:,:,offset+handles.maxFrames+1) = rotI(sz+1:end-sz,sz+1:end-sz,:);
  end
  
  if pointer.movie_fid,
    fclose(pointer.movie_fid);
  end

  
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

  
function axesH = initAxes(handles,ax,sz,row,col)
  axesH.image = []; 
  hold(ax,'on');
%   set(ax,'XLimMode','manual');  xlim(ax,[1 2*sz+1]);
%   set(ax,'YLimMode','manual');  ylim(ax,[1 2*sz+1]);

%   set(ax,'XLimMode','manual');  xlim(ax,[1 sz/5*2+1]);
%   set(ax,'YLimMode','manual');  ylim(ax,[1 sz/5*2+1]);
  axesH.image = imagesc(zeros(2*sz+1),'ButtonDownFcn',{@displayOnJLabel,handles,row,col},'Parent',ax,[0,255]);
  axesH.labelNeg = plot(ax,nan,nan,'Linestyle','-','Color',[0.0 0.0 0.7],'Linewidth',3);
  axesH.labelPos = plot(ax,nan,nan,'Linestyle','-','Color',[0.7 0.0 0.0],'Linewidth',3);
  axesH.predNeg = plot(ax,nan,nan,'Linestyle','-','Color',[0.0 0.0 0.7],'Linewidth',3);
  axesH.predPos = plot(ax,nan,nan,'Linestyle','-','Color',[0.7 0.0 0.0],'Linewidth',3);
  axesH.trax = plot(ax,nan,nan,'Linestyle','-','Marker','.',...
    'Color',[0.1 0.1 0.1],'MarkerSize',4,'Linewidth',0.1);
  axesH.fly = plot(ax,nan,nan,'Linestyle','-','Color',[0.7 0.2 0.2]);
  axesH.flyPred = plot(ax,nan,nan,'Linestyle','-','Color',[0.7 0.2 0.2],...
    'Linewidth',2);
  colormap(ax,'gray');
  axis(ax,'off','image');

function tlH = initTimeline(tl,sz,tcurr)
  tlH.image = []; 
  hold(tl,'on');
  %set(tl,'XLimMode','auto');  xlim(tl,[1 2*sz+1]);
  %set(tl,'YLimMode','auto');  ylim(tl,[1 2*sz+1]);
  tlH.image = imagesc(zeros(2,2*sz+1,3),'Parent',tl,[0,255]);
  tlH.horig = plot(tl,[tcurr,tcurr],[0, 3],'g--','HitTest','off','linewidth',2);
  tlH.hcurr = plot(tl,[tcurr,tcurr],[0, 3],'y-','HitTest','off','linewidth',2);
  colormap(tl,'gray');
  axis(tl,'off','image');
  
function updatefly(h1,h2,trx,t,displayChoice)
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

switch displayChoice
    case 'automatic'
        if isnan(trx.predictions(t)), colr = [0.2 0.2 0.2];
        elseif trx.predictions(t)<1.5, colr = [0.7 0.0 0.0];
        else colr = [0.0 0.0 0.7];
        end
    case 'manual'
        colr = [0.2 0.2 0.2];
        if ~isnan(trx.labels.ptrx(t)); colr = [0.7 0.0 0.0]; end
        if ~isnan(trx.labels.ntrx(t)); colr = [0.0 0.0 0.7]; end
    otherwise
        colr = [0.2 0.2 0.2];
end
    
% plot for labels.
set(h1,'xdata',pts([1 2],1),'ydata',pts([1 2],2));
set(h1,'Color',colr);
set(h2,'xdata',pts([2 3 1], 1),'ydata',pts([2 3 1],2));
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


% --- Executes when user attempts to close figure_showSimilarFrames.
function figure_showSimilarFrames_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure_showSimilarFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure

handles.JLDobj.frameFig = [];
delete(hObject);




function fps_edit_Callback(hObject, eventdata, handles)
% hObject    handle to fps_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of fps_edit as text
%        str2double(get(hObject,'String')) returns contents of fps_edit as a double
v = str2num(get(hObject,'String'));
if (isempty(v))||(v >handles.maxfps|| v < 0)
    error('frame range requires 1 argument that is a number between 0 and maxfps ');
end

handles.fps = round(v);
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function fps_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fps_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function frame_range_edit_Callback(hObject, eventdata, handles)
% hObject    handle to frame_range_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of frame_range_edit as text
%        str2double(get(hObject,'String')) returns contents of frame_range_edit as a double

frameRanges = str2double(get(hObject,'String'));
if (isempty(frameRanges))||(frameRanges >handles.maxFrames||frameRanges < 0)
    warnStr = sprintf('Frame range requires 1 argument that is a number between 0 and %d, please reinput a value.\n',handles.maxFrames);
    warndlg(warnStr,'frame range is out of range');
else
    
    handles.revLimit = round(frameRanges);
    handles.fwdLimit = handles.revLimit;
    
    guidata(hObject,handles);
end


% --- Executes during object creation, after setting all properties.
function frame_range_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frame_range_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in timeline_properties.
function timeline_properties_Callback(hObject, eventdata, handles)
% hObject    handle to timeline_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns timeline_properties contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        timeline_properties

v = get(hObject,'Value');
s = handles.timeline_prop_options{v};

prop = find(strcmpi(s,handles.JLDobj.allperframefns),1);
handles.perframepropis = prop;
handles.perframeprops = s;

[perframedata,T0,T1] = handles.JLDobj.GetPerFrameData(handles.frames{2}.expNum,handles.frames{2}.flyNum,prop);
TStart = handles.frames{2}.curTime - handles.maxFrames;
IStart = TStart - T0 + 1;
TEnd = handles.frames{2}.curTime + handles.maxFrames;
IEnd = TEnd - T0 + 1;
perframedata = perframedata(IStart:IEnd);
set(handles.htimeline_propdata{2},'XData',-handles.maxFrames:handles.maxFrames,'YData',perframedata);

%reset Ymax and Ymin value for current time and original time
set(handles.hPropOrig,'YData',[0,0]);
set(handles.hPropCurr,'XData',[0,0],'YData',[0,0]);

%for row 1 and 3
for row = [1 3]
    for ndx = 1:handles.numSimilar
        [perframedata,T0,T1] = handles.JLDobj.GetPerFrameData(handles.frames{row}(ndx).expNum,handles.frames{row}(ndx).flyNum,prop);
        TStart = handles.frames{row}(ndx).curTime - handles.maxFrames;
        IStart = TStart - T0 + 1;
        TEnd = handles.frames{row}(ndx).curTime + handles.maxFrames;
        IEnd = TEnd - T0 + 1;
        perframedata = perframedata(IStart:IEnd);
        set(handles.htimeline_propdata{row}(ndx),'XData',-handles.maxFrames:handles.maxFrames,'YData',perframedata);
    end
    fprintf('\n');
end

%plot hcurr for axes_timeline_properties
yLimVal=ylim(handles.axes_timeline_properties);
frameNo = handles.frameNo - handles.maxFrames - 1;
set(handles.hPropOrig,'YData',[yLimVal(1),yLimVal(2)]);
set(handles.hPropCurr,'XData',[frameNo,frameNo],'YData',[yLimVal(1),yLimVal(2)]);

%update the text_timeline_prop value
curExperiment = handles.frames{handles.cur_rb_row}(handles.cur_rb_col).expNum;
curFly = handles.frames{ handles.cur_rb_row}(handles.cur_rb_col).flyNum; 
curTime = handles.frames{ handles.cur_rb_row}(handles.cur_rb_col).curTime+handles.frameNo-handles.maxFrames-1;
curProp= handles.perframepropis;
perframedata = handles.JLDobj.GetPerFrameData1(curExperiment,curFly,curProp,curTime);
s = sprintf('%.3f',perframedata);
set(handles.text_timeline_prop,'String',s);

guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function timeline_properties_CreateFcn(hObject, eventdata, handles)
% hObject    handle to timeline_properties (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in uipanel6.
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel6 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'rb_automatic'
        for row = 1:3
            for ndx = 1:handles.numSimilar
                
                % Code for when rb_top1 is selected.
                handles.display = 'automatic';                
                set(handles.hAxes(row,ndx).predPos,'visible', 'on');
                set(handles.hAxes(row,ndx).predNeg,'visible', 'on');
                set(handles.hAxes(row,ndx).labelPos,'visible', 'off');
                set(handles.hAxes(row,ndx).labelNeg,'visible', 'off');

            end
        end
        
    case 'rb_manual'
        for row = 1:3
            for ndx = 1:handles.numSimilar
                handles.display = 'manual';
                set(handles.hAxes(row,ndx).labelPos,'visible', 'on');
                set(handles.hAxes(row,ndx).labelNeg,'visible', 'on');
                set(handles.hAxes(row,ndx).predPos,'visible', 'off');
                set(handles.hAxes(row,ndx).predNeg,'visible', 'off');
            end
        end
    otherwise
        % Code for when there is no match.
end
updatePlots(hObject,handles,handles.frameNo);
guidata(hObject,handles);


function displayOnJLabel(obj,eventdata,handles,row,col)

switch get(gcf,'selectiontype')
    case 'normal'
        % LEFT CLICK
    case 'extend'
        % SHIFT-CLICK LEFT (or L/R simultaneous)
    case 'alt'
        % CTRL-CLICK LEFT (or RIGHT-CLICK)
    case 'open'
        % DOUBLE-CLICK
        
        handles=guidata(obj);
        JLabelHandles = handles.JLabelHandles;
        figure(JLabelHandles.figure_JLabel);
        if row==2
            JLabelHandles.guidata.ts(1) =999; % in order to update plot in SetCurrentFrame
            handles.curExp = handles.frames{2}.expNum;
            [JLabelHandles,~] = JLabel('SetCurrentMovie',JLabelHandles,handles.curExp);
            handles.curFly = handles.frames{2}.flyNum;
            JLabelHandles = JLabel('SetCurrentFlies',JLabelHandles,handles.curFly);
            handles.curFrame = handles.frames{2}.curTime;
            JLabelHandles = JLabel('SetCurrentFrame',JLabelHandles,1,handles.curFrame,obj);
            guidata(obj, handles);
        else
            handles.curExp = handles.frames{row}(col).expNum;
            [JLabelHandles,~] = JLabel('SetCurrentMovie',JLabelHandles,handles.curExp);
            handles.curFly = handles.frames{row}(col).flyNum;
            JLabelHandles = JLabel('SetCurrentFlies',JLabelHandles,handles.curFly);
            handles.curFrame = handles.frames{row}(col).curTime;
            JLabelHandles = JLabel('SetCurrentFrame',JLabelHandles,1,handles.curFrame,obj);
            guidata(obj, handles);
        end
        
end

% --- Executes when selected object is changed in chooseFrame.
function chooseFrame_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in chooseFrame 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

old_rb_row = handles.cur_rb_row;
old_rb_col = handles.cur_rb_col;

switch get(eventdata.NewValue,'Tag') % Get Tag of selected object.
    case 'rb_top1'
        % Code for when rb_top1 is selected.
        handles.cur_rb_row = 1;
        handles.cur_rb_col = 1;
    case 'rb_top2'
        % Code for when radiobutton2 is selected.
        handles.cur_rb_row = 1;
        handles.cur_rb_col = 2;
    case 'rb_top3'
        handles.cur_rb_row = 1;
        handles.cur_rb_col = 3;
        
    case 'rb_top4'
        handles.cur_rb_row = 1;
        handles.cur_rb_col = 4;
        
    case 'rb_middle1'
        handles.cur_rb_row = 2;
        handles.cur_rb_col = 1;
    case 'rb_bottom1'
        handles.cur_rb_row = 3;
        handles.cur_rb_col = 1;
    case 'rb_bottom2'
        handles.cur_rb_row = 3;
        handles.cur_rb_col = 2;
        
    case 'rb_bottom3'
        handles.cur_rb_row = 3;
        handles.cur_rb_col = 3;
    case 'rb_bottom4'
        handles.cur_rb_row = 3;
        handles.cur_rb_col = 4;
        
    otherwise
        % Code for when there is no match.
end

%change the color of the old lines to blue
set(handles.htimeline_propdata{old_rb_row}(old_rb_col),'LineWidth',0.5);

%change the color of the new line to red
set(handles.htimeline_propdata{handles.cur_rb_row}(handles.cur_rb_col),'LineWidth',2);

%update the text_timeline_prop value
curExperiment = handles.frames{handles.cur_rb_row}(handles.cur_rb_col).expNum;
curFly = handles.frames{ handles.cur_rb_row}(handles.cur_rb_col).flyNum; 
curTime = handles.frames{ handles.cur_rb_row}(handles.cur_rb_col).curTime+handles.frameNo-handles.maxFrames-1;
curProp= handles.perframepropis;
perframedata = handles.JLDobj.GetPerFrameData1(curExperiment,curFly,curProp,curTime);
s = sprintf('%.3f',perframedata);
set(handles.text_timeline_prop,'String',s);

guidata(hObject,handles);

function plotPredictionsInAutoMode(handles)
%for row 2

relTrx = getRelTrx(handles,handles.frames{2});
[labels predictions] = getLabels(handles,handles.frames{2});
ptrx = nan(size(predictions));
ntrx = nan(size(predictions));

for i=1:length(ptrx)
    if predictions(i)<1.5
        ptrx(i) = 1;
    else
        ntrx(i) = 1;
    end
end

%To fill in the gap between the lablePos and labelNeg
diffPred = diff(predictions);
changeIdex = find(diffPred);

for i = 1:length(changeIdex)
    if diffPred(changeIdex(i))<0 %negative -> positive
        ntrx(changeIdex(i)+1) = 1;
    else %positive -> negative
        ptrx(changeIdex(i)+1) = 1;
    end
end

for ndx = 1:handles.numSimilar
    handles.trx{2,ndx} = relTrx;
    set(handles.hAxes(2,ndx).predPos,'XData',relTrx.X.*ptrx,'Ydata',relTrx.Y.*ptrx);
    set(handles.hAxes(2,ndx).predNeg,'XData',relTrx.X.*ntrx,'Ydata',relTrx.Y.*ntrx);
end

%for row 1 and 3
for row = [1 3]
    for ndx = 1:handles.numSimilar
        relTrx = getRelTrx(handles,handles.frames{row}(ndx));
        handles.trx{row,ndx} = relTrx;
        [labels predictions] = getLabels(handles,handles.frames{row}(ndx));
        
        ptrx = nan(size(predictions));
        ntrx = nan(size(predictions));
        
        for i=1:length(ptrx)
            if predictions(i)<1.5
                ptrx(i) = 1;
            else
                ntrx(i) = 1;
            end
        end
        
        diffPred = diff(predictions);
        changeIdex = find(diffPred);
        
        for i = 1:length(changeIdex)
            if diffPred(changeIdex(i))<0 %negative -> positive
                ntrx(changeIdex(i)+1) = 1;
            else %positive -> negative
                ptrx(changeIdex(i)+1) = 1;
            end
        end
        
        set(handles.hAxes(row,ndx).predPos,'XData',relTrx.X.*ptrx,'YData',relTrx.Y.*ptrx);
        set(handles.hAxes(row,ndx).predNeg,'XData',relTrx.X.*ntrx,'YData',relTrx.Y.*ntrx);
    end
end


function poltLabelsInManualMode(handles)

%for row 2
relTrx = getRelTrx(handles,handles.frames{2});
[labels predictions] = getLabels(handles,handles.frames{2});
ptrx = labels.ptrx;
ntrx = labels.ntrx;

for ndx = 1:handles.numSimilar
    handles.trx{2,ndx} = relTrx;
    set(handles.hAxes(2,ndx).labelPos,'XData',relTrx.X.*ptrx,'Ydata',relTrx.Y.*ptrx);
    set(handles.hAxes(2,ndx).labelNeg,'XData',relTrx.X.*ntrx,'Ydata',relTrx.Y.*ntrx);
end

%for row 1 and 3
for row = [1 3]
    for ndx = 1:handles.numSimilar
        relTrx = getRelTrx(handles,handles.frames{row}(ndx));
        handles.trx{row,ndx} = relTrx;
        [labels predictions] = getLabels(handles,handles.frames{row}(ndx));
        ptrx = labels.ptrx;
        ntrx = labels.ntrx;
        
        set(handles.hAxes(row,ndx).labelPos,'XData',relTrx.X.*ptrx,'YData',relTrx.Y.*ptrx);
        set(handles.hAxes(row,ndx).labelNeg,'XData',relTrx.X.*ntrx,'YData',relTrx.Y.*ntrx);
    end
end



% --- Executes on key press with focus on figure_showSimilarFrames and none of its controls.
function figure_showSimilarFrames_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure_showSimilarFrames (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)
switch eventdata.Key,
  
  case 'leftarrow',
      go_previous_frame(hObject, handles);
     
  case 'rightarrow',
      go_next_frame(hObject, handles);
      
end

function go_next_frame(hObject, handles)
% hObject    handle to menu_go_next_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TODO: make this work with multiple preview axes

frameNo = handles.frameNo;
if(frameNo-handles.maxFrames>handles.fwdLimit)
    frameNo = handles.maxFrames+1-handles.revLimit;
else
    frameNo = frameNo+1;
end

handles.frameNo = frameNo;

updatePlots(hObject,handles,frameNo);
set(handles.frameSlider,'Value', frameNo/(2*handles.maxFrames+1));

guidata(hObject,handles);

% --------------------------------------------------------------------
function go_previous_frame(hObject, handles)
% hObject    handle to menu_go_previous_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

frameNo = handles.frameNo;
if(handles.maxFrames - frameNo>= handles.revLimit - 1)
    frameNo = handles.maxFrames + handles.fwdLimit + 1;
else
    frameNo = frameNo-1;
end

handles.frameNo = frameNo;

updatePlots(hObject,handles,frameNo);
set(handles.frameSlider,'Value', frameNo/(2*handles.maxFrames+1));

guidata(hObject,handles);


% --- Executes on mouse press over figure background.
function figure_showSimilarFrames_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure_showSimilarFrames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function zoomOut_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to zoomOut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldLimVal = get(handles.middle1,'xLim');
mag = (oldLimVal(2) - oldLimVal(1))/(handles.halfSize*0.4);

if mag<=4
    mag = mag +1;
else
    mag=5;
end
limVal = [handles.halfSize*(1-mag/5), handles.halfSize*(1+mag/5)];
set(handles.middle1,'xLim',limVal,'yLim',limVal);


% --------------------------------------------------------------------
function zoomIn_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to zoomIn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
oldLimVal = get(handles.middle1,'xLim');
mag = (oldLimVal(2) - oldLimVal(1))/(handles.halfSize*0.4);

if mag>=2
    mag = mag-1;
else
    mag=1;
end
limVal = [handles.halfSize*(1-mag/5), handles.halfSize*(1+mag/5)];
set(handles.middle1,'xLim',limVal,'yLim',limVal);
