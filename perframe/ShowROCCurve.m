function varargout = ShowROCCurve(varargin)
% SHOWROCCURVE MATLAB code for ShowROCCurve.fig
%      SHOWROCCURVE, by itself, creates a new SHOWROCCURVE or raises the existing
%      singleton*.
%
%      H = SHOWROCCURVE returns the handle to a new SHOWROCCURVE or the handle to
%      the existing singleton*.
%
%      SHOWROCCURVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SHOWROCCURVE.M with the given input arguments.
%
%      SHOWROCCURVE('Property','Value',...) creates a new SHOWROCCURVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ShowROCCurve_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ShowROCCurve_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ShowROCCurve

% Last Modified by GUIDE v2.5 25-Jan-2012 10:03:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ShowROCCurve_OpeningFcn, ...
                   'gui_OutputFcn',  @ShowROCCurve_OutputFcn, ...
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

% --- Executes just before ShowROCCurve is made visible.
function ShowROCCurve_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ShowROCCurve (see VARARGIN)

% Choose default command line output for ShowROCCurve
handles.output = hObject;
handles.labels = varargin{1};
handles.scores = varargin{2};
handles.JLDObj = varargin{3};

numBins = 21;
handles.scoreNorm = handles.JLDObj.windowdata.scoreNorm;
pos = handles.labels==1;
neg = ~pos;
bins = linspace(-handles.scoreNorm,handles.scoreNorm,numBins);
bins = [-inf bins(2:end-1) inf];
histPos = histc(handles.scores(pos),bins);
histNeg = histc(handles.scores(neg),bins);
histPos(end) = []; histNeg(end) = [];

handles.thres1 = round(handles.JLDObj.GetConfidenceThreshold(1)*numBins/2)/(numBins/2);
handles.thres2 = -round(handles.JLDObj.GetConfidenceThreshold(2)*numBins/2)/(numBins/2);


xLocs = linspace(-1+1/(numBins-1),1-1/(numBins-1),numBins-1);
hBar = bar(handles.axes1, xLocs,[histPos histNeg],'BarWidth',1.5);
ylim = get(handles.axes1,'ylim');
set(handles.axes1,'xlim',[-1 1])
shortylim = ylim;
shortylim(1) = ylim(1) + 0.001;
shortylim(2) = ylim(2) - 0.001;
linepos = imline(handles.axes1,[handles.thres1 handles.thres1],shortylim);
lineneg = imline(handles.axes1,[handles.thres2 handles.thres2],shortylim);
setColor(linepos,[0 0 1]);
setColor(lineneg,[1 0 0]);

possibleXLocs = linspace(-1,1,numBins);
possibleXLocs = possibleXLocs(2:end-1);
api = iptgetapi(linepos);
fcn = GetLineConstraintFcn(get(handles.axes1,'XLim'),...
   get(handles.axes1,'YLim'),[[handles.thres1 handles.thres1]; ylim],...
   possibleXLocs);
api.setPositionConstraintFcn(fcn);

api = iptgetapi(lineneg);
fcn = GetLineConstraintFcn(get(handles.axes1,'XLim'),...
   get(handles.axes1,'YLim'),[[handles.thres2 handles.thres2]; ylim],...
   possibleXLocs);
api.setPositionConstraintFcn(fcn);
% UpdateText(handles);
% % Update handles structure
guidata(hObject, handles);

% UIWAIT makes ShowROCCurve wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function UpdateText(handles)
pos = handles.labels>0;
posCorrect = sum(handles.scores(pos)>=handles.T1);
posIncorrect = sum(handles.scores(pos)<=handles.T2);
posAbstain = sum(handles.scores(pos)<handles.T1 & handles.scores(pos)>handles.T2);
negCorrect = sum(handles.scores(~pos)<=handles.T2);
negIncorrect = sum(handles.scores(~pos)>=handles.T1);
negAbstain = sum(handles.scores(~pos)<handles.T1 & handles.scores(~pos)>handles.T2);
textStr = '';
textStr = sprintf('                      Predicted Positive      Abstained        Predicted Negative\n');
textStr = sprintf('%s Labeled Pos      %5d                 %5d                %5d      \n',...
  textStr,posCorrect, posAbstain,posIncorrect);
textStr = sprintf('%s Labeled Neg      %5d                 %5d                %5d      \n',...
  textStr,negIncorrect, negAbstain,negCorrect);
textStr = sprintf('%sAUC:%.4f',textStr,handles.AUC);
set(handles.text1,'String',textStr);

% --- Outputs from this function are returned to the command line.
function varargout = ShowROCCurve_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
