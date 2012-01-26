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
[xdata,ydata,T,AUC] = perfcurve(handles.labels,handles.scores,1);
handles.xdata = xdata;
handles.ydata = ydata;
handles.T = T;
handles.AUC = AUC

hLine = plot(handles.axes1,xdata, ydata);
% First get the figure's data-cursor mode, activate it, and set some of its properties
cursorMode = datacursormode(hObject);
set(cursorMode, 'enable','on','NewDataCursorOnClick',false);
% Note: the optional @setDataTipTxt is used to customize the data-tip's appearance
% Note: the following code was adapted from %matlabroot%\toolbox\matlab\graphics\datacursormode.m
% Create a new data tip
hTarget = handle(hLine);
hDatatip = cursorMode.createDatatip(hTarget);
% Create a copy of the context menu for the datatip:
% set(hDatatip,'UIContextMenu',get(cursorMode,'UIContextMenu'));
set(hDatatip,'UIContextMenu',[]);
set(hDatatip,'HandleVisibility','off');
set(hDatatip,'Host',hTarget);
set(hDatatip,'ViewStyle','datatip');
% Set the data-tip orientation to top-right rather than auto
% set(hDatatip,'OrientationMode','manual');
set(hDatatip,'Orientation','bottom-right');
set(hDatatip,'UpdateFcn',@setDataTipTxtFN);
% Update the datatip marker appearance
set(hDatatip, 'MarkerSize',15, 'MarkerFaceColor','none', ...
'MarkerEdgeColor','k', 'Marker','.', 'HitTest','off');
% Move the datatip to the right-most data vertex point
position = [xdata(end-10),ydata(end-10),1; xdata(end-10),ydata(end-10),-1];
handles.T1 = T(end-10);
update(hDatatip, position);

hDatatip2 = cursorMode.createDatatip(hTarget);
% Create a copy of the context menu for the datatip:
set(hDatatip2,'UIContextMenu',[]);
set(hDatatip2,'HandleVisibility','off');
set(hDatatip2,'Host',hTarget);
set(hDatatip2,'ViewStyle','datatip');
% Set the data-tip orientation to top-right rather than auto
set(hDatatip2,'OrientationMode','manual');
set(hDatatip2,'Orientation','bottom-right');
set(hDatatip2,'UpdateFcn',@setDataTipTxtFP);
% Update the datatip marker appearance
set(hDatatip2, 'MarkerSize',5, 'MarkerFaceColor','none', ...
'MarkerEdgeColor','k', 'Marker','o', 'HitTest','off');
% Move the datatip to the right-most data vertex point
position = [xdata(10),ydata(10),1; xdata(10),ydata(10),-1];
handles.T2 = T(10);
update(hDatatip2, position);

UpdateText(handles);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ShowROCCurve wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function output_txt = setDataTipTxtFP(~,eventData)
% handles = guidata(eventData.Target);
handles = guidata(get(eventData,'Target'));
position = get(eventData,'Position');

curNdx = min(sum(handles.ydata<=position(2)),sum(handles.xdata<=position(1)));
handles.T2 = handles.T(curNdx);
output_txt = sprintf('TP:%.2f',position(2));
UpdateText(handles);
guidata(handles.output,handles);


function output_txt = setDataTipTxtFN(~,eventData)
% handles = guidata(eventData.Target);
handles = guidata(get(eventData,'Target'));
position = get(eventData,'Position');

curNdx = min(sum(handles.ydata<=position(2)),sum(handles.xdata<=position(1)));
handles.T1 = handles.T(curNdx);
output_txt = sprintf('FP:%.2f',position(1));
UpdateText(handles);
guidata(handles.output,handles);

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
textStr = sprintf('%s,T1:%.2f  T2:%.2f',textStr,handles.T1,handles.T2);
set(handles.text1,'String',textStr);

% --- Outputs from this function are returned to the command line.
function varargout = ShowROCCurve_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
