function varargout = PostProcess(varargin)
% POSTPROCESS MATLAB code for PostProcess.fig
%      POSTPROCESS, by itself, creates a new POSTPROCESS or raises the existing
%      singleton*.
%
%      H = POSTPROCESS returns the handle to a new POSTPROCESS or the handle to
%      the existing singleton*.
%
%      POSTPROCESS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in POSTPROCESS.M with the given input arguments.
%
%      POSTPROCESS('Property','Value',...) creates a new POSTPROCESS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PostProcess_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PostProcess_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PostProcess

% Last Modified by GUIDE v2.5 24-Jul-2012 16:04:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PostProcess_OpeningFcn, ...
                   'gui_OutputFcn',  @PostProcess_OutputFcn, ...
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


% --- Executes just before PostProcess is made visible.
function PostProcess_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PostProcess (see VARARGIN)

% Choose default command line output for PostProcess
handles.output = hObject;

set(handles.uipanel_method,'SelectionChangeFcn',@methodchange)

handles.data = varargin{1};
handles.JLabelHandle = varargin{2};
params = handles.data.GetPostprocessingParams();
% params = [];

if isempty(params)
  params.method = 'Hysteresis';
  params.hystopts(1) = struct('name','High Threshold','tag','hthres','value',0);
  params.hystopts(2) = struct('name','Low Threshold','tag','lthres','value',0);
  params.filtopts(1) = struct('name','Size','tag','size','value',1);
%   params.filtopts(2) = struct('name','High Threshold','tag','filt_hthres','value',0);
%   params.filtopts(3) = struct('name','Low Threshold','tag','filt_lthres','value',0);
  params.blen = 1;
end

if strcmp(params.method,'Hysteresis')
  set(handles.uipanel_method,'SelectedObject',handles.radiobutton_hyst);
else
  set(handles.uipanel_method,'SelectedObject',handles.radiobutton_filt);
end  

handles.params = params;
handles.origparams = params;

handles.paramsui = [];
handles = UpdateParams(handles);

set(handles.edit_blen,'String',sprintf('%d',params.blen));
pushbutton_update_Callback(hObject,eventdata,handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PostProcess wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PostProcess_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton_update.
function pushbutton_update_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.SetPostprocessingParams(handles.params);
blen = handles.data.GetPostprocessedBoutLengths();

JLabel('UpdateTimelineIms',handles.JLabelHandle);
JLabel('UpdatePlots',handles.JLabelHandle,'refreshim',false,'refreshflies',false,...
  'refreshtrx',false,'refreshlabels',true,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_auto',true,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);


if isempty(blen), 
  axis(handles.axes1,'off');
else
  
  % wsize = handles.data.GetFeatureWindowSize();
  %wsize = handles.data.featureWindowSize();
  wsize=10;
  edges = unique([1 2 round(wsize/2) wsize 2*wsize 4*wsize 8*wsize inf]);
  
  vals = histc(blen,edges);
  vals(end) = [];
  bar(handles.axes1,1:numel(edges)-1,vals);
  xlocs = 1:numel(edges)-1;
  xlabels = {}; xlabels{1} = '1';
  for ndx = 2:numel(edges)-1
    xlabels{ndx} = sprintf('%d-%d',edges(ndx),edges(ndx+1)-1);
  end
  set(handles.axes1,'XTick',xlocs,'XTickLabel',xlabels);
  xlabel(handles.axes1,'Bout Length');
  ylabel(handles.axes1,'Bouts');
  
end

[labels,lscores,allscores,scoreNorm] = handles.data.GetAllLabelsAndScores();
if isempty(allscores), 
  axis(handles.axes2,'off'); 
  return; 
end
numBins = 21;
if isnan(scoreNorm); scoreNorm = 1; end;
pos = labels==1;
neg = ~pos;
bins = linspace(-scoreNorm,scoreNorm,numBins);
bins = [-inf bins(2:end-1) inf];
histPos = histc(lscores(pos),bins);
histNeg = histc(lscores(neg),bins);
if ~isempty(histPos),histPos(end) = []; 
else histPos = zeros(numBins-1,1);end
if ~isempty(histNeg),histNeg(end) = []; 
else histNeg = zeros(numBins-1,1);end

handles.posColor = handles.JLabelHandle.guidata.labelcolors(1,:);
handles.negColor = handles.JLabelHandle.guidata.labelcolors(2,:);

% For axis 1
xLocs = linspace(-1+1/(numBins-1),1-1/(numBins-1),numBins-1);
hold(handles.axes2,'off');
hBar = bar(handles.axes2, xLocs,[histPos histNeg],'BarWidth',1.5);
set(hBar(1),'FaceColor',handles.posColor);
set(hBar(2),'FaceColor',handles.negColor);
ylim = get(handles.axes2,'ylim');
if(ylim(1)<0); 
  ylim(1) = 0; 
  set(handles.axes2,'ylim',ylim);
end
set(handles.axes2,'xlim',[-1 1]);
xlabel(handles.axes2,'Scores');
ylabel(handles.axes2,'Frames');

hold(handles.axes2,'on');
histscores = histc(allscores,bins);
if ~isempty(histscores),histscores(end) = [];
histscores = histscores./max(histscores)*0.9*ylim(2);
plot(handles.axes2,xLocs,histscores,'Color',max(handles.posColor,handles.negColor));
end;
handles.JLabelHandle = JLabel('SetNeedSave',handles.JLabelHandle);
guidata(handles.JLabelHandle.figure_JLabel,handles.JLabelHandle);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);

% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.data.SetPostprocessingParams(handles.params);
%handles.data.ApplyPostprocessing();
  % now done inside SetPostprocessingParams()
JLabel('UpdateTimelineIms',handles.JLabelHandle);
JLabel('UpdatePlots',handles.JLabelHandle,'refreshim',false,'refreshflies',false,...
  'refreshtrx',false,'refreshlabels',true,...
  'refresh_timeline_manual',false,...
  'refresh_timeline_auto',true,...
  'refresh_timeline_xlim',false,...
  'refresh_timeline_hcurr',false,...
  'refresh_timeline_selection',false,...
  'refresh_curr_prop',false);
handles.JLabelHandle = JLabel('SetNeedSave',handles.JLabelHandle);
guidata(handles.JLabelHandle.figure_JLabel,handles.JLabelHandle);
delete(handles.figure1);


function methodchange(hObject,eventdata)

method = get(eventdata.NewValue,'String');
handles = guidata(hObject);

handles.params.method = method;

handles = UpdateParams(handles);
guidata(hObject,handles);

function handles = UpdateParams(handles)

for ndx = 1:numel(handles.paramsui)
  delete(handles.paramsui(ndx));
end

if strcmp(handles.params.method,'Hysteresis');
  opts = handles.params.hystopts;
else
  opts = handles.params.filtopts;
end

psz = get(handles.uipanel_params,'Position');
ui = [];
for ndx = 1:numel(opts)
  ypos = psz(4) - 20 - ndx*50;
  ui(end+1) = uicontrol(handles.uipanel_params,'Style','text','String',opts(ndx).name);
  set(ui(end),'Units','pixels');
  set(ui(end),'Position',[5 ypos psz(3)/2-15 40]);

  ui(end+1) = uicontrol(handles.uipanel_params,'Style','edit','String',opts(ndx).value);
  set(ui(end),'Units','pixels');
  set(ui(end),'Position',[psz(3)/2+5 ypos psz(3)/2-15 40]);
  set(ui(end),'Value',ndx,'Callback',@paramsEdit);
end
handles.paramsui = ui;

function paramsEdit(hObject,eventdata)

val = str2num(get(hObject,'String'));
if isempty(val) || isnan(val) 
  uiwait(warndlg('Enter non-negative numerical values'));
end

ndx = get(hObject,'Value');
handles = guidata(hObject);
if strcmp(handles.params.method,'Hysteresis')
  mtype = 'hystopts';
else
  mtype = 'filtopts';
end

handles.params.(mtype)(ndx).value = val;
guidata(hObject,handles);



function edit_blen_Callback(hObject, eventdata, handles)
% hObject    handle to edit_blen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_blen as text
%        str2double(get(hObject,'String')) returns contents of edit_blen as a double
val = str2double(get(hObject,'String'));
if isempty(val) || isnan(val) || (round(val)-val)>0 || val<0
  uiwait(warndlg('Minimum bout length should be a positive integer value'))
end
handles.params.blen = val;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function edit_blen_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_blen (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
