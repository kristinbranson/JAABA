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

% Last Modified by GUIDE v2.5 12-Jul-2012 14:02:14

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
params = handles.data.GetPostprocessingParams();
% params = [];

if isempty(params)
  params.method = 'Hysteresis';
  params.hystopts(1) = struct('name','High Threshold','tag','hthres','value',0);
  params.hystopts(2) = struct('name','Low Threshold','tag','lthres','value',0);
  params.filtopts(1) = struct('name','Size','tag','size','value',1);
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

if isempty(blen), axis(handles.axes1,'off'); return; end

% wsize = handles.data.GetFeatureWindowSize();
wsize = handles.data.featureWindowSize();
edges = [1 2 round(wsize/2) wsize 2*wsize 4*wsize 8*wsize inf];


vals = histc(blen,edges);
vals(end) = [];
bar(handles.axes1,1:numel(edges)-1,vals);
xlocs = 1:numel(edges)-1;
xlabels = {}; xlabels{1} = '1';
for ndx = 2:numel(edges)-1
  xlabels{ndx} = sprintf('%d-%d',edges(ndx),edges(ndx+1)-1);
end
set(handles.axes1,'XTick',xlocs,'XTickLabel',xlabels);

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
