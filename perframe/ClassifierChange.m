function varargout = ClassifierChange(varargin)
% CLASSIFIERCHANGE MATLAB code for ClassifierChange.fig
%      CLASSIFIERCHANGE, by itself, creates a new CLASSIFIERCHANGE or raises the existing
%      singleton*.
%
%      H = CLASSIFIERCHANGE returns the handle to a new CLASSIFIERCHANGE or the handle to
%      the existing singleton*.
%
%      CLASSIFIERCHANGE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CLASSIFIERCHANGE.M with the given input arguments.
%
%      CLASSIFIERCHANGE('Property','Value',...) creates a new CLASSIFIERCHANGE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ClassifierChange_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ClassifierChange_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ClassifierChange

% Last Modified by GUIDE v2.5 24-Jan-2012 12:28:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ClassifierChange_OpeningFcn, ...
                   'gui_OutputFcn',  @ClassifierChange_OutputFcn, ...
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


% --- Executes just before ClassifierChange is made visible.
function ClassifierChange_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ClassifierChange (see VARARGIN)

% Choose default command line output for ClassifierChange
handles.output = hObject;

handles.JLabelhObject = varargin{1};
JLabelHandles = guidata(handles.JLabelhObject);
handles.JLDObj = JLabelHandles.data;
% set(handles.axes1,'axes',off);
handles.curExp = JLabelHandles.expi;
handles.curFly = JLabelHandles.flies;


% Initialize the table
set(handles.table,'ColumnName',...
  {'Experiment Name','Fly Number','Trajectory Length',...
  'Start frame','Bouts labeled','Sex (% male)',...
  sprintf('%s predicted',handles.JLDObj.labelnames{1}),...
  sprintf('%s predicted',handles.JLDObj.labelnames{2}),...
  'total predicted frames',...
  sprintf('%s errors',handles.JLDObj.labelnames{1}),...
  sprintf('%s errors',handles.JLDObj.labelnames{2}),...
  sprintf('%s switched to %s',handles.JLDObj.labelnames{1},handles.JLDObj.labelnames{2}),...
  sprintf('%s switched to %s',handles.JLDObj.labelnames{2},handles.JLDObj.labelnames{1}),...
  });
set(handles.table,'ColumnFormat',...
  {'char','numeric','numeric',...
  'numeric','numeric','numeric',...
  'numeric','numeric','numeric','numeric',...
  'numeric','numeric'...
  });

tableData = {};
count = 1;
for selExp = 1:handles.JLDObj.nexps
  for flyNum = 1:handles.JLDObj.nflies_per_exp(selExp)
    flyStats = handles.JLDObj.GetFlyStats(selExp,flyNum);
    tableData{count,1} = handles.JLDObj.expnames{selExp};
    tableData{count,2} = flyNum;
    tableData{count,3} = flyStats.trajLength;
    tableData{count,4} = flyStats.firstframe;
    tableData{count,5} = flyStats.nbouts;
    if flyStats.hassex
      if ~isempty(flyStats.sexfrac),
        tableData{count,6} = flyStats.sexfrac.M*100;
      else
        if strcmpi(flyStats.sex,'M')
          tableData{count,6} = 100;
        else
          tableData{count,6} = 0;
        end
      end
    else
      tableData{count,6} = 0;
    end
    if ~isempty(flyStats.nscoreframes),
      tableData{count,7} = flyStats.nscorepos;
      tableData{count,8} = flyStats.nscoreframes-flyStats.nscorepos;
      tableData{count,9} = flyStats.nscoreframes;
    else
      tableData{count,7} = 0;
      tableData{count,8} = 0;
      tableData{count,9} = 0;
    end
    
    if ~isempty(flyStats.errorsPos),
      tableData{count,10} = flyStats.errorsPos;
      tableData{count,11} = flyStats.errorsNeg;
    else
      tableData{count,10} = 0;
      tableData{count,11} = 0;
    end
    
    if ~isempty(flyStats.one2two),
      tableData{count,12} = flyStats.one2two;
      tableData{count,13} = flyStats.two2one;
    else
      tableData{count,12} = 0;
      tableData{count,13} = 0;
    end
    count = count+1;
  end
end
handles.tableData = tableData;
set(handles.table,'Data',...
  tableData);
set(handles.table,'ColumnEditable',false);
set(handles.table,'CellSelectionCallback',@tableSelect);




guidata(hObject, handles);
% Update handles structure

% UIWAIT makes ClassifierChange wait for user response (see UIRESUME)
% uiwait(handles.figure1);


function initTable(hObject)
% Use java objects to tweak the table. Found this online at 
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/298335

handles = guidata(hObject);
jscrollpane = findjobj(handles.table);
jtable = jscrollpane.getViewport.getView;
jtable.setSortable(false);	

colWidth = repmat({70},1,13);
colWidth{1} = 300;
set(handles.table,'ColumnWidth',colWidth);


% rowHeaderViewport=jscrollpane.getComponent(4);
% rowHeader=rowHeaderViewport.getComponent(0);
% rowHeader.setSize(80,360);
% 
% %resize the row header
% newWidth=0; %100 pixels.
% rowHeaderViewport.setPreferredSize(java.awt.Dimension(newWidth,0));
% height=rowHeader.getHeight;
% rowHeader.setPreferredSize(java.awt.Dimension(newWidth,height));
% rowHeader.setSize(newWidth,height); 


% --- Outputs from this function are returned to the command line.
function varargout = ClassifierChange_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function tableSelect(hObject,eventData)

handles = guidata(hObject);
jscrollpane = findjobj(handles.table);
jtable = jscrollpane.getViewport.getView;

if(size(eventData.Indices,1)==1)
  ndx = jtable.getActualRowAt(eventData.Indices(1,1)-1)+1;
  handles.curExp = find(strcmp(handles.JLDObj.expnames,handles.tableData{ndx,1}));
  handles.curFly = handles.tableData{ndx,2};
end

guidata(hObject,handles);


% --- Executes on button press in pushSwitchfly.
function pushSwitchfly_Callback(hObject, eventdata, handles)
% hObject    handle to pushSwitchfly (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
JLabelHandles = guidata(handles.JLabelhObject);
[JLabelHandles,~] = JLabel('SetCurrentMovie',JLabelHandles,handles.curExp);
guidata(handles.JLabelhObject,JLabelHandles);
JLabelHandles = JLabel('SetCurrentFlies',JLabelHandles,handles.curFly);
guidata(handles.JLabelhObject,JLabelHandles);


% --- Executes on button press in pushClose.
function pushClose_Callback(hObject, eventdata, handles)
% hObject    handle to pushClose (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);
