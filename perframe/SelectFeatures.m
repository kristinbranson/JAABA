function varargout = SelectFeatures(varargin)
% SELECTFEATURES MATLAB code for SelectFeatures.fig
%      SELECTFEATURES, by itself, creates a new SELECTFEATURES or raises the existing
%      singleton*.
%
%      H = SELECTFEATURES returns the handle to a new SELECTFEATURES or the handle to
%      the existing singleton*.
%
%      SELECTFEATURES('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECTFEATURES.M with the given input arguments.
%
%      SELECTFEATURES('Property','Value',...) creates a new SELECTFEATURES or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SelectFeatures_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SelectFeatures_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SelectFeatures

% Last Modified by GUIDE v2.5 29-Nov-2011 12:54:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SelectFeatures_OpeningFcn, ...
                   'gui_OutputFcn',  @SelectFeatures_OutputFcn, ...
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


% --- Executes just before SelectFeatures is made visible.
function SelectFeatures_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectFeatures (see VARARGIN)

% Choose default command line output for SelectFeatures
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SelectFeatures wait for user response (see UIRESUME)


% --- Outputs from this function are returned to the command line.
function varargout = SelectFeatures_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function createWindowTable(hObject)
% Sets values for the window table.

handles = guidata(hObject);

% Deal with windowTable
handles.windowComp = {'Default','Mean','Min','Max','Histogram','Prctile',...
   'Change','Std-dev','Harmonic','Diff Neighbor Mean',...
   'Diff Neighbor Min','Diff Neighbor Max','Zscore Neighbors'};

 set(handles.windowTable,'RowName', handles.windowComp,...
    'ColumnName',{'Computation Type','Select'});
 selVals = [true true true true false false false false false false false false false];
 tableData = {};
 for ndx = 1:numel(handles.windowComp)
   tableData{ndx,1} = handles.windowComp{ndx};
   tableData{ndx,2} = selVals(ndx);
 end
 set(handles.windowTable,'Data',...
 tableData);
set(handles.windowTable','ColumnWidth',{135 'auto'});
set(handles.pfTable,'ColumnEditable',[false,true]);

% Use java objects to tweak the table. Found this online at 
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/298335

jscroll=findjobj(handles.windowTable);
rowHeaderViewport=jscroll.getComponent(4);
rowHeader=rowHeaderViewport.getComponent(0);
rowHeader.setSize(80,360);

%resize the row header
newWidth=0; %100 pixels.
rowHeaderViewport.setPreferredSize(java.awt.Dimension(newWidth,0));
height=rowHeader.getHeight;
rowHeader.setPreferredSize(java.awt.Dimension(newWidth,height));
rowHeader.setSize(newWidth,height); 

guidata(hObject,handles);

function createPfTable(hObject)
% Sets the values for feature table.

handles = guidata(hObject);
pfList = handles.JLDobj.GetAllPerframeFeatures();
[params cellParams] = handles.JLDobj.GetPerframeParams();
handles.params = params;
handles.cellParams = cellParams;

curParams = fieldnames(params);

tableData = cell(length(pfList),2);
for ndx = 1:numel(pfList)
  tableData{ndx,1} = pfList{ndx};
  paramNdx = find(strcmp(pfList{ndx},curParams));
  if paramNdx
    tableData{ndx,2} = true;
  else
    tableData{ndx,2} = false;
  end
end

handles.windowData = windowData;

set(handles.pfTable,'Data',tableData);
set(handles.pfTable,'ColumnName',{'Features','Select'});
set(handles.pfTable,'ColumnEditable',[false,true]);

% Make the table sortable. Use underlying java objects to do that. Found
% this at http://undocumentedmatlab.com/blog/uitable-sorting/
jscrollpane = findjobj(handles.pfTable);
jtable = jscrollpane.getViewport.getView;
jtable.setSortable(true);	
jtable.setAutoResort(false);
jtable.setMultiColumnSortable(true);

% Set the size for the row headers.
rowHeaderViewport=jscrollpane.getComponent(4);
rowHeader=rowHeaderViewport.getComponent(0);
newWidth=0; 
rowHeaderViewport.setPreferredSize(java.awt.Dimension(newWidth,0));
height=rowHeader.getHeight;
rowHeader.setPreferredSize(java.awt.Dimension(newWidth,height));
rowHeader.setSize(newWidth,height); 
set(handles.pfTable,'ColumnWidth',{200,85});

set(handles.pfTable,'CellSelectionCallback',@pfSelect);

guidata(hObject,handles);


function pfSelect(hObject,eventData)
% Called when user selects cells in pfTable.

handles = guidata(hObject);
jscrollpane = findjobj(handles.pfTable);
jtable = jscrollpane.getViewport.getView;
pfData = get(handles.pfTable,'Data');

if(size(eventData.Indices,1)>1)
  handles.multiSelect = true;
  handles.disableValues = true;
else
  handles.multiSelect = false;
end

if ~handles.multiSelect
  ndx = jtable.getActualRowAt(eventData.Indices(1,1));
  handles.pfNdx = ndx;
  if(pfData{ndx,2})
    handles.disableValues = false;
    setWindowTable(handles);
  end
  
end
guidata(hObject,handles);
  


function createWindowData(handles)
handles = guidata(hObject);
pfData = get(handles.pfTable,'Data');
winFunc = handles.windowComp;
origSelPf = fieldnames(handles.params);

winData = cell(size(pfData,1),numel(winFunc));
for ndx = 1:numel(pfData)
  pfFunc = pfData{ndx,1};
  paramNdx = find(strcmp(pfList{ndx},origSelPf));
  if ~paramNdx, continue; end;
  curCellParams = handles.cellParams.(origSelPf(paramNdx));
  names = curCellParams{1:2:end};
  
  for winNdx = 1:numel(winFunc)
    curFunc = winFunc(winNdx);
    curStr = [curFunc '_params'] % The names are different!!
    
  end
  
end
  

function setWindowTable(handles)
% Sets the data in window table based on selections made in pfTable.
curData = handles.windowData(handles.pfNdx);

curFeatureTypes = curData{1}{end};
windowData = {};
for ndx = 1:length(handles.windowComp)
  curFunc = handles.windowComp{ndx};
  windowData{ndx,1} = curFunc;
  
  if any(strcmp(curFunc,curFeatureTypes))
    
  end
  
end

function setJLDobj(hObject,JLDobj)
% Set the JLabelDataObject.
handles = guidata(hObject);
handles.JLDobj = JLDobj;
guidata(hObject,handles);


function MinWindow_Callback(hObject, eventdata, handles)
% hObject    handle to MinWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinWindow as text
%        str2double(get(hObject,'String')) returns contents of MinWindow as a double


% --- Executes during object creation, after setting all properties.
function MinWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MinWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MaxWindow_Callback(hObject, eventdata, handles)
% hObject    handle to MaxWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MaxWindow as text
%        str2double(get(hObject,'String')) returns contents of MaxWindow as a double


% --- Executes during object creation, after setting all properties.
function MaxWindow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MaxWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WindowStep_Callback(hObject, eventdata, handles)
% hObject    handle to WindowStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowStep as text
%        str2double(get(hObject,'String')) returns contents of WindowStep as a double


% --- Executes during object creation, after setting all properties.
function WindowStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function WindowOffsets_Callback(hObject, eventdata, handles)
% hObject    handle to WindowOffsets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowOffsets as text
%        str2double(get(hObject,'String')) returns contents of WindowOffsets as a double


% --- Executes during object creation, after setting all properties.
function WindowOffsets_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WindowOffsets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in TransFlip.
function TransFlip_Callback(hObject, eventdata, handles)
% hObject    handle to TransFlip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransFlip


% --- Executes on button press in TransAbs.
function TransAbs_Callback(hObject, eventdata, handles)
% hObject    handle to TransAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransAbs


% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to push_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
delete(handles.figure1);


% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);
delete(handles.figure1);


% --- Executes on button press in pushbutton_applydefault.
function pushbutton_applydefault_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_applydefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
