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

% Last Modified by GUIDE v2.5 03-Jan-2012 14:11:19

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
handles.windowComp = {'default','mean','min','max','histogram','prctile',...
   'change','std','harmonic','diff_neighbor_mean',...
   'diff_neighbor_min','diff_neighbor_max','zscore_neighbors'};

handles.winParams = {'max_window_radius','min_window_radius','nwindow_radii',...
  'trans_types','window_offsets'};
handles.defaultWinParams = {0,0,0,'',0};
 
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
set(handles.windowTable,'ColumnEditable',[false,true]);
set(handles.windowTable,'CellSelectionCallback',@windowSelect);
set(handles.windowTable,'CellEditCallback',@windowEdit);

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

initData(hObject);
handles = guidata(hObject);

pfList = handles.pfList;
tableData = cell(length(pfList),2);
for ndx = 1:numel(pfList)
  tableData{ndx,1} = pfList{ndx};
  tableData{ndx,2} = handles.data{ndx}.valid;
end

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
set(handles.pfTable,'CellEditCallback',@pfEdit);
disableWindowTable(handles);
guidata(hObject,handles);


% Initialize the data structure.
function initData(hObject)
handles = guidata(hObject);
pfList =  handles.JLDobj.GetAllPerframeFeatures();
handles.pfList = pfList;
data = {};
[params,~] = handles.JLDobj.GetPerframeParams();
validPfs = fieldnames(params);
winParams = handles.winParams;
for ndx = 1:numel(pfList)
  curPfName = pfList{ndx};
  data{ndx}.name = curPfName;
  pNdx = find(strcmp(curPfName,validPfs));
  if pNdx
    data{ndx}.valid = true;
    data{ndx}.sanitycheck = params.(curPfName).sanitycheck;
    
    % Fill the default values.
    for winParamsNdx = 1:numel(winParams)
      curType = winParams{winParamsNdx};
      if isfield(params.(curPfName),curType)
        data{ndx}.default.values.(curType) = params.(curPfName).(curType);
      else % Fill in the default value
        % TODO.
        msgbox(sprintf('%s default window parameters value were not defined',curType));
      end
      data{ndx}.default.valid = true;
      data{ndx}.default.values.extra ='';
    end
    
    % Fill for different window function type.
    
    % Find which ones are valid.
    curParam = params.(curPfName);
    validWinfn = fieldnames(curParam);
    for winfnNdx = 2:numel(handles.windowComp)
      curFn = handles.windowComp{winfnNdx};
      wNdx = find(strcmp(curFn,validWinfn));
      
      if wNdx,
        
        data{ndx}.(curFn).valid = true;
        curWinFnParams = params.(curPfName).(curFn);
        for winParamsNdx = 1:numel(winParams)
          curType = winParams{winParamsNdx};
          if isfield(curWinFnParams,curType)
            data{ndx}.(curFn).values.(curType) = curWinFnParams.(curType);
          else % fill in the default values
            data{ndx}.(curFn).values.(curType) = data{ndx}.default.values.(curType);
          end
        end
        
        data{ndx}.(curFn).values.extra = findExtraParams(curWinFnParams,winParams);
        
      else % Values for window comp hasn't been defined.
        
        data{ndx}.(curFn).valid = false;
        for winParamsNdx = 1:numel(winParams)
          curType = winParams{winParamsNdx};
          data{ndx}.(curFn).values.(curType) = data{ndx}.default.values.(curType);
        end
        
      end
      
    end
  
  else % Default values for invalid pf's.
    data{ndx}.valid = false;
    data{ndx}.sanitycheck = false;
  end
end

handles.data = data;
guidata(hObject,handles);



function str = findExtraParams(curParam,winComp)
% Find feature specific params.
allFields = fieldnames(curParam);
str = '';
for ndx = 1:numel(allFields)
  curF = allFields{ndx};
  if any(strcmp(curF,winComp)); continue;end
  vals = mat2str(curParam.(curF));
  str = sprintf('%s%s:%s\n',str,curF,vals);
end

function [names vals] = convertStrToExtraParams(str)
names = {}; vals = {};
str = strtrim(str);
if isempty(str); return; end
s = regexp(str,'\n','split');
for ndx = 1:numel(s)
  ss = regexp(strtrim(s{ndx}),':','split');
  if (numel(ss)~=2)
    msgbox(sprintf('Invalid syntax in feature specific params:%s',s{ndx}));
  end
  names{ndx} = strtrim(ss{1});
  vals{ndx} = str2num(ss{2});
end

function pfSelect(hObject,eventData)
% Called when user selects cells in pfTable.

% When the table is sorted without any cell selected
if isempty(eventData.Indices)
  return;
end

handles = guidata(hObject);
jscrollpane = findjobj(handles.pfTable);
jtable = jscrollpane.getViewport.getView;
pfData = get(handles.pfTable,'Data');

if(size(eventData.Indices,1)>1)
  disableWindowTable(handles);
else
  ndx = jtable.getActualRowAt(eventData.Indices(1,1)-1)+1;
  handles.pfNdx = ndx;
  if(pfData{ndx,2})
    setWindowTable(handles,handles.pfNdx);
  else
    disableWindowTable(handles);
  end
end

guidata(hObject,handles);


function pfEdit(hObject,eventData)
% When a perframe feature is added or removed.

handles = guidata(hObject);
jscrollpane = findjobj(handles.pfTable);
jtable = jscrollpane.getViewport.getView;
pfNdx = jtable.getActualRowAt(eventData.Indices(1,1)-1)+1;
handles.pfNdx = pfNdx;

handles.data{pfNdx}.valid = eventData.NewData;

% When a feature is unchecked.
if ~eventData.NewData, 
    disableWindowTable(handles);  
    return;
end

% When it already has values
if isfield(handles.data{pfNdx},'default')
  enableWindowTable(handles);
  return;
end

% Else fill in the default values.
handles.data{pfNdx}.sanitycheck = false;
handles.data{pfNdx}.default.valid = true;
for winParamsNdx = 1:numel(handles.winParams)
  curType = handles.winParams{winParamsNdx};
  handles.data{pfNdx}.default.values.(curType) = ...
    handles.defaultWinParams{winParamsNdx};
end
handles.data{pfNdx}.default.values.extra = '';

% Copy the default values into the other window params.
for winfnNdx = 2:numel(handles.windowComp)
  curFn = handles.windowComp{winfnNdx};
  handles.data{pfNdx}.(curFn).valid = false;
  for winParamsNdx = 1:numel(handles.winParams)
    curType = handles.winParams{winParamsNdx};
    handles.data{pfNdx}.(curFn).values.(curType) = ...
      handles.data{pfNdx}.default.values.(curType);
  end
  handles.data{pfNdx}.(curFn).values.extra = '';
end
guidata(hObject,handles);
setWindowTable(handles,handles.pfNdx);


function setWindowTable(handles,pfNdx)
% Sets the data in window table based on selections made in pfTable.
enableWindowTable(handles);

curData = handles.data{pfNdx};

windowData = {};
for ndx = 1:length(handles.windowComp)
  curFunc = handles.windowComp{ndx};
  windowData{ndx,1} = curFunc;
  windowData{ndx,2} = curData.(curFunc).valid;
end
set(handles.windowTable,'Data',windowData);

jscrollpane = findjobj(handles.windowTable);
jtable = jscrollpane.getViewport.getView;
jtable.changeSelection(0,0, false, false);


function windowSelect(hObject,eventData)
% Called when user selects cells in windowTable.

% When something changes in the pfTable window.
if isempty(eventData.Indices)
  return;
end

handles = guidata(hObject);
winData = get(handles.windowTable,'Data');

if(size(eventData.Indices,1)>1)
  disableWinParams(handles);
else
  ndx = eventData.Indices(1,1);
  handles.winNdx = ndx;
  if winData{ndx,2},
    setWinParams(handles,handles.winNdx);
  else
    disableWinParams(handles);
  end
end

guidata(hObject,handles);


function windowEdit(hObject,eventData)
% Called when a window function is edited.

handles = guidata(hObject);
winNdx = eventData.Indices(1);
handles.winNdx = winNdx;
curFn = handles.windowComp{winNdx};
handles.data{handles.pfNdx}.(curFn).valid = eventData.NewData;
guidata(hObject,handles);
if eventData.NewData
  setWinParams(handles,winNdx);
else
  disableWinParams(handles);
end


function setWinParams(handles,winNdx)
enableWinParams(handles);
curFn = handles.windowComp{winNdx};
curParams = handles.data{handles.pfNdx}.(curFn);
set(handles.MinWindow,'String',num2str(curParams.values.min_window_radius));
set(handles.MaxWindow,'String',num2str(curParams.values.max_window_radius));
set(handles.WindowStep,'String',num2str(curParams.values.nwindow_radii));
set(handles.WindowOffsets,'String',num2str(curParams.values.window_offsets));
set(handles.TransFlip,'Value',...
  any(strcmp('flip',curParams.values.trans_types)));
set(handles.TransAbs,'Value',...
  any(strcmp('abs',curParams.values.trans_types)));
% TODO: For relative trans.
set(handles.TransRel,'Value',...
  any(strcmp('relative',curParams.values.trans_types)));
set(handles.ExtraParams,'String',curParams.values.extra);

function setJLDobj(hObject,JLDobj)
% Set the JLabelDataObject.
handles = guidata(hObject);
handles.JLDobj = JLDobj;
guidata(hObject,handles);


function [params, cellparams] = convertData(handles)
% Converts the data into format used by JLabelData.

params = struct;

for ndx = 1:numel(handles.pfList)
  curPf = handles.pfList{ndx};
  if ~handles.data{ndx}.valid; continue;end
  curD = handles.data{ndx};
  params.(curPf).sanitycheck = curD.sanitycheck;
  
  for winParamsNdx = 1:numel(handles.winParams)
    curType = handles.winParams{winParamsNdx};
    params.(curPf).(curType) = curD.default.values.(curType);
  end
  
  for winfnNdx = 2:numel(handles.windowComp)
    curFn = handles.windowComp{winfnNdx};
    
    if curD.(curFn).valid,
      for winParamsNdx = 1:numel(handles.winParams)
        curType = handles.winParams{winParamsNdx};
        params.(curPf).(curFn).(curType) =...
            curD.(curFn).values.(curType);
      end
      [names vals] = convertStrToExtraParams(curD.(curFn).values.extra);
      for extraNdx = 1:numel(names)
        params.(curPf).(curFn).(names{extraNdx}) = vals{extraNdx};
      end
    end
    
  end

end

% Code from ReadPerFrameParams.
cellparams = struct;
fns1 = fieldnames(params);
for i1 = 1:numel(fns1),
  fn1 = fns1{i1};
  fns2 = fieldnames(params.(fn1));
  cellparams.(fn1) = {};
  feature_types = {};
  for i2 = 1:numel(fns2),
    fn2 = fns2{i2};
    if ~isstruct(params.(fn1).(fn2)),
      cellparams.(fn1)(end+1:end+2) = {fn2,params.(fn1).(fn2)};
    else
      cellparams.(fn1)(end+1:end+2) = {[fn2,'_params'],struct2paramscell(params.(fn1).(fn2))};
      feature_types{end+1} = fn2; %#ok<AGROW>
    end
  end
  cellparams.(fn1)(end+1:end+2) = {'feature_types',feature_types};
end

function createParamsXML(params)
docNode = com.mathworks.xml.XMLUtils.createDocument('params');
toc = docNode.getDocumentElement;
toc.setAttribute('version','1.0');



function disableWindowTable(handles)
set(handles.windowTable,'enable','off');
disableWinParams(handles);

function enableWindowTable(handles)
set(handles.windowTable,'enable','on');

function disableWinParams(handles)
set(handles.MinWindow,'enable','off');
set(handles.MaxWindow,'enable','off');
set(handles.WindowStep,'enable','off');
set(handles.WindowOffsets,'enable','off');
set(handles.TransFlip,'enable','off');
set(handles.TransAbs,'enable','off');
set(handles.TransRel,'enable','off');
set(handles.ExtraParams,'enable','off');

function enableWinParams(handles)
defBack = get(handles.text5,'BackgroundColor');
defFore = [0 0 0];
set(handles.MinWindow,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.MaxWindow,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.WindowStep,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.WindowOffsets,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.ExtraParams,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.TransFlip,'enable','on');
set(handles.TransAbs,'enable','on');
set(handles.TransRel,'enable','on');

function MinWindow_Callback(hObject, eventdata, handles)
% hObject    handle to MinWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinWindow as text
%        str2double(get(hObject,'String')) returns contents of MinWindow as a double
curVal = str2double(get(hObject,'String'));
if isnan(curVal)
  msgbox('Enter numerical values');
end
handles = guidata(hObject);
curFn = handles.windowComp{handles.winNdx};
handles.data{handles.pfNdx}.(curFn).values.min_window_radius = curVal;
guidata(hObject,handles);


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
curVal = str2double(get(hObject,'String'));
if isnan(curVal)
  msgbox('Enter numerical values');
end
handles = guidata(hObject);
curFn = handles.windowComp{handles.winNdx};
handles.data{handles.pfNdx}.(curFn).values.max_window_radius = curVal;
guidata(hObject,handles);


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
curVal = str2double(get(hObject,'String'));
if isnan(curVal)
  msgbox('Enter numerical values');
end
handles = guidata(hObject);
curFn = handles.windowComp{handles.winNdx};
handles.data{handles.pfNdx}.(curFn).values.nwindow_radii = curVal;
guidata(hObject,handles);


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
curVal = str2num(get(hObject,'String'));
if isempty(curVal)
  msgbox('Enter numerical values. eg: "-1 0 1" (without with quotes)');
end
handles = guidata(hObject);
curFn = handles.windowComp{handles.winNdx};
handles.data{handles.pfNdx}.(curFn).values.window_offsets = curVal;
guidata(hObject,handles);


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
curFn = handles.windowComp{handles.winNdx};
handles = guidata(hObject);
curT = handles.data{handles.pfNdx}.(curFn).values.trans_types;
if get(hObject,'Value')
  if ~any(strcmp('flip',curT))
    handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'flip';
  end
else
  allNdx = find(strcmp('flip',curT));
    handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
end
guidata(hObject,handles);


% --- Executes on button press in TransAbs.
function TransAbs_Callback(hObject, eventdata, handles)
% hObject    handle to TransAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransAbs

curFn = handles.windowComp{handles.winNdx};
handles = guidata(hObject);
curT = handles.data{handles.pfNdx}.(curFn).values.trans_types;
if get(hObject,'Value')
  if ~any(strcmp('abs',curT))
    handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'abs';
  end
else
  allNdx = find(strcmp('abs',curT));
    handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
end
guidata(hObject,handles);

% --- Executes on button press in TransRel.
function TransRel_Callback(hObject, eventdata, handles)
% hObject    handle to TransRel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransRel

curFn = handles.windowComp{handles.winNdx};
handles = guidata(hObject);
curT = handles.data{handles.pfNdx}.(curFn).values.trans_types;
if get(hObject,'Value')
  if ~any(strcmp('relative',curT))
    handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'relative';
  end
else
  allNdx = find(strcmp('relative',curT));
    handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
end
guidata(hObject,handles);



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
handles = guidata(hObject);
[params,cellParams] = convertData(handles);
handles.JLDobj.UpdatePerframeParams(params,cellParams);
uiresume(handles.figure1);
delete(handles.figure1);


% --- Executes on button press in pushbutton_applydefault.
function pushbutton_applydefault_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_applydefault (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

function ExtraParams_Callback(hObject, eventdata, handles)
% hObject    handle to ExtraParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ExtraParams as text
%        str2double(get(hObject,'String')) returns contents of ExtraParams as a double


% --- Executes during object creation, after setting all properties.
function ExtraParams_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ExtraParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
