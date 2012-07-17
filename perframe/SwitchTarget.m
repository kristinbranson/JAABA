function varargout = SwitchTarget(varargin)
% SWITCHTARGET MATLAB code for SwitchTarget.fig
%      SWITCHTARGET, by itself, creates a new SWITCHTARGET or raises the existing
%      singleton*.
%
%      H = SWITCHTARGET returns the handle to a new SWITCHTARGET or the handle to
%      the existing singleton*.
%
%      SWITCHTARGET('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SWITCHTARGET.M with the given input arguments.
%
%      SWITCHTARGET('Property','Value',...) creates a new SWITCHTARGET or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SwitchTarget_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SwitchTarget_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SwitchTarget

% Last Modified by GUIDE v2.5 23-May-2012 10:45:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SwitchTarget_OpeningFcn, ...
                   'gui_OutputFcn',  @SwitchTarget_OutputFcn, ...
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


% --- Executes just before SwitchTarget is made visible.
function SwitchTarget_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SwitchTarget (see VARARGIN)

% Choose default command line output for SwitchTarget
handles.output = hObject;

handles.JLabelhObject = varargin{1};
[handles.ntargets_per_page] = myparse(varargin(2:end),'ntargets_per_page',100);

JLabelHandles = guidata(handles.JLabelhObject);
handles.JLDObj = JLabelHandles.guidata.data;
% set(handles.axes1,'axes',off);
handles.tablePos = get(handles.table,'Position');
handles.switchFlyPos = get(handles.pushSwitchfly,'Position');
handles.updatePos = get(handles.pushbutton_update,'Position');
handles.closePos = get(handles.pushClose,'Position');
handles.figurePos = get(handles.figure1,'Position');
handles.pagePos = get(handles.popupmenu_Page,'Position');

% initialize page number
handles.page_number = [];

% added an extra guidata here, as sometimes I was getting an error that
% figurePos wasn't getting set
guidata(hObject,handles);
handles = updateTable(handles);
guidata(hObject,handles);

function [exp,target] = GlobalTarget2LocalTarget(handles,globaltarget)

exp = find(handles.cs_ntargets >= globaltarget,1,'first');
if exp == 1,
  target = globaltarget;
else
  target = globaltarget - handles.cs_ntargets(exp-1);
end

function [globaltarget] = LocalTarget2GlobalTarget(handles,exp,target)

if exp == 1,
  globaltarget = target;
else
  globaltarget = handles.cs_ntargets(exp-1)+target;
end

function [page,row] = GlobalTarget2PageRow(handles,globaltarget)

page = ceil(globaltarget/handles.ntargets_per_page);
row = globaltarget - (page-1)*handles.ntargets_per_page;

function [globaltarget] = PageRow2GlobalTarget(handles,page,row)

globaltarget = handles.ntargets_per_page*(page-1)+row;

function [exp,target] = PageRow2LocalTarget(handles,page,row)

[exp,target] = GlobalTarget2LocalTarget(handles,PageRow2GlobalTarget(handles,page,row));

function handles = updateTable(handles)
% Initialize the table

handles.curExp = handles.JLDObj.GetExp();
handles.curFly = handles.JLDObj.GetFlies();

% which fly out of all flies
handles.cs_ntargets = cumsum(handles.JLDObj.nflies_per_exp);
handles.curFlyGlobal = LocalTarget2GlobalTarget(handles,handles.curExp,handles.curFly);

% which page, row
[handles.current_target_page_number,handles.current_target_row] = GlobalTarget2PageRow(handles,handles.curFlyGlobal);

% set page number if not already set
if isempty(handles.page_number),
  handles.page_number = handles.current_target_page_number;
end

% set page number options
handles.npages = ceil(handles.cs_ntargets(end)/handles.ntargets_per_page);
handles.pagestrs = cell(1,handles.npages);
for i = 1:handles.npages,
  exp1 = PageRow2LocalTarget(handles,i,1);
  tmp = min(PageRow2GlobalTarget(handles,i,handles.ntargets_per_page),handles.cs_ntargets(end));
  exp2 = GlobalTarget2LocalTarget(handles,tmp);
  if exp1 >= exp2,
    handles.pagestrs{i} = sprintf('%d (Exp %d)',i,exp1);
  else
    handles.pagestrs{i} = sprintf('%d (Exps %d-%d)',i,exp1,exp2);
  end
  if i == handles.current_target_page_number,
    handles.pagestrs{i} = ['* ',handles.pagestrs{i},' *'];
  end
  
end

if handles.page_number > handles.npages,
  handles.page_number = handles.npages;
end

set(handles.popupmenu_Page,'String',handles.pagestrs,'Value',handles.page_number);

fieldList = {};
if ~handles.JLDObj.IsGTMode(),
  fieldList(end+1,:) = {  'trajLength','Trajectory| Length'};
  fieldList(end+1,:) = {  'firstframe','Start frame'};
  fieldList(end+1,:) = {  'nbouts','Bouts|labeled'};
  fieldList(end+1,:) = {  'totalframes','Frames|labeled'};
  fieldList(end+1,:) = {  'posframes',sprintf('Frames|labeled|%s',handles.JLDObj.labelnames{1})};
  fieldList(end+1,:) = {  'negframes',sprintf('Frames|labeled|%s',handles.JLDObj.labelnames{2})};
  fieldList(end+1,:) = {  'sexfrac','Sex|(% male)'};
  fieldList(end+1,:) = {  'nscorepos',  sprintf('%s| predicted | (current)',handles.JLDObj.labelnames{1})};
  fieldList(end+1,:) = {  'nscoreneg',  sprintf('%s| predicted | (current)',handles.JLDObj.labelnames{2})};
  fieldList(end+1,:) = {  'nscoreframes',  'total|predicted | frames | (current)'};
  fieldList(end+1,:) = {  'nscorepos_loaded',  sprintf('%s| predicted | (loaded)',handles.JLDObj.labelnames{1})};
  fieldList(end+1,:) = {  'nscoreneg_loaded',  sprintf('%s| predicted | (loaded)',handles.JLDObj.labelnames{2})};
  fieldList(end+1,:) = {  'nscoreframes_loaded',  'total|predicted | frames | (loaded)'};
  fieldList(end+1,:) = {  'errorsPos',  sprintf('%s|errors',handles.JLDObj.labelnames{1})};
  fieldList(end+1,:) = {  'errorsNeg',  sprintf('%s|errors',handles.JLDObj.labelnames{2})};
  fieldList(end+1,:) = {  'validatedErrorsPos',  sprintf('%s|validation|errors',handles.JLDObj.labelnames{1})};
  fieldList(end+1,:) = {  'validatedErrorsNeg',  sprintf('%s|validation|errors',handles.JLDObj.labelnames{2})};
  fieldList(end+1,:) = {  'one2two',  sprintf('%s|switched to|%s',handles.JLDObj.labelnames{1},handles.JLDObj.labelnames{2})};
  fieldList(end+1,:) = {  'two2one',  sprintf('%s|switched to|%s',handles.JLDObj.labelnames{2},handles.JLDObj.labelnames{1})};
else
  fieldList(end+1,:) = {  'trajLength','Trajectory| Length'};
  fieldList(end+1,:) = {  'firstframe','Start frame'};
  fieldList(end+1,:) = {  'sexfrac','Sex|(% male)'};
  fieldList(end+1,:) = {  'gt_nbouts','Ground Truthing|Bouts|labeled'};
  fieldList(end+1,:) = {  'gt_totalframes','Ground Truthing|Frames|labeled'};
  fieldList(end+1,:) = {  'gt_posframes',sprintf('Ground Truthing|Frames|labeled|%s',handles.JLDObj.labelnames{1})};
  fieldList(end+1,:) = {  'gt_negframes',sprintf('Ground Truthing|Frames|labeled|%s',handles.JLDObj.labelnames{2})};
  fieldList(end+1,:) = {  'gt_suggestion_frames','Ground Truthing|Frames|Suggested'};
end  

colFormat{1} = 'char';
colFormat{2} = 'numeric';
for ndx = 1:size(fieldList,1)
  colFormat{2+ndx} = 'numeric';
end

% reorganize cached data
newCachedTableData = cell(handles.cs_ntargets(end),size(fieldList,1)+2);
newCachedDataExpi = zeros(1,handles.cs_ntargets(end));
newCachedDataTarget = zeros(1,handles.cs_ntargets(end));
newCachedExpDirs = handles.JLDObj.expdirs;

if isfield(handles,'cachedTableData'),

  % reindexing for experiments
  [ism,old2newexpi] = ismember(handles.cachedExpDirs,newCachedExpDirs);

  % move data from old to new
  for oldexpi = find(ism(:)'),
    newexpi = old2newexpi(oldexpi);
    oldidx = find(handles.cachedDataExpi == oldexpi);
    if isempty(oldidx),
      continue;
    end
    [targets_curr,order] = sort(handles.cachedDataTarget(oldidx));
    oldidx = oldidx(order);
    newidx = LocalTarget2GlobalTarget(handles,newexpi,targets_curr);
    newCachedTableData(newidx,:) = handles.cachedTableData(oldidx,:);
    newCachedDataExpi(newidx) = newexpi;
    newCachedDataTarget(newidx) = handles.cachedDataTarget(oldidx);
  end
  
end
  
% overwrite
handles.cachedTableData = newCachedTableData;
handles.cachedDataExpi = newCachedDataExpi;
handles.cachedDataTarget = newCachedDataTarget;
handles.cachedExpDirs = newCachedExpDirs;

tableData = {};
tableExpi = [];
tableTarget = [];
count = 1;

[exp_start,target_start] = PageRow2LocalTarget(handles,handles.page_number,1);
tmp = min(PageRow2GlobalTarget(handles,handles.page_number,handles.ntargets_per_page),handles.cs_ntargets(end));
[exp_end,target_end] = GlobalTarget2LocalTarget(handles,tmp);

for selExp = exp_start:exp_end,
  if selExp == exp_start,
    startcurr = target_start;
  else
    startcurr = 1;
  end
  endcurr = handles.JLDObj.nflies_per_exp(selExp);
  if selExp == exp_end,
    endcurr = min(endcurr,target_end);
  end
  for flyNum = startcurr:endcurr,
    flyStats = handles.JLDObj.GetFlyStats(selExp,flyNum);
    tableData{count,1} = handles.JLDObj.expnames{selExp};
    tableData{count,2} = flyNum;
    for ndx = 1:size(fieldList,1)
      tableData{count,2+ndx} = flyStats.(fieldList{ndx,1});
    end
    tableExpi(count) = selExp;
    tableTarget(count) = flyNum;
    count = count+1;
  end
end
handles.tableData = tableData;

% store in global cache
start_target = PageRow2GlobalTarget(handles,handles.page_number,1);
end_target = start_target + numel(tableExpi) - 1;
handles.cachedTableData(start_target:end_target,:) = tableData;
handles.cachedDataExpi(start_target:end_target) = tableExpi;
handles.cachedDataTarget(start_target:end_target) = tableTarget;

set(handles.table,'Data',...
  tableData);
set(handles.table,'ColumnName',...
  {'Experiment Name','Target|Number',fieldList{:,2}});
set(handles.table,'ColumnEditable',false);
set(handles.table,'ColumnFormat',colFormat);
set(handles.table,'CellSelectionCallback',@tableSelect);


% Update handles structure

% UIWAIT makes SwitchTarget wait for user response (see UIRESUME)
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
function varargout = SwitchTarget_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes when figure1 is resized.
function figure1_ResizeFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'figurePos'),
  return;
end

oldClosePosRight = handles.figurePos(3)-handles.closePos(1);
guiPos = get(handles.figure1,'Position');
newClosePos = [guiPos(3)-oldClosePosRight handles.closePos(2:4)];
set(handles.pushClose,'Position',newClosePos);

oldSwitchPosRight = handles.figurePos(3)-handles.switchFlyPos(1);
newSwitchFlyPos = [guiPos(3)-oldSwitchPosRight handles.switchFlyPos(2:4)];
set(handles.pushSwitchfly,'Position',newSwitchFlyPos);

oldUpdatePosRight = handles.figurePos(3)-handles.updatePos(1);
newUpdatePos = [guiPos(3)-oldUpdatePosRight handles.updatePos(2:4)];
set(handles.pushbutton_update,'Position',newUpdatePos);

newTableWidth = guiPos(3)- handles.figurePos(3)+handles.tablePos(3);
newTableHeight = guiPos(4)- handles.figurePos(4)+handles.tablePos(4);
set(handles.table,'Position',[handles.tablePos(1:2) newTableWidth newTableHeight]);


% --- Executes on button press in pushbutton_update.
function pushbutton_update_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_update (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = updateTable(handles);
guidata(hObject,handles);


% --- Executes on selection change in popupmenu_Page.
function popupmenu_Page_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_Page (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_Page contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_Page

handles.page_number = get(hObject,'Value');
start_target = PageRow2GlobalTarget(handles,handles.page_number,1);
end_target = min(handles.cs_ntargets(end),start_target+handles.ntargets_per_page-1);
if all(handles.cachedDataExpi(start_target:end_target) ~= 0),
  handles.tableData = handles.cachedTableData(start_target:end_target,:);  
  set(handles.table,'Data',handles.tableData);
else
  handles = updateTable(handles);
end
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function popupmenu_Page_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_Page (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
