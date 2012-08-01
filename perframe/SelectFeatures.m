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

% Last Modified by GUIDE v2.5 09-Mar-2012 10:29:28

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
set(hObject,'Visible','off');
JLDobj = varargin{1};
handles.output = hObject;

handles.mode = 'basic';
curPos = get(handles.figure1,'Position');
tablePos = get(handles.basicTable,'Position');
reducedWidth = tablePos(1)+tablePos(3) + 15;
reducedHeight = curPos(4);
handles.advancedSize = curPos(3:4);
handles.basicSize = [reducedWidth reducedHeight];
set(handles.figure1,'Position',[curPos(1:2) handles.basicSize]);
% Update handles structure
guidata(hObject, handles);

set(handles.pfTable,'UserData',0);
setJLDobj(hObject,JLDobj);
set(hObject,'Visible','on');
removeRowHeaders(hObject); pause(0.5);
% uiwait(hObject);
% UIWAIT makes SelectFeatures wait for user response (see UIRESUME)


% --- Outputs from this function are returned to the command line.
function varargout = SelectFeatures_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function setJLDobj(hObject,JLDobj)
% Set the JLabelDataObject.
handles = guidata(hObject);
handles.JLDobj = JLDobj;

handles.windowComp = {'default','mean','min','max','hist','prctile',...
   'change','std','harmonic','diff_neighbor_mean',...
   'diff_neighbor_min','diff_neighbor_max','zscore_neighbors'};

handles.winextraParams = {'','','','','hist_edges','prctiles','change_window_radii',...
  '','num_harmonic','','','',''};

handles.winParams = {'max_window_radius','min_window_radius','nwindow_radii',...
  'trans_types','window_offsets'};

handles.winextraDefaultParams = {[],[],[],[],[-400000 0 40000],[5 10 30 50 70 90 95],[1 3],...
  [],[2],[],[],[],[]};
handles.defaultWinParams = {10,1,3,{'none'},0};

guidata(hObject,handles);

[params,~] = handles.JLDobj.GetPerframeParams();
initData(hObject,params);


function createWindowTable(hObject)
% Sets values for the window table.

handles = guidata(hObject);

% Deal with windowTable
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

guidata(hObject,handles);


function createPfTable(hObject)
% Sets the values for feature table.

handles = guidata(hObject);

pfList = handles.pfList;
tableData = cell(length(pfList),2);
for ndx = 1:numel(pfList)
  tableData{ndx,1} = pfList{ndx};
  tableData{ndx,2} = handles.data{ndx}.valid;
  str = sprintf('%s',handles.pftype.(pfList{ndx}){1});
  for sndx = 2:numel(handles.pftype.(pfList{ndx}))
    str = sprintf('%s,%s',str,handles.pftype.(pfList{ndx}){sndx});
  end
  tableData{ndx,3} = str;
end

set(handles.pfTable,'Data',tableData);
set(handles.pfTable,'ColumnName',{'Features','Select','Category'});
set(handles.pfTable,'ColumnEditable',[false,true]);

set(handles.pfTable,'ColumnWidth',{190,50,95});
set(handles.pfTable,'CellSelectionCallback',@pfSelect);
set(handles.pfTable,'CellEditCallback',@pfEdit);
disableWindowTable(handles);
guidata(hObject,handles);

function createFeatureTable(hObject)
handles = guidata(hObject);

% Adding drop down menu kind of thing.

if ~isempty(handles.JLDobj.basicFeatureTable)
  tableData = handles.JLDobj.basicFeatureTable;
else
  tableData = {};
  for ndx = 1:numel(handles.pftypeList)
    tableData{ndx,1} = handles.pftypeList{ndx};
    tableData{ndx,2} = 'Custom';
    tableData{ndx,3} = 'normal';
  end
end
set(handles.basicTable,'Data',tableData);
set(handles.basicTable,'ColumnName',{'Categories','Select','Amount'});
set(handles.basicTable,'ColumnFormat',{'char',...
  {'All' 'None' 'Custom'}, fieldnames(handles.categ)'});
set(handles.basicTable,'ColumnEditable',[false,true,true]);

set(handles.basicTable,'ColumnWidth',{85,65,75});
set(handles.basicTable,'CellSelectionCallback',@basicSelect);
set(handles.basicTable,'CellEditCallback',@basicEdit);

function removeRowHeaders(hObject)
% Tweaking the table. Use underlying java objects to do that. Found
% this at http://undocumentedmatlab.com/blog/uitable-sorting/

handles = guidata(hObject);

% Basic Table
jscrollpane = findjobj(handles.basicTable);
jtable = jscrollpane.getViewport.getView;
jtable.setSortable(false);	
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

% Pf Table.
jscrollpane = findjobj(handles.pfTable);
jtable = jscrollpane.getViewport.getView;
jtable.setSortable(false);	
jtable.setAutoResort(false);
jtable.setMultiColumnSortable(false);

% Set the size for the row headers.
rowHeaderViewport=jscrollpane.getComponent(4);
rowHeader=rowHeaderViewport.getComponent(0);
newWidth=0; 
rowHeaderViewport.setPreferredSize(java.awt.Dimension(newWidth,0));
height=rowHeader.getHeight;
rowHeader.setPreferredSize(java.awt.Dimension(newWidth,height));
rowHeader.setSize(newWidth,height); 

% Window Table.
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


function createCopyFromMenus(hObject)

handles = guidata(hObject);
% can copy from any window feature type
set(handles.popupmenu_copy_windowparams,'String',...
  [handles.windowComp,{'---'}],...
  'Value',numel(handles.windowComp)+1);
% can copy from any per-frame feature
set(handles.popupmenu_copy_windowtypes,'String',...
  [handles.pfList;{'---'}],...
  'Value',numel(handles.pfList)+1);

function createDescriptionPanels(hObject)

handles = guidata(hObject);

% which tab to show
handles.currentTab = 'description';
set(handles.togglebutton_tabdescription,'Value',1);

% set visibility
uipanel_tabs_SelectionChangeFcn(handles.uipanel_tabs, struct('NewValue',handles.togglebutton_tabdescription), handles);

guidata(hObject,handles);

% Initialize the data structure.
function initData(hObject,params)
readFeatureConfiguration(hObject);
createFeatureTable(hObject);

handles = guidata(hObject);
if ~isempty(handles.JLDobj.featureWindowSize)
  set(handles.editSize,'String',num2str(handles.JLDobj.featureWindowSize));
end

pfList = fieldnames(handles.pftype);
handles.pfList = pfList;
data = {};
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
        data{pfNdx}.default.values.(curType) = ...
            handles.defaultWinParams{winParamsNdx};
      end
      data{ndx}.default.valid = true;
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
        
        if ~isempty(handles.winextraParams{winfnNdx})
          extraParam = handles.winextraParams{winfnNdx};
          data{ndx}.(curFn).values.(extraParam) = curWinFnParams.(extraParam);
        end
        
      else % Values for window comp haven't been defined.
        
        data{ndx}.(curFn).valid = false;
        for winParamsNdx = 1:numel(winParams)
          curType = winParams{winParamsNdx};
          data{ndx}.(curFn).values.(curType) = data{ndx}.default.values.(curType);
        end
        if ~isempty(handles.winextraParams{winfnNdx})
          extraParam = handles.winextraParams{winfnNdx};
          data{ndx}.(curFn).values.(extraParam) = handles.winextraDefaultParams{winfnNdx};
        end
        
      end
      
    end
  
  else % Default values for invalid pf's.
    
    data{ndx}.valid = false;
    data{ndx}.sanitycheck = false;

    data{ndx}.default.valid = true;
    for winParamsNdx = 1:numel(handles.winParams)
      curType = handles.winParams{winParamsNdx};
      data{ndx}.default.values.(curType) = ...
        handles.defaultWinParams{winParamsNdx};
    end
    
    % Copy the default values into the other window params.
    for winfnNdx = 2:numel(handles.windowComp)
      curFn = handles.windowComp{winfnNdx};
      data{ndx}.(curFn).valid = false;
      for winParamsNdx = 1:numel(handles.winParams)
        curType = handles.winParams{winParamsNdx};
        data{ndx}.(curFn).values.(curType) = ...
          data{ndx}.default.values.(curType);
      end
        if ~isempty(handles.winextraParams{winfnNdx})
          extraParam = handles.winextraParams{winfnNdx};
          data{ndx}.(curFn).values.(extraParam) = handles.winextraDefaultParams{winfnNdx};
        end
    end

    
  end
end

% initialize histogramData
handles.histogramData = struct;
handles.histogramData.lastPfNdx = nan;
handles.histogramData.lastType = '';
handles.histogramData.perframe_idx = [];
handles.histogramData.hhist = [];
handles.histogramData.frac = {};
handles.histogramData.frac_outside = {};
handles.histogramData.edges = {};
handles.histogramData.centers_plot = {};

handles.data = data;
handles.pfNdx = [];
handles.winNdx = [];
guidata(hObject,handles);
createPfTable(hObject);
createWindowTable(hObject);
createCopyFromMenus(hObject);
createDescriptionPanels(hObject);
compatibleBasicAdvanced(handles);

function readFeatureConfiguration(hObject)
% Reads the configuration file that sets the default parameters.

handles = guidata(hObject);

configfile = handles.JLDobj.featureConfigFile;
settings = ReadXMLParams(configfile);

% Read the default parameters for different categories.
categories = fieldnames(settings.defaults);

if isempty(categories) 
  % If no default have been specified in the config file.
  
  categories{1} = 'default';
  curParams = struct;
  
  % Default values.
  for ndx = 1:numel(handles.winParams)
    curwinpname = handles.winParams{ndx};
    curParams.default.values.(curwinpname) = handles.defaultWinParams{ndx};
  end
  curParams.default.valid = true;
  curParams.default.values.sanitycheck = false;
  
  % For each defined window computation.
  for ndx = 2:numel(handles.windowComp)
    % Copy the default values.
    curwname = handles.windowComp{ndx};
    for pndx = 1:numel(handles.winParams)
      curwinpname = handles.winParams{pndx};
      curParams.(curwname).values.(curwinpname) = handles.defaultWinParams{pndx};
    end
    
    if ~isempty(handles.winextraParams{ndx})
      extraParam = handles.winextraParams{ndx};
      curParams.(curwname).values.(extraParam) = handles.winextraParams{ndx};
    end
    curParams.(curwname).valid = false;
  end
  handles.categ.(categories{1}) = curParams;
  
  
else
  % fill the window params for different feature categories.
  
  for j = 1:numel(categories)
    curParams = struct;
    cur = settings.defaults.(categories{j});
    
    % Default values.
    for ndx = 1:numel(handles.winParams)
      curwinpname = handles.winParams{ndx};
      curParams.default.values.(curwinpname) = cur.(curwinpname);
    end
    curParams.default.valid = true;
    curParams.default.values.sanitycheck = false;
    
    % For each defined window computation.
    for ndx = 2:numel(handles.windowComp)
      % Copy the default values.
      curwname = handles.windowComp{ndx};
      for pndx = 1:numel(handles.winParams)
        curwinpname = handles.winParams{pndx};
        curParams.(curwname).values.(curwinpname) = curParams.default.values.(curwinpname);
      end
      
      % Override the default window params for the current window computation type
      if isfield(cur,curwname)
        curParams.(curwname).valid = true;
        diffFields = fieldnames(cur.(curwname));
        for dndx = 1:numel(diffFields)
          if any(strcmp(handles.winParams,diffFields{dndx}))
            curwinpname = diffFields{dndx};
            curParams.(curwname).values.(curwinpname) = cur.(curwname).(curwinpname);
          end
        end
        
        if ~isempty(handles.winextraParams{ndx})
          extraParam = handles.winextraParams{ndx};
          if isfield(cur.(curwname),extraParam)
            curParams.(curwname).values.(extraParam) = cur.(curwname).(extraParam);
          else
            curParams.(curwname).values.(extraParam) = '';
          end
        end
      else
        curParams.(curwname).valid = false;
      end
    end
    handles.categ.(categories{j}) = curParams;
  end
end

% Now for each different perframe feature read the type default trans_types.
perframeL = handles.JLDobj.allperframefns;

transType = struct;
pftype = struct;

% perframeL might not contain all the perframe features.
for pfndx = 1:numel(perframeL)
  curpf = perframeL{pfndx};
  transType.(curpf) = settings.perframe.(curpf).trans_types;
  curtypes  = settings.perframe.(curpf).type; 
  if ischar(curtypes)
    pftype.(curpf)  = {curtypes}; 
  else    
    pftype.(curpf)  = curtypes; 
  end
end

fallpf = fieldnames(settings.perframe);
pftypeList = {};
for pfndx = 1:numel(fallpf)
  curpf = fallpf{pfndx};
  curtypes  = settings.perframe.(curpf).type; 
  if ischar(curtypes)
    curT = curtypes;
    if ~any(strcmp(pftypeList,curT))
      pftypeList{end+1} = curT;
    end
  else    
    for tndx = 1:numel(curtypes)
      curT = curtypes{tndx};
      if ~any(strcmp(pftypeList,curT))
        pftypeList{end+1} = curT;
      end
    end
  end
end


handles.transType = transType;
handles.pftype = pftype;
handles.pftypeList = pftypeList;
guidata(hObject,handles);
set(handles.editSize,'String',...
  num2str(handles.categ.(categories{1}).(handles.windowComp{1}).values.max_window_radius));

function basicSelect(hObject,eventData)

function basicEdit(hObject,eventData)

% the user selects the category
if isempty(eventData.Indices); return; end

handles = guidata(hObject);
if eventData.Indices(2) ==2
% Select-deselect the whole category.  
  switch eventData.NewData

    case 'None'
      h = waitbar(0,'Unselecting the perframe features');
      basicTable = get(handles.basicTable,'Data');
      for ndx = 1:numel(handles.pfList)
        disable = true;
        for tndx = 1:size(basicTable,1)
          if ~any(strcmp(handles.pftype.(handles.pfList{ndx}),basicTable{tndx,1}));
            continue;
          end
          if ~strcmp(basicTable{tndx,2},'None'); 
            disable = false; break; 
          end
        end
        if disable
          handles.data{ndx}.valid = false;
        end
      end
      guidata(hObject,handles);
      createPfTable(handles.pfTable);
      close(h);
    case 'All'
      h = waitbar(0,'Selecting the perframe features');
      handles = applyCategoryType(handles,eventData.Indices(1));
      guidata(hObject,handles);
      createPfTable(handles.pfTable);
      close(h);
  end
elseif eventData.Indices(2) == 3
% Choose the category.
  basicTable = get(handles.basicTable,'Data');
  if strcmp(basicTable{eventData.Indices(1),2},'All') 
    handles=applyCategoryType(handles,eventData.Indices(1));
    guidata(hObject,handles);
    createPfTable(handles.pfTable);
  end
end

function handles = applyCategoryType(handles,basicNdx)
basicData = get(handles.basicTable,'Data');
curType = handles.pftypeList{basicNdx};
% Copy the parameters.
for ndx = 1:numel(handles.pfList)
  if any(strcmp(handles.pftype.(handles.pfList{ndx}),curType))
    handles.data{ndx}.valid = true;
    for winfnNdx = 1:numel(handles.windowComp)
      category = basicData{basicNdx,3};
      curFn = handles.windowComp{winfnNdx};
      if ~handles.categ.(category).(curFn).valid
        handles.data{ndx}.(curFn).valid = false;
        continue;
      end
      handles = CopyDefaultWindowParams(handles,...
        category, ndx,winfnNdx);
      curFn = handles.windowComp{winfnNdx};
      handles.data{ndx}.(curFn).values.trans_types = handles.transType.(handles.pfList{ndx});
    end
  end
end

function compatibleBasicAdvanced(handles)
basicTable = get(handles.basicTable,'Data');
incompatible = '';
for ndx = 1:numel(handles.pfList)
  selected = false;
  for bndx = 1:size(basicTable,1)
    if ~any(strcmp(handles.pftype.(handles.pfList{ndx}),basicTable{bndx,1})),
      continue;
    end
    if strcmpi(basicTable{bndx,2},'all')
       selected = true;
    end
    if selected && ~handles.data{ndx}.valid
      incompatible = sprintf('%s %s', incompatible,handles.pfList{ndx});
    end
    
  end
end
if numel(incompatible)>0
  uiwait(warndlg(sprintf('Perframe feature(s) %s should have been selected but are not',incompatible),...
    'Mismatch in categories and perframe features'));
end

function pfSelect(hObject,eventData)
% Called when user selects cells in pfTable.

% When the table is sorted without any cell selected

% while( get(hObject,'UserData')~=0)
%   pause(0.2);
% end
% set(hObject,'UserData',1);
if isempty(eventData.Indices)
  return;
end

handles = guidata(hObject);
% jscrollpane = findjobj(handles.pfTable);
% jtable = jscrollpane.getViewport.getView;
pfData = get(handles.pfTable,'Data');

if(size(eventData.Indices,1)>1)
  disableWindowTable(handles);
else
%   ndx = jtable.getActualRowAt(eventData.Indices(1,1)-1)+1;
  ndx = eventData.Indices(1,1);
  handles.pfNdx = ndx;
  if(pfData{ndx,2})
    setWindowTable(handles,handles.pfNdx);
    handles.winNdx = 1;
    winData = get(handles.windowTable,'Data');

    if winData{handles.winNdx,2},
      setWinParams(handles,handles.winNdx);
      enableWinParams(handles);
    else
      disableWinParams(handles);
    end

    setWinParams(handles,handles.winNdx);
  else
    disableWindowTable(handles);
  end
end
handles = UpdateDescriptionPanels(handles);

guidata(hObject,handles);
% set(hObject,'UserData',0);


function pfEdit(hObject,eventData)
% When a perframe feature is added or removed.

% while( get(hObject,'UserData')~=0)
%   pause(0.2);
% end
% set(hObject,'UserData',2);

handles = guidata(hObject);
% jscrollpane = findjobj(handles.pfTable);
% jtable = jscrollpane.getViewport.getView;
% pfNdx = jtable.getActualRowAt(eventData.Indices(1,1)-1)+1;
pfNdx = eventData.Indices(1,1);
handles.pfNdx = pfNdx;

handles.data{pfNdx}.valid = eventData.NewData;
setCategoryToCustom(handles);

% When a feature is unchecked, disable and return
if ~eventData.NewData, 
  handles.data{pfNdx}.valid = false;
  disableWindowTable(handles);  
  guidata(hObject,handles);
  return;
end

handles.data{pfNdx}.valid = true;
curType = handles.pftype.(handles.pfList{pfNdx}){1};
basicNdx = find(strcmp(handles.pftypeList,curType));
basicData = get(handles.basicTable,'Data');
handles.data{pfNdx}.valid = true;
for winfnNdx = 1:numel(handles.windowComp)
  category = basicData{basicNdx,3};
  curFn = handles.windowComp{winfnNdx};
  if ~handles.categ.(category).(curFn).valid
    handles.data{pfNdx}.(curFn).valid = false;
    continue;
  end
  handles = CopyDefaultWindowParams(handles,...
    category, pfNdx,winfnNdx);
  curFn = handles.windowComp{winfnNdx};
  handles.data{pfNdx}.(curFn).values.trans_types = handles.transType.(handles.pfList{pfNdx});
end


% 
% % When it already has values
% if isfield(handles.data{pfNdx},'default')
%   guidata(hObject,handles);
%   setWindowTable(handles,handles.pfNdx);
%   return;
% end
% 
% % Else fill in the default values.
% handles.data{pfNdx}.sanitycheck = false;
% handles.data{pfNdx}.default.valid = true;
% for winParamsNdx = 1:numel(handles.winParams)
%   curType = handles.winParams{winParamsNdx};
%   handles.data{pfNdx}.default.values.(curType) = ...
%     handles.defaultWinParams{winParamsNdx};
% end
% 
% % Copy the default values into the other window params.
% for winfnNdx = 2:numel(handles.windowComp)
%   curFn = handles.windowComp{winfnNdx};
%   handles.data{pfNdx}.(curFn).valid = false;
%   for winParamsNdx = 1:numel(handles.winParams)
%     curType = handles.winParams{winParamsNdx};
%     handles.data{pfNdx}.(curFn).values.(curType) = ...
%       handles.data{pfNdx}.default.values.(curType);
%   end
%   if ~isempty(handles.winextraParams{winfnNdx})
%     extraParam = handles.winextraParams{winfnNdx};
%     handles.data{pfNdx}.(curFn).values.(extraParam) = '';
%   end
%   
% end


guidata(hObject,handles);
setWindowTable(handles,handles.pfNdx);
% set(hObject,'UserData',0);


function setWindowTable(handles,pfNdx)
% Sets the data in window table based on selections made in pfTable.
enableWindowTable(handles);

curData = handles.data{pfNdx};

windowData = {};
for ndx = 1:length(handles.windowComp)
  curFunc = handles.windowComp{ndx};
  if ~isfield(curData,curFunc),
    warning('This error is occurring!');
    windowData{ndx,1} = curFunc;
    windowData{ndx,2} = false;
  else
    windowData{ndx,1} = curFunc;
    windowData{ndx,2} = curData.(curFunc).valid;
  end
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
    enableWinParams(handles);
  else
    disableWinParams(handles);
  end
end

guidata(hObject,handles);


function windowEdit(hObject,eventData)
% Called when a window function is edited.

handles = guidata(hObject);
setCategoryToCustom(handles);
winNdx = eventData.Indices(1);
handles.winNdx = winNdx;
curFn = handles.windowComp{winNdx};
handles.data{handles.pfNdx}.(curFn).valid = eventData.NewData;
guidata(hObject,handles);
if eventData.NewData
  setWinParams(handles,winNdx);
  enableWinParams(handles);
else
  disableWinParams(handles);
end

function setCategoryToCustom(handles)
curTypes = find(ismember(handles.pftypeList,handles.pftype.(handles.pfList{handles.pfNdx})));
basicData = get(handles.basicTable,'Data');
for ndx = 1:numel(curTypes)
  basicData{curTypes(ndx),2} = 'Custom';
end
set(handles.basicTable,'Data',basicData);


function setWinParams(handles,winNdx)

curFn = handles.windowComp{winNdx};
curParams = handles.data{handles.pfNdx}.(curFn);
set(handles.MinWindow,'String',num2str(curParams.values.min_window_radius));
set(handles.MaxWindow,'String',num2str(curParams.values.max_window_radius));
set(handles.WindowStep,'String',num2str(curParams.values.nwindow_radii));
set(handles.WindowOffsets,'String',num2str(curParams.values.window_offsets));
set(handles.TransNone,'Value',...
  any(strcmp('none',curParams.values.trans_types)));
  %bitand(1,curParams.values.trans_types));
set(handles.TransFlip,'Value',...
  any(strcmp('flip',curParams.values.trans_types)));
  %bitand(4,curParams.values.trans_types));
set(handles.TransAbs,'Value',...
  any(strcmp('abs',curParams.values.trans_types)));
  %bitand(2,curParams.values.trans_types));
set(handles.TransRel,'Value',...
  any(strcmp('relative',curParams.values.trans_types)));
  %bitand(8,curParams.values.trans_types));
if isfield(curParams.values,handles.winextraParams{winNdx})
  extraParam = handles.winextraParams{winNdx};
  set(handles.ExtraParams,'String',curParams.values.(extraParam));
  set(handles.extraParamStatic,'String',extraParam);
else
  set(handles.ExtraParams,'Enable','off','String','--');
  set(handles.extraParamStatic,'String','Feature Specific');
end

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
      if ~isempty(handles.winextraParams{winfnNdx})
        extraParam = handles.winextraParams{winfnNdx};
        extraParamVal = curD.(curFn).values.(extraParam);
        params.(curPf).(curFn).(extraParam) = extraParamVal;
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

function docNode = createParamsXML(params,basicData,inFeatureWindowSize)
docNode = com.mathworks.xml.XMLUtils.createDocument('configuration');
toc = docNode.getDocumentElement;
basicStruct = struct;
for ndx = 1:size(basicData,1)
  curCat = basicData{ndx,1};
  basicStruct.(curCat).mode = basicData{ndx,2};
  basicStruct.(curCat).selection = basicData{ndx,3};
end
featureWindowSize = struct;
featureWindowSize.size = inFeatureWindowSize;
toc.appendChild(createXMLNode(docNode,'basicParams',basicStruct));
toc.appendChild(createXMLNode(docNode,'featureWindowSize',featureWindowSize));
toc.appendChild(createXMLNode(docNode,'params',params));

% att = fieldnames(params);
% for ndx = 1:numel(att)
%   toc.appendChild(createXMLNode(docNode,att{ndx},params.(att{ndx})));
% end


function disableWindowTable(handles)
set(handles.windowTable,'enable','off');
set(handles.pushbutton_copy_windowtypes,'enable','off');
disableWinParams(handles);

function enableWindowTable(handles)
set(handles.windowTable,'enable','on');
set(handles.pushbutton_copy_windowtypes,'enable','on');

function disableWinParams(handles)
set(handles.MinWindow,'enable','off');
set(handles.MaxWindow,'enable','off');
set(handles.WindowStep,'enable','off');
set(handles.WindowOffsets,'enable','off');
set(handles.TransNone,'enable','off');
set(handles.TransFlip,'enable','off');
set(handles.TransAbs,'enable','off');
set(handles.TransRel,'enable','off');
set(handles.ExtraParams,'enable','off');
set(handles.pushbutton_copy_windowparams,'enable','off');

function enableWinParams(handles)
defBack = get(handles.text5,'BackgroundColor');
defFore = [0 0 0];
set(handles.MinWindow,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.MaxWindow,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.WindowStep,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.WindowOffsets,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
if ~isempty(handles.winextraParams{handles.winNdx})
  set(handles.ExtraParams,'enable','on');
  set(handles.extraParamStatic,'String',handles.winextraParams{handles.winNdx});
else
  set(handles.ExtraParams,'enable','off');
  set(handles.extraParamStatic,'String','Feature Specific');
end
set(handles.TransNone,'enable','on');
set(handles.TransFlip,'enable','on');
set(handles.TransAbs,'enable','on');
set(handles.TransRel,'enable','on');
set(handles.pushbutton_copy_windowparams,'enable','on');

function MinWindow_Callback(hObject, eventdata, handles)
% hObject    handle to MinWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MinWindow as text
%        str2double(get(hObject,'String')) returns contents of MinWindow as a double
curVal = str2double(get(hObject,'String'));
if isnan(curVal)||(round(curVal)-curVal)~=0
  msgbox('Enter numerical values');
end
handles = guidata(hObject);
curFn = handles.windowComp{handles.winNdx};
handles.data{handles.pfNdx}.(curFn).values.min_window_radius = curVal;
guidata(hObject,handles);
setCategoryToCustom(handles);


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
if isnan(curVal)||(round(curVal)-curVal)~=0
  msgbox('Enter numerical values');
end
handles = guidata(hObject);
curFn = handles.windowComp{handles.winNdx};
handles.data{handles.pfNdx}.(curFn).values.max_window_radius = curVal;
guidata(hObject,handles);
setCategoryToCustom(handles);

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
if isnan(curVal) || (round(curVal)-curVal)~=0
  msgbox('Enter an integer value');
end
handles = guidata(hObject);
curFn = handles.windowComp{handles.winNdx};
handles.data{handles.pfNdx}.(curFn).values.nwindow_radii = curVal;
guidata(hObject,handles);
setCategoryToCustom(handles);


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
setCategoryToCustom(handles);


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

% --- Executes on button press in TransNone.
function TransNone_Callback(hObject, eventdata, handles)
% hObject    handle to TransNone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransNone
curFn = handles.windowComp{handles.winNdx};
handles = guidata(hObject);
curT = handles.data{handles.pfNdx}.(curFn).values.trans_types;
if get(hObject,'Value')
  %curT=bitor(1,curT);
  if ~any(strcmp('none',curT))
   handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'none';
  end
else
%   curT=bitand(14,curT);
  allNdx = strcmp('none',curT);
  handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
  if isempty(handles.data{handles.pfNdx}.(curFn).values.trans_types),
%  if handles.data{handles.pfNdx}.(curFn).values.trans_types==0,
    warndlg('Select at least one transformation type');
  end
end
guidata(hObject,handles);
setCategoryToCustom(handles);


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
%   curT=bitor(1,curT);
  if ~any(strcmp('flip',curT))
   handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'flip';
  end
else
%   curT=bitand(14,curT);
  allNdx = find(strcmp('flip',curT));
  handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
  if isempty(handles.data{handles.pfNdx}.(curFn).values.trans_types),
%   if handles.data{handles.pfNdx}.(curFn).values.trans_types==0,
    warndlg('Select at least one transformation type');
  end
end
guidata(hObject,handles);
setCategoryToCustom(handles);


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
%   curT=bitor(1,curT);
  if ~any(strcmp('abs',curT))
   handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'abs';
  end
else
%   curT=bitand(14,curT);
  allNdx = find(strcmp('abs',curT));
  handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
  if isempty(handles.data{handles.pfNdx}.(curFn).values.trans_types),
%   if handles.data{handles.pfNdx}.(curFn).values.trans_types==0,
    warndlg('Select at least one transformation type');
  end
end
guidata(hObject,handles);
setCategoryToCustom(handles);


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
%   curT=bitor(1,curT);
  if ~any(strcmp('relative',curT))
   handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'relative';
  end
else
%   curT=bitand(14,curT);
  allNdx = find(strcmp('relative',curT));
  handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
  if isempty(handles.data{handles.pfNdx}.(curFn).values.trans_types),
%   if handles.data{handles.pfNdx}.(curFn).values.trans_types==0,
    warndlg('Select at least one transformation type');
  end
end
guidata(hObject,handles);
setCategoryToCustom(handles);


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
basicData = get(handles.basicTable,'Data');
featureWindowSize = str2double(get(handles.editSize,'String'));
[params,cellParams] = convertData(handles);
set(handles.output,'Visible','off');
handles.JLDobj.UpdatePerframeParams(params,cellParams,basicData,featureWindowSize);



function ExtraParams_Callback(hObject, eventdata, handles)
% hObject    handle to ExtraParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ExtraParams as text
%        str2double(get(hObject,'String')) returns contents of ExtraParams as a double

if isempty(handles.winextraParams{handles.winNdx}), return; end

str = get(hObject,'String');
str = strtrim(str);
vals = str2num(str);
extraParam = handles.winextraParams{handles.winNdx};
curFn = handles.windowComp{handles.winNdx};
handles.data{handles.pfNdx}.(curFn).values.(extraParam) = vals;
guidata(hObject,handles);
setCategoryToCustom(handles);


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


% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Use project name.
[fName,pName] = uiputfile('params/*.xml','Save feature configurations to..');
if ~fName
  return;
end

[params,~] = convertData(handles);
basicData = get(handles.basicTable,'Data');
featureWindowSize = round(str2double(get(handles.editSize,'String')));
docNode = createParamsXML(params,basicData,featureWindowSize);
fName = fullfile(pName,fName);
xmlwrite(fName,docNode);


% --- Executes on button press in Load.
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
[fname,pname]= uigetfile('*.xml');
if ~fname; return; end;

featureparamsfilename = fullfile(pname,fname);
[params,~,basicTable,windowSize] = ...
  ReadPerFrameParams(featureparamsfilename,handles.JLDobj.featureConfigFile);
guidata(hObject,handles);
initData(hObject,params);
set(handles.basicTable,'Data',basicTable);
set(handles.editSize,'String',num2str(windowSize));


% --- Executes on selection change in popupmenu_copy_windowparams.
function popupmenu_copy_windowparams_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_copy_windowparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_copy_windowparams contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_copy_windowparams


% --- Executes during object creation, after setting all properties.
function popupmenu_copy_windowparams_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_copy_windowparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_copy_windowparams.
function pushbutton_copy_windowparams_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy_windowparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.pfNdx) || isempty(handles.winNdx),
  return;
end

windowComp = handles.windowComp;

% which perframe fn are we on?
pfNdx = handles.pfNdx;

% which window fn are we copying to?
winNdxTo = handles.winNdx;

% which windwo fn are we copying from?
winNdxFrom = get(handles.popupmenu_copy_windowparams,'Value');

% --- selected?
if winNdxTo > numel(windowComp),
  return;
end

handles = CopyWindowParams(handles,pfNdx,winNdxFrom,pfNdx,winNdxTo);
guidata(hObject,handles);
setCategoryToCustom(handles);


% --- Executes on selection change in popupmenu_copy_windowtypes.
function popupmenu_copy_windowtypes_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_copy_windowtypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_copy_windowtypes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_copy_windowtypes


% --- Executes during object creation, after setting all properties.
function popupmenu_copy_windowtypes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_copy_windowtypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_copy_windowtypes.
function pushbutton_copy_windowtypes_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy_windowtypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.pfNdx),
  return;
end

pfList = handles.pfList;
windowComp = handles.windowComp;

% which perframe fn are we copying to?
pfNdxTo = handles.pfNdx;
if ~handles.data{pfNdxTo}.valid,
  return;
end

% which perframe fn are we copying from?
pfNdxFrom = get(handles.popupmenu_copy_windowtypes,'Value');

% --- selected?
if pfNdxFrom > numel(pfList),
  return;
end

% if ~handles.data{pfNdxFrom}.valid,
%   return;
% end

for winfnNdx = 1:numel(windowComp),
  handles = CopyWindowParams(handles,pfNdxFrom,winfnNdx,pfNdxTo,winfnNdx);
end

setWindowTable(handles,pfNdxTo);

guidata(hObject,handles);
setCategoryToCustom(handles);


function handles = CopyWindowParams(handles,pfNdxFrom,winfnNdxFrom,pfNdxTo,winfnNdxTo)

curFnFrom = handles.windowComp{winfnNdxFrom};
% something to copy from?
if ~isfield(handles.data{pfNdxFrom},curFnFrom),
  return;
end
curFnTo = handles.windowComp{winfnNdxTo};
handles.data{pfNdxTo}.(curFnTo).valid = handles.data{pfNdxFrom}.(curFnFrom).valid;
for winParamsNdx = 1:numel(handles.winParams),
  curType = handles.winParams{winParamsNdx};
  handles.data{pfNdxTo}.(curFnTo).values.(curType) = ...
    handles.data{pfNdxFrom}.(curFnFrom).values.(curType);
end
if ~isempty(handles.winextraParams{winfnNdxTo})
  extraParam = handles.winextraParams{winfnNdxTo};
  if isfield(handles.data{pfNdxFrom}.(curFnFrom).values,extraParam)
    handles.data{pfNdxTo}.(curFnTo).values.(extraParam) = ...
      handles.data{pfNdxFrom}.(curFnFrom).values.(extraParam);
  else
    handles.data{pfNdxTo}.(curFnTo).values.(extraParam) = '';
  end
end
setWinParams(handles,winfnNdxTo);
if handles.data{pfNdxTo}.(curFnTo).valid,
  enableWinParams(handles);
end

function handles = CopyDefaultWindowParams(handles,category,pfNdxTo,winfnNdx)

curFn = handles.windowComp{winfnNdx};
% something to copy from?
if ~isfield(handles.categ.(category),curFn) &&...
    isfield(handles.data{pfNdxTo},curFn),
    handles.data{pfNdxTo}.(curFn).valid = false;
  return;
end
handles.data{pfNdxTo}.(curFn).valid = handles.categ.(category).(curFn).valid;
for winParamsNdx = 1:numel(handles.winParams),
  curType = handles.winParams{winParamsNdx};
  handles.data{pfNdxTo}.(curFn).values.(curType) = ...
    handles.categ.(category).(curFn).values.(curType);
end
if ~isempty(handles.winextraParams{winfnNdx})
  extraParam = handles.winextraParams{winfnNdx};
  handles.data{pfNdxTo}.(curFn).values.(extraParam) = ...
    handles.categ.(category).(curFn).values.(extraParam);
end

% --- Executes when selected object is changed in uipanel_tabs.
function uipanel_tabs_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel_tabs 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

if eventdata.NewValue == handles.togglebutton_tabdescription,
  handles.currentTab = 'description';
elseif eventdata.NewValue == handles.togglebutton_tabperframehistogram,
  handles.currentTab = 'perframehistogram';
end

handles = UpdateDescriptionPanels(handles);
guidata(hObject,handles);

function handles = UpdateDescriptionPanels(handles)

if isempty(handles.pfNdx),
  return;
end

% update visibility
if strcmpi(handles.currentTab,'description')
  set(handles.uipanel_description,'Visible','on');
  set(handles.uipanel_histogram,'Visible','off');
else
  set(handles.uipanel_description,'Visible','off');
  set(handles.uipanel_histogram,'Visible','on');
end

% histogram if necessary
if strcmpi(handles.currentTab,'perframehistogram'),
  if handles.histogramData.lastPfNdx ~= handles.pfNdx || ...
    ~strcmpi(handles.histogramData.lastType,'perframe'),
    
    i = find(handles.histogramData.perframe_idx == handles.pfNdx,1);
    if isempty(i),
      i = numel(handles.histogramData.perframe_idx)+1;
      [handles.histogramData.hhist,...
        ~,~,hleg,hxlabel,hylabel,...
        handles.histogramData.frac{i},handles.histogramData.frac_outside{i},...
        handles.histogramData.edges{i},handles.histogramData.centers_plot{i}] = ...
        HistogramPerFrameFeature(handles.JLDobj,handles.pfList{handles.pfNdx},...
        'axes',handles.axes_histogram,...
        'unknowncolor','w',...
        'labelcolors',jet(handles.JLDobj.nbehaviors)*.7);
      handles.histogramData.perframe_idx(i) = handles.pfNdx;
    else
      [handles.histogramData.hhist,...
        ~,~,hleg,hxlabel,hylabel,...
        handles.histogramData.frac{i},handles.histogramData.frac_outside{i},...
        handles.histogramData.edges{i},handles.histogramData.centers_plot{i}] = ...
        HistogramPerFrameFeature(handles.JLDobj,handles.pfList{handles.pfNdx},...
        'axes',handles.axes_histogram,...
        'edges',handles.histogramData.edges{i},...
        'frac',handles.histogramData.frac{i},...
        'frac_outside',handles.histogramData.frac_outside{i},...
        'unknowncolor','w',...
        'labelcolors',jet(handles.JLDobj.nbehaviors)*.7);
    end
    handles.histogramData.lastPfNdx = handles.pfNdx;
    handles.histogramData.lastType = 'perframe';
    textcolor = get(handles.togglebutton_tabdescription,'ForegroundColor');
    set(handles.axes_histogram,'XColor',textcolor,'YColor',textcolor,'Color','k','Box','off');
    set(hxlabel,'Color',textcolor);
    set(hylabel,'Color',textcolor);
    set(hleg,'Color','k','TextColor',textcolor,'Box','off','Location','Best');
  end
end


function editSize_Callback(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSize as text
%        str2double(get(hObject,'String')) returns contents of editSize as a double
curVal = str2double(get(hObject,'String'));
if isempty(curVal)||(round(curVal)-curVal)~=0
  msgbox('Enter numerical values. eg: "-1 0 1" (without with quotes)');
  return;
end

% First set the default category values

categories = fieldnames(handles.categ);
winComp = handles.windowComp;

for cndx = 1:numel(categories)
  curCat = categories{cndx};
  for wndx = 1:numel(winComp)
    if ~isfield(handles.categ.(curCat),winComp{wndx}); continue; end
    handles.categ.(curCat).(winComp{wndx}).values.max_window_radius = curVal;
  end
end

% Now copy the default values to the perframe features.
basicData = get(handles.basicTable,'Data');
for ndx = 1:size(basicData,1)
  if ~strcmp(basicData{ndx,2},'All'); continue;end
  handles = applyCategoryType(handles,ndx);
end
guidata(hObject,handles);
createPfTable(handles.pfTable);


% --- Executes during object creation, after setting all properties.
function editSize_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in togglebuttonMode.
function togglebuttonMode_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
curLoc = get(handles.figure1,'Position');

if get(hObject,'Value')
  handles.mode = 'advanced';
  set(handles.figure1,'Position',[curLoc(1:2) handles.advancedSize]);
  set(hObject,'String','Basic <');
else
  handles.mode = 'basic';
  set(handles.figure1,'Position',[curLoc(1:2) handles.basicSize]);  
  set(hObject,'String','Advanced >');
end


% --- Executes on button press in pushbutton_hist.
function pushbutton_hist_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prcEdges = [5 15 30 50 70 85 95];

histfnNdx = find(strcmp('hist',handles.windowComp));
histExtraName = handles.winextraParams{histfnNdx};

for ndx = 1:numel(handles.pfList)
  curPf = handles.pfList{ndx};
  
  if ~handles.data{ndx}.valid || ~handles.data{ndx}.hist.valid; 
    continue;
  end
  
  allData = [];
  for expi = 1:handles.JLDobj.nexps,
    
    % load per-frame data for this experiment
    perframedir = handles.JLDobj.GetFile('perframedir',expi);
    file = fullfile(perframedir,[curPf,'.mat']);
    if ~exist(file,'file'),
      warning('Per-frame data file %s does not exist',file);
      continue;
    end
    
    perframedata = load(file);
    
    for fly = 1:handles.JLDobj.nflies_per_exp(expi),
      
      x = perframedata.data{fly};
      allData = [allData ; x(:)];
    end
  end
  bins = prctile(allData,prcEdges);
  minD = min(allData);
  maxD = max(allData);
  binMin = minD - (maxD-minD);
  binMax = maxD + (maxD-minD);
  bins = [binMin bins binMax];
  handles.data{ndx}.hist.values.(histExtraName) = bins;
end
guidata(hObject,handles);
if ~isempty(handles.pfNdx)
  setWindowTable(handles,handles.pfNdx);
end


function CloseRequestFcn(hObject,eventdata,handles)
push_cancel_callback(hObject,eventdata,handles);


% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

configfile = handles.JLDobj.configfilename;
configparams = ReadXMLParams(configfile);

if ~isfield(configparams.file,'featureparamfilename') || isempty(configparams.file.featureparamfilename)
  behaviorname = configparams.behaviors.names;
  if iscell(behaviorname),
    defaultname = sprintf('WindowFeatures_%s.xml',behaviorname{1});
  else
    defaultname = sprintf('WindowFeatures_%s.xml',behaviorname);
  end
  [fname,fpath]= uiputfile(fullfile('params','*.xml'),'Enter a name for feature config file',defaultname);
  if isempty(fname),return, end
  featureconfigfile = fullfile(fpath,fname);
  configparams.file.featureparamfilename = featureconfigfile;
  docNode = com.mathworks.xml.XMLUtils.createDocument('params');
  toc = docNode.getDocumentElement;
  fnames = fieldnames(configparams);
  for ndx = 1:numel(fnames)
    toc.appendChild(createXMLNode(docNode,fnames{ndx},configparams.(fnames{ndx})));
  end
  xmlwrite(configfile,docNode);
end

featureconfigfile = configparams.file.featureparamfilename;
[params,~] = convertData(handles);
basicData = get(handles.basicTable,'Data');
featureWindowSize = round(str2double(get(handles.editSize,'String')));
docNode = createParamsXML(params,basicData,featureWindowSize);
xmlwrite(featureconfigfile,docNode);

pushbutton_done_Callback(hObject,eventdata,handles);
push_cancel_Callback(hObject,eventdata,handles);
