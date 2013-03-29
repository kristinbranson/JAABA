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


% --- Executes just before SelectFeatures is made visible.-----------------
function SelectFeatures_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectFeatures (see VARARGIN)

% Change a few things so they still work well on Mac
adjustColorsIfMac(hObject);

% Choose default command line output for SelectFeatures
figureJLabel = varargin{1};
handles.figureJLabel=figureJLabel;
jld = JLabel('getJLabelData',figureJLabel);
handles.output = hObject;

handles.mode = 'basic';
curPos = get(handles.figure1,'Position');
tablePos = get(handles.basicTable,'Position');
reducedWidth = tablePos(1)+tablePos(3) + 15;
reducedHeight = curPos(4);
handles.advancedSize = curPos(3:4);
handles.basicSize = [reducedWidth reducedHeight];
set(handles.figure1,'Position',[curPos(1:2) handles.basicSize]);

%handles.scoresasinput = varargin{2};
%handles.scoresasinput = JLDobj.scoresasinput;

% create the object that will serve as the model, or at least part of the
% model
handles.featureVocabulary=FeatureVocabularyForSelectFeatures(jld);
handles.jld=jld;

% Update handles structure
guidata(hObject, handles);

set(handles.pfTable,'UserData',0);
setJLDobj(hObject,jld);
set(hObject,'Visible','on');  % have to do this b/c of Java hacking
removeRowHeaders(hObject);
pause(0.5);

return


% --- Outputs from this function are returned to the command line.
function varargout = SelectFeatures_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
return


% -------------------------------------------------------------------------
function setJLDobj(hObject,JLDobj)
% Set the JLabelDataObject.
handles = guidata(hObject);
% handles.JLDobj = JLDobj;

% handles.windowComp = {'default','mean','min','max','hist','prctile',...
%    'change','std','harmonic','diff_neighbor_mean',...
%    'diff_neighbor_min','diff_neighbor_max','zscore_neighbors'};
% 
% handles.winextraParams = {'','','','','hist_edges','prctiles','change_window_radii',...
%   '','num_harmonic','','','',''};
% 
% handles.winParams = {'max_window_radius','min_window_radius','nwindow_radii',...
%   'trans_types','window_offsets'};
% 
% handles.winextraDefaultParams = {[],[],[],[],[-400000 0 40000],[5 10 30 50 70 90 95],[1 3],...
%   [],[2],[],[],[],[]};
% handles.defaultWinParams = {10,1,3,{'none'},0};

guidata(hObject,handles);

windowFeatureParams = JLDobj.GetPerframeParams();
featureLexicon= JLDobj.featureLexicon;
initData(hObject,windowFeatureParams,featureLexicon,JLDobj);

return


% -------------------------------------------------------------------------
function createWindowTable(hObject)
% Sets values for the window table.

handles = guidata(hObject);

% Deal with windowTable
featureVocabulary=handles.featureVocabulary;
wfTypes=featureVocabulary.wfTypes;
set(handles.windowTable,'RowName', wfTypes,...
    'ColumnName',{'Computation Type','Select'});
selVals = [true true true true false false false false false false false false false];
tableData = {};
for ndx = 1:numel(wfTypes)
  tableData{ndx,1} = wfTypes{ndx};
  tableData{ndx,2} = selVals(ndx);
end
set(handles.windowTable,'Data',tableData);
set(handles.windowTable','ColumnWidth',{135 'auto'});
set(handles.windowTable,'ColumnEditable',[false,true]);
set(handles.windowTable,'CellSelectionCallback',@windowSelect);
set(handles.windowTable,'CellEditCallback',@windowEdit);

guidata(hObject,handles);
return


% -------------------------------------------------------------------------
function createPfTable(hObject)
% Sets the values for feature table.

handles = guidata(hObject);

featureVocabulary=handles.featureVocabulary;  % a ref
pfList = featureVocabulary.pfNameList;
tableData = cell(length(pfList),2);
for ndx = 1:numel(pfList)
  tableData{ndx,1} = pfList{ndx};
  tableData{ndx,2} = featureVocabulary.vocabulary{ndx}.valid;
  pfCategories=featureVocabulary.pfCategoriesFromName.(pfList{ndx});
  str = sprintf('%s',pfCategories{1});
  for sndx = 2:numel(pfCategories)
    str = sprintf('%s,%s',str,pfCategories{sndx});
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


% -------------------------------------------------------------------------
function createFeatureTable(hObject,JLDobj)
handles = guidata(hObject);

% Adding drop down menu kind of thing.

pfCategoryList=handles.featureVocabulary.pfCategoryList;
if ~isempty(JLDobj.basicFeatureTable)
  oldtableData = JLDobj.basicFeatureTable;
  tableData = {};
  for ndx = 1:numel(pfCategoryList)
    curcateg = pfCategoryList{ndx};
    tableData{ndx,1} = curcateg;
    ondx = find(strcmpi(oldtableData(:,1),curcateg));
    if isempty(ondx),
      tableData{ndx,2} = 'none';
      tableData{ndx,3} = 'normal';
    else
      tableData(ndx,2:3) = oldtableData(ondx,2:3);
    end
  end  
else
  tableData = {};
  for ndx = 1:numel(pfCategoryList)
    tableData{ndx,1} = pfCategoryList{ndx};
    tableData{ndx,2} = 'none';
    tableData{ndx,3} = 'normal';
  end
end
set(handles.basicTable,'Data',tableData);
set(handles.basicTable,'ColumnName',{'Categories','Select','Amount'});
set(handles.basicTable,'ColumnFormat',{'char',...
  {'all' 'none' 'custom'}, handles.featureVocabulary.wfAmounts'});
set(handles.basicTable,'ColumnEditable',[false,true,true]);

set(handles.basicTable,'ColumnWidth',{85,65,75});
set(handles.basicTable,'CellSelectionCallback',@basicSelect);
set(handles.basicTable,'CellEditCallback',@basicEdit);


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function createCopyFromMenus(hObject)

handles = guidata(hObject);
% can copy from any window feature type
fv=handles.featureVocabulary;
wfTypes=fv.wfTypes;
pfNameList=fv.pfNameList;
set(handles.popupmenu_copy_windowparams, ...
    'String',[wfTypes,{'---'}],...
    'Value',numel(wfTypes)+1);
% can copy from any per-frame feature
set(handles.popupmenu_copy_windowtypes, ...
    'String',[pfNameList;{'---'}],...
    'Value',numel(pfNameList)+1);
return


% -------------------------------------------------------------------------
function createDescriptionPanels(hObject)

handles = guidata(hObject);

% which tab to show
handles.currentTab = 'perframehistogram';
set(handles.togglebutton_tabperframehistogram,'Value',1);

% set visibility
uipanel_tabs_SelectionChangeFcn(handles.uipanel_tabs, struct('NewValue',handles.togglebutton_tabdescription), handles);

guidata(hObject,handles);


% -------------------------------------------------------------------------
% Initialize the data structure.
function initData(hObject,params,featureLexicon,JLDobj)
readFeatureLexicon(hObject,featureLexicon);
createFeatureTable(hObject,JLDobj);

handles = guidata(hObject);
if ~isempty(JLDobj.featureWindowSize)
  set(handles.editSize,'String',num2str(JLDobj.featureWindowSize));
  
%   curVal = handles.JLDobj.featureWindowSize;
%   wfAmounts = fieldnames(handles.wfParamsFromAmount);
%   winComp = handles.windowComp;
% 
%   for cndx = 1:numel(wfAmounts)
%     curCat = wfAmounts{cndx};
%     for wndx = 1:numel(winComp)
%       if ~isfield(handles.wfParamsFromAmount.(curCat),winComp{wndx}); continue; end
%       handles.wfParamsFromAmount.(curCat).(winComp{wndx}).values.max_window_radius = curVal;
%     end
%   end

end

% pfList = fieldnames(handles.pfCategoriesFromName);
% handles.pfList = pfList;
% data = {};
% validPfs = fieldnames(params);
% winParams = handles.winParams;
% for ndx = 1:numel(pfList)
%   curPfName = pfList{ndx};
%   data{ndx}.name = curPfName;
%   pNdx = find(strcmp(curPfName,validPfs));
%   if pNdx
%     data{ndx}.valid = true;
%     data{ndx}.sanitycheck = params.(curPfName).sanitycheck;
%     
%     % Fill the default values.
%     for winParamsNdx = 1:numel(winParams)
%       curType = winParams{winParamsNdx};
%       if isfield(params.(curPfName),curType)
%         data{ndx}.default.values.(curType) = params.(curPfName).(curType);
%       else % Fill in the default value
%         data{pfNdx}.default.values.(curType) = ...
%             handles.defaultWinParams{winParamsNdx};
%       end
%       data{ndx}.default.valid = true;
%     end
%     
%     % Fill for different window function type.
%     
%     % Find which ones are valid.
%     curParam = params.(curPfName);
%     validWinfn = fieldnames(curParam);
%     for winfnNdx = 2:numel(handles.windowComp)
%       curFn = handles.windowComp{winfnNdx};
%       wNdx = find(strcmp(curFn,validWinfn));
%       
%       if wNdx,
%         
%         data{ndx}.(curFn).valid = true;
%         curWinFnParams = params.(curPfName).(curFn);
%         for winParamsNdx = 1:numel(winParams)
%           curType = winParams{winParamsNdx};
%           if isfield(curWinFnParams,curType)
%             data{ndx}.(curFn).values.(curType) = curWinFnParams.(curType);
%           else % fill in the default values
%             data{ndx}.(curFn).values.(curType) = data{ndx}.default.values.(curType);
%           end
%         end
%         
%         if ~isempty(handles.winextraParams{winfnNdx})
%           extraParam = handles.winextraParams{winfnNdx};
%           data{ndx}.(curFn).values.(extraParam) = curWinFnParams.(extraParam);
%         end
%         
%       else % Values for window comp haven't been defined.
%         
%         data{ndx}.(curFn).valid = false;
%         for winParamsNdx = 1:numel(winParams)
%           curType = winParams{winParamsNdx};
%           data{ndx}.(curFn).values.(curType) = data{ndx}.default.values.(curType);
%         end
%         if ~isempty(handles.winextraParams{winfnNdx})
%           extraParam = handles.winextraParams{winfnNdx};
%           data{ndx}.(curFn).values.(extraParam) = handles.winextraDefaultParams{winfnNdx};
%         end
%         
%       end
%       
%     end
%   
%   else % Default values for invalid pf's.
%     
%     data{ndx}.valid = false;
%     data{ndx}.sanitycheck = false;
% 
%     data{ndx}.default.valid = true;
%     for winParamsNdx = 1:numel(handles.winParams)
%       curType = handles.winParams{winParamsNdx};
%       data{ndx}.default.values.(curType) = ...
%         handles.defaultWinParams{winParamsNdx};
%     end
%     
%     % Copy the default values into the other window params.
%     for winfnNdx = 2:numel(handles.windowComp)
%       curFn = handles.windowComp{winfnNdx};
%       data{ndx}.(curFn).valid = false;
%       for winParamsNdx = 1:numel(handles.winParams)
%         curType = handles.winParams{winParamsNdx};
%         data{ndx}.(curFn).values.(curType) = ...
%           data{ndx}.default.values.(curType);
%       end
%         if ~isempty(handles.winextraParams{winfnNdx})
%           extraParam = handles.winextraParams{winfnNdx};
%           data{ndx}.(curFn).values.(extraParam) = handles.winextraDefaultParams{winfnNdx};
%         end
%     end
% 
%     
%   end
% end

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

% handles.data = data;
handles.pfNdx = [];
handles.winNdx = [];
guidata(hObject,handles);
createPfTable(hObject);
createWindowTable(hObject);
createCopyFromMenus(hObject);
createDescriptionPanels(hObject);
compatibleBasicAdvanced(handles);
return


% -------------------------------------------------------------------------
function readFeatureLexicon(hObject,featureLexicon)
% Reads the configuration file that sets the default parameters.

handles = guidata(hObject);

%configfile = handles.JLDobj.featureConfigFile;
%settings = ReadXMLParams(configfile);
settings=featureLexicon;

% Read the default parameters for different categories.
categories = fieldnames(settings.defaults);

% if isempty(categories) 
%   % If no default have been specified in the config file.
%   
%   categories{1} = 'default';
%   curParams = struct;
%   
%   % Default values.
%   for ndx = 1:numel(handles.winParams)
%     curwinpname = handles.winParams{ndx};
%     curParams.default.values.(curwinpname) = handles.defaultWinParams{ndx};
%   end
%   curParams.default.valid = true;
%   curParams.default.values.sanitycheck = false;
%   
%   % For each defined window computation.
%   for ndx = 2:numel(handles.windowComp)
%     % Copy the default values.
%     curwname = handles.windowComp{ndx};
%     for pndx = 1:numel(handles.winParams)
%       curwinpname = handles.winParams{pndx};
%       curParams.(curwname).values.(curwinpname) = handles.defaultWinParams{pndx};
%     end
%     
%     if ~isempty(handles.winextraParams{ndx})
%       extraParam = handles.winextraParams{ndx};
%       curParams.(curwname).values.(extraParam) = handles.winextraParams{ndx};
%     end
%     curParams.(curwname).valid = false;
%   end
%   handles.wfParamsFromAmount.(categories{1}) = curParams;
%   
%   
% else
%   % fill the window params for different feature categories.
%   
%   for j = 1:numel(categories)
%     curParams = struct;
%     cur = settings.defaults.(categories{j});
%     
%     % Default values.
%     for ndx = 1:numel(handles.winParams)
%       curwinpname = handles.winParams{ndx};
%       curParams.default.values.(curwinpname) = cur.(curwinpname);
%     end
%     curParams.default.valid = true;
%     curParams.default.values.sanitycheck = false;
%     
%     % For each defined window computation.
%     for ndx = 2:numel(handles.windowComp)
%       % Copy the default values.
%       curwname = handles.windowComp{ndx};
%       for pndx = 1:numel(handles.winParams)
%         curwinpname = handles.winParams{pndx};
%         curParams.(curwname).values.(curwinpname) = curParams.default.values.(curwinpname);
%       end
%       
%       % Override the default window params for the current window computation type
%       if isfield(cur,curwname)
%         curParams.(curwname).valid = true;
%         diffFields = fieldnames(cur.(curwname));
%         for dndx = 1:numel(diffFields)
%           if any(strcmp(handles.winParams,diffFields{dndx}))
%             curwinpname = diffFields{dndx};
%             curParams.(curwname).values.(curwinpname) = cur.(curwname).(curwinpname);
%           end
%         end
%         
%         if ~isempty(handles.winextraParams{ndx})
%           extraParam = handles.winextraParams{ndx};
%           if isfield(cur.(curwname),extraParam)
%             curParams.(curwname).values.(extraParam) = cur.(curwname).(extraParam);
%           else
%             curParams.(curwname).values.(extraParam) = '';
%           end
%         end
%       else
%         curParams.(curwname).valid = false;
%       end
%     end
%     handles.wfParamsFromAmount.(categories{j}) = curParams;
%   end
% end
% 
% % Now for each different perframe feature read the type default trans_types.
% perframeL = handles.JLDobj.allperframefns;
% scores_perframe = {handles.scoresasinput(:).scorefilename};
% transType = struct;
% pfCategoriesFromName = struct;
% 
% % perframeL might not contain all the perframe features.
% for pfndx = 1:numel(perframeL)
%   curpf = perframeL{pfndx};
%   
%   if any(strcmp(curpf ,scores_perframe)), % This is a score perframe.
%     transType.(curpf) = {'none'};
%     pfCategoriesFromName.(curpf) = {'scores'};
%     continue;
%   end
%   
%   transType.(curpf) = settings.perframe.(curpf).trans_types;
%   if ischar(transType.(curpf))
%     transType.(curpf) = {transType.(curpf)};
%   end
%   curtypes  = settings.perframe.(curpf).type; 
%   if ischar(curtypes)
%     pfCategoriesFromName.(curpf)  = {curtypes}; 
%   else    
%     pfCategoriesFromName.(curpf)  = curtypes; 
%   end
% end
% 
% fallpf = fieldnames(settings.perframe);
% pfCategoryList = {};
% for pfndx = 1:numel(fallpf)
%   curpf = fallpf{pfndx};
%   curtypes  = settings.perframe.(curpf).type; 
%   if ischar(curtypes)
%     curT = curtypes;
%     if ~any(strcmp(pfCategoryList,curT))
%       pfCategoryList{end+1} = curT;
%     end
%   else    
%     for tndx = 1:numel(curtypes)
%       curT = curtypes{tndx};
%       if ~any(strcmp(pfCategoryList,curT))
%         pfCategoryList{end+1} = curT;
%       end
%     end
%   end
% end
% if ~any(strcmp(pfCategoryList,'scores')),
%   pfCategoryList{end+1} = 'scores';
% end
% 
% handles.transType = transType;
% handles.pfCategoriesFromName = pfCategoriesFromName;
% handles.pfCategoryList = pfCategoryList;

% basicData = get(handles.basicTable,'Data');
% scoresBasicNdx = find(strcmpi(basicData{:,1},'scores'));
% handles = applyCategoryType(handles,scoresBasicNdx);

guidata(hObject,handles);
featureVocabulary=handles.featureVocabulary;  % a ref
wfType=featureVocabulary.wfTypes{1};
set(handles.editSize,'String',...
  num2str(featureVocabulary.wfParamsFromAmount.(categories{1}).(wfType).values.max_window_radius));
return


% -------------------------------------------------------------------------
function basicSelect(hObject,eventData)  %#ok
return


% -------------------------------------------------------------------------
function basicEdit(hObject,eventData)

% the user selects the category
if isempty(eventData.Indices); return; end

handles = guidata(hObject);
fv=handles.featureVocabulary;  % a ref
if eventData.Indices(2)==2
  % Select-deselect the whole category.  
  switch eventData.NewData
    case 'none'
      h = waitbar(0,'Unselecting the perframe features');
      basicData = get(handles.basicTable,'Data');
      pfNameList=fv.pfNameList;
      for ndx = 1:numel(pfNameList)
        pfName=pfNameList{ndx};
        pfCategories=fv.pfCategoriesFromName.(pfName);
        disable = true;
        for tndx = 1:size(basicData,1)
          if ~any(strcmp(pfCategories,basicData{tndx,1}));
            continue;
          end
          if ~strcmp(basicData{tndx,2},'none'); 
            disable = false; break; 
          end
        end
        if disable
          %handles.data{ndx}.valid = false;
          fv.disablePerframeFeature(pfName);
        end
      end
      close(h);
    case 'all'
      h = waitbar(0,'Selecting the perframe features');
      categoryIndex=eventData.Indices(1);
      %handles = applyCategoryType(handles,categoryIndex);
      basicData = get(handles.basicTable,'Data');
      handles.featureVocabulary.setPFCategoryToWFAmount(categoryIndex,basicData{categoryIndex,3});
      close(h);
  end
  guidata(hObject,handles);
  createPfTable(handles.pfTable);
elseif eventData.Indices(2) == 3
  % Set the window-feature level appropriately
  basicData = get(handles.basicTable,'Data');
  pfLevel=basicData{eventData.Indices(1),2};
  if isequal(pfLevel,'all')
    categoryIndex=eventData.Indices(1);
    %handles=applyCategoryType(handles,categoryIndex);
    handles.featureVocabulary.setPFCategoryToWFAmount(categoryIndex,basicData{categoryIndex,3});
    guidata(hObject,handles);
    createPfTable(handles.pfTable);
  end
end
return


% % -------------------------------------------------------------------------
% function handles = applyCategoryType(handles,iSelectedCategory)
% % get some variables out of handles
% basicData = get(handles.basicTable,'Data');
% wfAmount = basicData{iSelectedCategory,3};  
%   % The amount of window features to use.  Can be 'normal', 'more', or
%   % 'less'
% fv=handles.featureVocabulary;
% selectedCategory = handles.pfCategoryList{iSelectedCategory};
% pfList=handles.pfList;
% pfCategoriesFromName=handles.pfCategoriesFromName;
% % Copy the parameters.
% % Iterate over the per-frame features, looking for ones that are
% % within the selected category.
% for iPF = 1:numel(pfList)
%   thisPF=pfList{iPF};
%   categoriesThisPFIsIn=pfCategoriesFromName.(thisPF);
%   if ismember(selectedCategory,categoriesThisPFIsIn)
%     % if the selected category contains this per-frame feature,
%     % do something...
%     handles.data{iPF}.valid = true;
%     for winfnNdx = 1:numel(fv.wfTypes)
%       curFn = fv.wfTypes{winfnNdx};
%       if ~handles.wfParamsFromAmount.(wfAmount).(curFn).valid
%         handles.data{iPF}.(curFn).valid = false;
%         continue;
%       end
%       handles = CopyDefaultWindowParams(handles,...
%                                         wfAmount, iPF, winfnNdx);
%       curFn = fv.wfTypes{winfnNdx};
%       handles.data{iPF}.(curFn).values.trans_types = handles.transType.(pfList{iPF});
%     end
%   end
% end


% -------------------------------------------------------------------------
function compatibleBasicAdvanced(handles)
basicTable = get(handles.basicTable,'Data');
incompatible = '';
fv=handles.featureVocabulary;
pfNameList=fv.pfNameList;
for ndx = 1:numel(pfNameList)
  pfName=pfNameList{ndx};
  pfCategories=fv.pfCategoriesFromName.(pfName);
  selected = false;
  for bndx = 1:size(basicTable,1)
    if ~any(strcmp(pfCategories,basicTable{bndx,1})),
      continue;
    end
    if strcmpi(basicTable{bndx,2},'all')
       selected = true;
    end
    if selected && ~fv.pfIsInVocabulary(ndx)
      incompatible = sprintf('%s %s', incompatible,pfNameList{ndx});
    end
    
  end
end
if numel(incompatible)>0
  uiwait(warndlg(sprintf('Perframe feature(s) %s are not selected but should have been based on the category selection',incompatible),...
    'Mismatch in categories and perframe features'));
end


% -------------------------------------------------------------------------
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
    updateWindowTable(handles);
    handles.winNdx = 1;
    winData = get(handles.windowTable,'Data');

    if winData{handles.winNdx,2},
      updateWinParams(handles);
      enableWinParams(handles);
    else
      disableWinParams(handles);
    end

    updateWinParams(handles);
  else
    disableWindowTable(handles);
  end
end
handles = UpdateDescriptionPanels(handles);

guidata(hObject,handles);
% set(hObject,'UserData',0);


% -------------------------------------------------------------------------
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

%handles.data{pfNdx}.valid = eventData.NewData;

% Update the feature vocabulary
fv=handles.featureVocabulary;  % a ref
% Get the name of the per-frame feature
pfData=get(handles.pfTable,'Data');
pfName=pfData{pfNdx,1};
if eventData.NewData
  % Get an amount to set the PF to, based on the Amount displayed in the
  % basic table for one of the PFs categories
  pfCategoryNames=fv.pfCategoriesFromName.(pfName);
  pfCategoryName=pfCategoryNames{1};  
    % Just take the first possible category.  This may cause the other
    % categories (if there are any) to have their Amount go to Custom, but
    % that's OK
  pfCategoryIndex = find(strcmp(pfCategoryName,fv.pfCategoryList));
  basicData = get(handles.basicTable,'Data');
  wfAmount = basicData{pfCategoryIndex,3};  %#ok
  fv.enablePerframeFeature(pfName,wfAmount);
  enableWindowTable(handles);  
else
  fv.disablePerframeFeature(pfName);
  disableWindowTable(handles);  
end


%guidata(hObject,handles);

% % When a feature is unchecked, disable and return
% if ~eventData.NewData, 
%   %handles.data{pfNdx}.valid = false;
%   disableWindowTable(handles);  
%   guidata(hObject,handles);
%   return;
% end
% 
% handles.data{pfNdx}.valid = true;
% curType = handles.pfCategoriesFromName.(handles.pfList{pfNdx}){1};
% pfCategoryIndex = find(strcmp(handles.pfCategoryList,curType));
% basicData = get(handles.basicTable,'Data');
% handles.data{pfNdx}.valid = true;
% for winfnNdx = 1:numel(fv.wfTypes)
%   wfAmount = basicData{pfCategoryIndex,3};  %#ok
%   curFn = fv.wfTypes{winfnNdx};
%   if ~handles.wfParamsFromAmount.(wfAmount).(curFn).valid
%     handles.data{pfNdx}.(curFn).valid = false;
%     continue;
%   end
%   handles = CopyDefaultWindowParams(handles,...
%     wfAmount, pfNdx,winfnNdx);
%   curFn = fv.wfTypes{winfnNdx};
%   handles.data{pfNdx}.(curFn).values.trans_types = handles.transType.(handles.pfList{pfNdx});
% end


% 
% % When it already has values
% if isfield(handles.data{pfNdx},'default')
%   guidata(hObject,handles);
%   updateWindowTable(handles);
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
% for winfnNdx = 2:numel(fv.wfTypes)
%   curFn = fv.wfTypes{winfnNdx};
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
updatePFCategoryAmountForCurrentPF(handles);
updateWindowTable(handles);
% set(hObject,'UserData',0);
return


% -------------------------------------------------------------------------
function updateWindowTable(handles)
% Sets the data in window table based on selections made in pfTable.
pfNdx=handles.pfNdx;
enableWindowTable(handles);
fv=handles.featureVocabulary;

%curData = handles.data{pfNdx};
curData = fv.vocabulary{pfNdx};

windowData = {};
for ndx = 1:length(fv.wfTypes)
  curFunc = fv.wfTypes{ndx};
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
return


% -------------------------------------------------------------------------
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
    updateWinParams(handles);
    enableWinParams(handles);
  else
    disableWinParams(handles);
  end
end

guidata(hObject,handles);


% -------------------------------------------------------------------------
function windowEdit(hObject,eventData)
% Called when a window function is edited.

handles = guidata(hObject);
fv=handles.featureVocabulary;
updatePFCategoryAmountForCurrentPF(handles);
winNdx = eventData.Indices(1);
handles.winNdx = winNdx;
wfType = fv.wfTypes{winNdx};
%handles.data{handles.pfNdx}.(wfType).valid = eventData.NewData;

% Enable/disable the selected window-feature type in the model
pfNdx=handles.pfNdx;
pfName=handles.featureVocabulary.pfNameList{pfNdx};
winNdx = eventData.Indices(1);
handles.winNdx = winNdx;
wfType = fv.wfTypes{winNdx};
handles.featureVocabulary.setWFTypeEnablement(pfName,wfType,eventData.NewData)

guidata(hObject,handles);
if eventData.NewData
  updateWinParams(handles);
  enableWinParams(handles);
else
  disableWinParams(handles);
end
return


% -------------------------------------------------------------------------
function updatePFCategoryAmountForCurrentPF(handles)
fv=handles.featureVocabulary;
pfNameList=fv.pfNameList;
pfName=pfNameList{handles.pfNdx};
pfCategoriesThisName=fv.pfCategoriesFromName.(pfName);
allPFCategories=fv.pfCategoryList;
pfCategoryIndicesThisName = find(ismember(allPFCategories,pfCategoriesThisName));
basicData = get(handles.basicTable,'Data');
for ndx = 1:numel(pfCategoryIndicesThisName)
  basicData{pfCategoryIndicesThisName(ndx),2} = ...
    fv.getPFCategoryAmount(pfCategoriesThisName{ndx});
end
set(handles.basicTable,'Data',basicData);
return


% -------------------------------------------------------------------------
function updateWinParams(handles)
% Update the GUI controls displaying window parameters to match the current
% model state.

winNdx=handles.winNdx;
fv=handles.featureVocabulary;
curFn = fv.wfTypes{winNdx};
curParams = fv.vocabulary{handles.pfNdx}.(curFn);
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
if isfield(curParams.values,fv.wfExtraParamNames{winNdx})
  extraParam = fv.wfExtraParamNames{winNdx};
  set(handles.ExtraParams,'String',curParams.values.(extraParam));
  set(handles.extraParamStatic,'String',extraParam);
else
  set(handles.ExtraParams,'Enable','off','String','--');
  set(handles.extraParamStatic,'String','Feature Specific');
end


% % -------------------------------------------------------------------------
% function windowFeatureParams = convertData(handles)
% % Converts the data into format used by JLabelData.
% 
% fv=handles.featureVocabulary;
% 
% windowFeatureParams = struct;
% 
% for ndx = 1:numel(handles.pfList)
%   curPf = handles.pfList{ndx};
%   if ~fv.vocabulary{ndx}.valid; continue;end
%   curD = fv.vocabulary{ndx};
%   windowFeatureParams.(curPf).sanitycheck = curD.sanitycheck;
%   
%   for winParamsNdx = 1:numel(handles.winParams)
%     curType = handles.winParams{winParamsNdx};
%     windowFeatureParams.(curPf).(curType) = curD.default.values.(curType);
%   end
%   
%   for winfnNdx = 2:numel(fv.wfTypes)
%     curFn = fv.wfTypes{winfnNdx};
%     
%     if curD.(curFn).valid,
%       for winParamsNdx = 1:numel(handles.winParams)
%         curType = handles.winParams{winParamsNdx};
%         windowFeatureParams.(curPf).(curFn).(curType) =...
%             curD.(curFn).values.(curType);
%       end
%       if ~isempty(fv.wfExtraParamNames{winfnNdx})
%         extraParam = fv.wfExtraParamNames{winfnNdx};
%         extraParamVal = curD.(curFn).values.(extraParam);
%         windowFeatureParams.(curPf).(curFn).(extraParam) = extraParamVal;
%       end
%     end
%     
%   end
% 
% end
% 
% return


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function disableWindowTable(handles)
set(handles.windowTable,'enable','off');
set(handles.pushbutton_copy_windowtypes,'enable','off');
disableWinParams(handles);


% -------------------------------------------------------------------------
function enableWindowTable(handles)
set(handles.windowTable,'enable','on');
set(handles.pushbutton_copy_windowtypes,'enable','on');


% -------------------------------------------------------------------------
function updateWinParamsEnablement(handles)
fv=handles.featureVocabulary;   % a ref
if fv.wfTypeIsInVocabulary(handles.pfNdx,handles.winNdx),
  enableWinParams(handles);
end
return


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function enableWinParams(handles)
fv=handles.featureVocabulary;   % a ref
defBack = get(handles.text5,'BackgroundColor');
defFore = [0 0 0];
set(handles.MinWindow,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.MaxWindow,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.WindowStep,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.WindowOffsets,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
if ~isempty(fv.wfExtraParamNames{handles.winNdx})
  set(handles.ExtraParams,'enable','on');
  set(handles.extraParamStatic,'String',fv.wfExtraParamNames{handles.winNdx});
else
  set(handles.ExtraParams,'enable','off');
  set(handles.extraParamStatic,'String','Feature Specific');
end
set(handles.TransNone,'enable','on');
set(handles.TransFlip,'enable','on');
set(handles.TransAbs,'enable','on');
set(handles.TransRel,'enable','on');
set(handles.pushbutton_copy_windowparams,'enable','on');


% -------------------------------------------------------------------------
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
fv=handles.featureVocabulary;
wfType = fv.wfTypes{handles.winNdx};
%handles.data{handles.pfNdx}.(wfType).values.min_window_radius = curVal;

% set in featureVocabulary
pfName=handles.featureVocabulary.pfNameList{handles.pfNdx};
wfType = FeatureVocabularyForSelectFeatures.wfTypes{handles.winNdx};
wfParamName='min_window_radius';
handles.featureVocabulary.setWFParam(pfName,wfType,wfParamName,curVal);

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


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


% -------------------------------------------------------------------------
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
fv=handles.featureVocabulary;
wfType = fv.wfTypes{handles.winNdx};
%handles.data{handles.pfNdx}.(wfType).values.max_window_radius = curVal;

% set in featureVocabulary
pfName=handles.featureVocabulary.pfNameList{handles.pfNdx};
wfType = FeatureVocabularyForSelectFeatures.wfTypes{handles.winNdx};
wfParamName='max_window_radius';
handles.featureVocabulary.setWFParam(pfName,wfType,wfParamName,curVal);

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


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


% -------------------------------------------------------------------------
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
fv=handles.featureVocabulary;
wfType = fv.wfTypes{handles.winNdx};
%handles.data{handles.pfNdx}.(wfType).values.nwindow_radii = curVal;

% set in featureVocabulary
pfName=handles.featureVocabulary.pfNameList{handles.pfNdx};
wfType = FeatureVocabularyForSelectFeatures.wfTypes{handles.winNdx};
wfParamName='nwindow_radii';
handles.featureVocabulary.setWFParam(pfName,wfType,wfParamName,curVal);

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


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


% -------------------------------------------------------------------------
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
fv=handles.featureVocabulary;
wfType = fv.wfTypes{handles.winNdx};
%handles.data{handles.pfNdx}.(wfType).values.window_offsets = curVal;

% set in featureVocabulary
pfName=handles.featureVocabulary.pfNameList{handles.pfNdx};
wfType = FeatureVocabularyForSelectFeatures.wfTypes{handles.winNdx};
wfParamName='window_offsets';
handles.featureVocabulary.setWFParam(pfName,wfType,wfParamName,curVal);

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


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


% -------------------------------------------------------------------------
% --- Executes on button press in TransNone.
function TransNone_Callback(hObject, eventdata, handles)
% hObject    handle to TransNone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransNone
fv=handles.featureVocabulary;
curFn = fv.wfTypes{handles.winNdx};
handles = guidata(hObject);
% curT = handles.data{handles.pfNdx}.(curFn).values.trans_types;
% if get(hObject,'Value')
%   %curT=bitor(1,curT);
%   if ~any(strcmp('none',curT))
%    handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'none';
%   end
% else
% %   curT=bitand(14,curT);
%   allNdx = strcmp('none',curT);
%   handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
%   if isempty(handles.data{handles.pfNdx}.(curFn).values.trans_types),
% %  if handles.data{handles.pfNdx}.(curFn).values.trans_types==0,
%     warndlg('Select at least one transformation type');
%   end
% end

% update featureVocabulary
pfName=handles.featureVocabulary.pfNameList{handles.pfNdx};
wfType = fv.wfTypes{handles.winNdx};
transformation='none';
if get(hObject,'Value')
  fv.addWFTransformation(pfName,wfType,transformation);
else
  fv.removeWFTransformation(pfName,wfType,transformation);
end
% Update the translation checkbox to reflect the model (the change can fail
% to happen if the user is trying to turn off the only transformation)
checked=fv.isWFTransformation(pfName,wfType,transformation);
set(hObject,'Value',checked);

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in TransFlip.
function TransFlip_Callback(hObject, eventdata, handles)
% hObject    handle to TransFlip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransFlip
fv=handles.featureVocabulary;
curFn = fv.wfTypes{handles.winNdx};
handles = guidata(hObject);
% curT = handles.data{handles.pfNdx}.(curFn).values.trans_types;
% if get(hObject,'Value')
% %   curT=bitor(1,curT);
%   if ~any(strcmp('flip',curT))
%    handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'flip';
%   end
% else
% %   curT=bitand(14,curT);
%   allNdx = find(strcmp('flip',curT));
%   handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
%   if isempty(handles.data{handles.pfNdx}.(curFn).values.trans_types),
% %   if handles.data{handles.pfNdx}.(curFn).values.trans_types==0,
%     warndlg('Select at least one transformation type');
%   end
% end

% update featureVocabulary
featureVocabulary=handles.featureVocabulary;  % a ref
pfName=handles.featureVocabulary.pfNameList{handles.pfNdx};
wfType = fv.wfTypes{handles.winNdx};
transformation='flip';
if get(hObject,'Value')
  featureVocabulary.addWFTransformation(pfName,wfType,transformation);
else
  featureVocabulary.removeWFTransformation(pfName,wfType,transformation);
end
% Update the translation checkbox to reflect the model (the change can fail
% to happen if the user is trying to turn off the only transformation)
checked=featureVocabulary.isWFTransformation(pfName,wfType,transformation);
set(hObject,'Value',checked);

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in TransAbs.
function TransAbs_Callback(hObject, eventdata, handles)
% hObject    handle to TransAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransAbs

fv=handles.featureVocabulary;
%curFn = fv.wfTypes{handles.winNdx};
handles = guidata(hObject);
% curT = handles.data{handles.pfNdx}.(curFn).values.trans_types;
% if get(hObject,'Value')
% %   curT=bitor(1,curT);
%   if ~any(strcmp('abs',curT))
%    handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'abs';
%   end
% else
% %   curT=bitand(14,curT);
%   allNdx = find(strcmp('abs',curT));
%   handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
%   if isempty(handles.data{handles.pfNdx}.(curFn).values.trans_types),
% %   if handles.data{handles.pfNdx}.(curFn).values.trans_types==0,
%     warndlg('Select at least one transformation type');
%   end
% end

% update featureVocabulary
featureVocabulary=handles.featureVocabulary;  % a ref
pfName=handles.featureVocabulary.pfNameList{handles.pfNdx};
wfType = fv.wfTypes{handles.winNdx};
transformation='abs';
if get(hObject,'Value')
  featureVocabulary.addWFTransformation(pfName,wfType,transformation);
else
  featureVocabulary.removeWFTransformation(pfName,wfType,transformation);
end
% Update the translation checkbox to reflect the model (the change can fail
% to happen if the user is trying to turn off the only transformation)
checked=featureVocabulary.isWFTransformation(pfName,wfType,transformation);
set(hObject,'Value',checked);

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in TransRel.
function TransRel_Callback(hObject, eventdata, handles)
% hObject    handle to TransRel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransRel

fv=handles.featureVocabulary;
curFn = fv.wfTypes{handles.winNdx};
handles = guidata(hObject);
% curT = handles.data{handles.pfNdx}.(curFn).values.trans_types;
% if get(hObject,'Value')
% %   curT=bitor(1,curT);
%   if ~any(strcmp('relative',curT))
%    handles.data{handles.pfNdx}.(curFn).values.trans_types{end+1} = 'relative';
%   end
% else
% %   curT=bitand(14,curT);
%   allNdx = find(strcmp('relative',curT));
%   handles.data{handles.pfNdx}.(curFn).values.trans_types(allNdx) = [];
%   if isempty(handles.data{handles.pfNdx}.(curFn).values.trans_types),
% %   if handles.data{handles.pfNdx}.(curFn).values.trans_types==0,
%     warndlg('Select at least one transformation type');
%   end
% end

% update featureVocabulary
featureVocabulary=handles.featureVocabulary;  % a ref
pfName=handles.featureVocabulary.pfNameList{handles.pfNdx};
wfType = fv.wfTypes{handles.winNdx};
transformation='relative';
if get(hObject,'Value')
  featureVocabulary.addWFTransformation(pfName,wfType,transformation);
else
  featureVocabulary.removeWFTransformation(pfName,wfType,transformation);
end
% Update the translation checkbox to reflect the model (the change can fail
% to happen if the user is trying to turn off the only transformation)
checked=featureVocabulary.isWFTransformation(pfName,wfType,transformation);
set(hObject,'Value',checked);

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, eventdata, handles)
% hObject    handle to push_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%uiresume(handles.figure1);
delete(handles.figure1);
return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
basicData = get(handles.basicTable,'Data');
featureWindowSize = str2double(get(handles.editSize,'String'));
%windowFeaturesParams = convertData(handles);
windowFeaturesParams = handles.featureVocabulary.getInJLabelDataFormat();
set(handles.output,'Visible','off');
% handles.jld.UpdatePerframeParams(windowFeaturesParams, ...
%                                  basicData, ...
%                                  featureWindowSize);
% notify the JLabel 'object' that SelectFeatures is done
% Currently, we pass basicData and featureWindowSize back to JLabel so 
% that they can be restored to their current state if user does Select
% Features... again.  This is a hack, and it would be nice to fix at some
% point.  That table and the featureWindowSize should be calculated from
% the feature vocabulary.
JLabel('selectFeaturesDone', ...
       handles.figureJLabel, ...
       windowFeaturesParams, ...
       basicData, ...
       featureWindowSize);
return


% -------------------------------------------------------------------------
function ExtraParams_Callback(hObject, eventdata, handles)
% hObject    handle to ExtraParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ExtraParams as text
%        str2double(get(hObject,'String')) returns contents of ExtraParams as a double

fv=handles.featureVocabulary;
if isempty(fv.wfExtraParamNames{handles.winNdx}), return; end

str = get(hObject,'String');
str = strtrim(str);
newValue = str2num(str);
extraParam = fv.wfExtraParamNames{handles.winNdx};
wfType = fv.wfTypes{handles.winNdx};
%handles.data{handles.pfNdx}.(wfType).values.(extraParam) = newValue;

% set in featureVocabulary
pfName=handles.featureVocabulary.pfNameList{handles.pfNdx};
wfType = handles.featureVocabulary.wfTypes{handles.winNdx};
wfParamName=handles.featureVocabulary.wfExtraParamNames{handles.winNdx};
handles.featureVocabulary.setWFParam(pfName,wfType,wfParamName,newValue);

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


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
return


% -------------------------------------------------------------------------
% --- Executes on button press in Save.
function Save_Callback(hObject, eventdata, handles)
% hObject    handle to Save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Use project name.
paramsdir = deployedRelative2Global('params');
[fName,pName] = uiputfile(fullfile(paramsdir,'*.xml'),'Save feature configurations to..');
if ~fName
  return;
end

params = convertData(handles);
basicData = get(handles.basicTable,'Data');
featureWindowSize = round(str2double(get(handles.editSize,'String')));
docNode = createParamsXML(params,basicData,featureWindowSize);
fName = fullfile(pName,fName);
xmlwrite(fName,docNode);
return


% -------------------------------------------------------------------------
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
  ReadPerFrameParams(featureparamsfilename,handles.jld.featureConfigFile);
guidata(hObject,handles);
initData(hObject,params);
set(handles.basicTable,'Data',basicTable);
set(handles.editSize,'String',num2str(windowSize));
return


% -------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_copy_windowparams.
function popupmenu_copy_windowparams_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_copy_windowparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_copy_windowparams contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_copy_windowparams


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_copy_windowparams.
function pushbutton_copy_windowparams_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_copy_windowparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.pfNdx) || isempty(handles.winNdx),
  return;
end
fv=handles.featureVocabulary;
windowComp = fv.wfTypes;

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

%handles = CopyWindowParams(handles,pfNdx,winNdxFrom,pfNdx,winNdxTo);
handles.featureVocabulary.copyWFParams(pfNdx,winNdxFrom,pfNdx,winNdxTo);

% Update the view
updateWinParams(handles);
updateWinParamsEnablement(handles);

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_copy_windowtypes.
function popupmenu_copy_windowtypes_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_copy_windowtypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_copy_windowtypes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_copy_windowtypes


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_copy_windowtypes.
function pushbutton_copy_windowtypes_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to pushbutton_copy_windowtypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.pfNdx),
  return;
end

fv=handles.featureVocabulary;
pfNameList = fv.pfNameList;
windowComp = fv.wfTypes;

% which perframe fn are we copying to?
pfNdxTo = handles.pfNdx;
if ~fv.pfIsInVocabulary(pfNdxTo),
  return;
end

% which perframe fn are we copying from?
pfNdxFrom = get(handles.popupmenu_copy_windowtypes,'Value');

% --- selected?
if pfNdxFrom > numel(pfNameList),
  return;
end

% if ~handles.data{pfNdxFrom}.valid,
%   return;
% end

for winfnNdx = 1:numel(windowComp),
  %handles = CopyWindowParams(handles,pfNdxFrom,winfnNdx,pfNdxTo,winfnNdx);
  handles.featureVocabulary.copyWFParams(pfNdxFrom,winfnNdx,pfNdxTo,winfnNdx);
end

% update the view
updateWindowTable(handles);
updateWinParams(handles);
updateWinParamsEnablement(handles);
% curFn = fv.wfTypes{handles.winNdx};
% if handles.data{handles.pfNdx}.(curFn).valid,
%   enableWinParams(handles);
% end

guidata(hObject,handles);
updatePFCategoryAmountForCurrentPF(handles);
return


% % -------------------------------------------------------------------------
% function handles = CopyWindowParams(handles,pfNdxFrom,winfnNdxFrom,pfNdxTo,winfnNdxTo)
% 
% fv=handles.featureVocabulary;
% curFnFrom = fv.wfTypes{winfnNdxFrom};
% % something to copy from?
% if ~isfield(handles.data{pfNdxFrom},curFnFrom),
%   return;
% end
% curFnTo = fv.wfTypes{winfnNdxTo};
% handles.data{pfNdxTo}.(curFnTo).valid = handles.data{pfNdxFrom}.(curFnFrom).valid;
% for winParamsNdx = 1:numel(handles.winParams),
%   curType = handles.winParams{winParamsNdx};
%   handles.data{pfNdxTo}.(curFnTo).values.(curType) = ...
%     handles.data{pfNdxFrom}.(curFnFrom).values.(curType);
% end
% if ~isempty(fv.wfExtraParamNames{winfnNdxTo})
%   extraParam = fv.wfExtraParamNames{winfnNdxTo};
%   if isfield(handles.data{pfNdxFrom}.(curFnFrom).values,extraParam)
%     handles.data{pfNdxTo}.(curFnTo).values.(extraParam) = ...
%       handles.data{pfNdxFrom}.(curFnFrom).values.(extraParam);
%   else
%     handles.data{pfNdxTo}.(curFnTo).values.(extraParam) = '';
%   end
% end
% % updateWinParams(handles,winfnNdxTo);
% % if handles.data{pfNdxTo}.(curFnTo).valid,
% %   enableWinParams(handles);
% % end
% return


% % -------------------------------------------------------------------------
% function handles = CopyDefaultWindowParams(handles,wfAmount,pfNdxTo,winfnNdx)
% % For the per-frame feature with index pfNDxTo, sets the window features 
% % with type given by index winfnNdx to the amount given by wfAmount (one of
% % 'normal', 'more', or 'less').
% 
% fv=handles.featureVocabulary;
% curFn = fv.wfTypes{winfnNdx};
% % something to copy from?
% if ~isfield(handles.wfParamsFromAmount.(wfAmount),curFn) &&...
%     isfield(handles.data{pfNdxTo},curFn),
%     handles.data{pfNdxTo}.(curFn).valid = false;
%   return;
% end
% handles.data{pfNdxTo}.(curFn).valid = handles.wfParamsFromAmount.(wfAmount).(curFn).valid;
% for winParamsNdx = 1:numel(handles.winParams),
%   curType = handles.winParams{winParamsNdx};
%   handles.data{pfNdxTo}.(curFn).values.(curType) = ...
%     handles.wfParamsFromAmount.(wfAmount).(curFn).values.(curType);
% end
% if ~isempty(fv.wfExtraParamNames{winfnNdx})
%   extraParam = fv.wfExtraParamNames{winfnNdx};
%   handles.data{pfNdxTo}.(curFn).values.(extraParam) = ...
%     handles.wfParamsFromAmount.(wfAmount).(curFn).values.(extraParam);
% end
% return


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
function handles = UpdateDescriptionPanels(handles)

if isempty(handles.pfNdx),
  return;
end

fv=handles.featureVocabulary;

% update visibility
if strcmpi(handles.currentTab,'description')
  set(handles.uipanel_description,'Visible','on');
  set(handles.uipanel_histogram,'Visible','off');
else
  set(handles.uipanel_description,'Visible','off');
  set(handles.uipanel_histogram,'Visible','on');
end

% histogram if necessary
if strcmpi(handles.currentTab,'perframehistogram') && ...
   (handles.histogramData.lastPfNdx ~= handles.pfNdx || ...
    ~strcmpi(handles.histogramData.lastType,'perframe') ) && ...
   (handles.jld.nexps>0),    
  i = find(handles.histogramData.perframe_idx == handles.pfNdx,1);
  if isempty(i),
    i = numel(handles.histogramData.perframe_idx)+1;
    [handles.histogramData.hhist,...
      ~,~,hleg,hxlabel,hylabel,...
      handles.histogramData.frac{i},handles.histogramData.frac_outside{i},...
      handles.histogramData.edges{i},handles.histogramData.centers_plot{i}] = ...
      HistogramPerFrameFeature(handles.jld,fv.pfNameList{handles.pfNdx},...
      'axes',handles.axes_histogram,...
      'unknowncolor','w',...
      'labelcolors',jet(handles.jld.nbehaviors)*.7);
    handles.histogramData.perframe_idx(i) = handles.pfNdx;
  else
    [handles.histogramData.hhist,...
      ~,~,hleg,hxlabel,hylabel,...
      handles.histogramData.frac{i},handles.histogramData.frac_outside{i},...
      handles.histogramData.edges{i},handles.histogramData.centers_plot{i}] = ...
      HistogramPerFrameFeature(handles.jld,fv.pfNameList{handles.pfNdx},...
      'axes',handles.axes_histogram,...
      'edges',handles.histogramData.edges{i},...
      'frac',handles.histogramData.frac{i},...
      'frac_outside',handles.histogramData.frac_outside{i},...
      'unknowncolor','w',...
      'labelcolors',jet(handles.jld.nbehaviors)*.7);
  end
  handles.histogramData.lastPfNdx = handles.pfNdx;
  handles.histogramData.lastType = 'perframe';
  textcolor = get(handles.togglebutton_tabdescription,'ForegroundColor');
  set(handles.axes_histogram,'XColor',textcolor,'YColor',textcolor,'Color','k','Box','off');
  set(hxlabel,'Color',textcolor);
  set(hylabel,'Color',textcolor);
  set(hleg,'Color','k','TextColor',textcolor,'Box','off','Location','Best');
end


% -------------------------------------------------------------------------
function editSize_Callback(hObject, eventdata, handles)
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSize as text
%        str2double(get(hObject,'String')) returns contents of editSize as a double
curVal = str2double(get(hObject,'String'));
if isempty(curVal)||(round(curVal)-curVal)~=0
  msgbox('Enter numerical values. eg: "10" (without with quotes)');
  return;
end

% First set the default window feature level values

fv=handles.featureVocabulary;
wfAmounts = fieldnames(handles.wfParamsFromAmount);
wfTypes = fv.wfTypes;

% Iterate over amounts and window-feature types, and set the max window
% radius for all of the pre-set amounts to whatever the user entered.
for i = 1:numel(wfAmounts)
  wfAmount = wfAmounts{i};
  for j = 1:numel(wfTypes)
    wfType=wfTypes{j};
    wfParams=handles.wfParamsFromAmount.(wfAmount)
    if isfield(wfParams,wfType),
      handles.wfParamsFromAmount.(wfAmount).(wfType).values.max_window_radius = curVal;
    end
  end
end

% Do the same for the featureVocabulary
handles.featureVocabulary.setMaxWindowRadiusForAllWFAmounts(curVal);

% Now copy the default values to the perframe features.
basicData = get(handles.basicTable,'Data');
for ndx = 1:size(basicData,1)
  if strcmp(basicData{ndx,2},'all')
    handles = applyCategoryType(handles,ndx);
    handles.featureVocabulary.setPFCategoryToWFAmount(ndx,basicData{ndx,3});
  end
end
guidata(hObject,handles);
createPfTable(handles.pfTable);
return


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
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


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_hist.
function pushbutton_hist_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to pushbutton_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prcEdges = [5 15 30 50 70 85 95];

fv=handles.featureVocabulary;
histfnNdx = find(strcmp('hist',fv.wfTypes));
histExtraName = fv.wfExtraParamNames{histfnNdx};  %#ok

h = waitbar(0,'Computing hist bins');
for ndx = 1:numel(fv.pfNameList)
  curPf = fv.pfNameList{ndx};
  waitbar(ndx/numel(fv.pfNameList),h);
  
  if ~handles.data{ndx}.valid || ~handles.data{ndx}.hist.valid; 
    continue;
  end
  
  allData = [];
  for expi = 1:handles.jld.nexps,
    
    % load per-frame data for this experiment
    perframedir = handles.jld.GetFile('perframedir',expi);
    file = fullfile(perframedir,[curPf,'.mat']);
    if ~exist(file,'file'),
      warning('Per-frame data file %s does not exist',file);
      continue;
    end
    
    perframedata = load(file);
    
    for fly = 1:handles.jld.nflies_per_exp(expi),
      
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
delete(h);
guidata(hObject,handles);
if ~isempty(handles.pfNdx)
  updateWindowTable(handles);
end
return


% -------------------------------------------------------------------------
function CloseRequestFcn(hObject,eventdata,handles)
push_cancel_callback(hObject,eventdata,handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_ok.
function pushbutton_ok_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ok (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% configfile = handles.jld.configfilename;
% [~,~,ext ] = fileparts(configfile);
% if strcmpi(ext,'.xml'),
%   configparams = ReadXMLParams(configfile);
% elseif strcmpi(ext,'.mat'),
%   configparams = load(configfile);
% else
%   errordlg('Project file %s is invalid. Cannot save the window features',configfile);
% end
% if ~isfield(configparams.file,'featureparamfilename') || isempty(configparams.file.featureparamfilename)
%   behaviorname = configparams.behaviors.names;
%   if iscell(behaviorname),
%     defaultname = sprintf('WindowFeatures_%s.xml',behaviorname{1});
%   else
%     defaultname = sprintf('WindowFeatures_%s.xml',behaviorname);
%   end
%   [fname,fpath]= uiputfile(fullfile('params','*.xml'),'Enter a name for feature config file',defaultname);
%   if isempty(fname),return, end
%   featureconfigfile = fullfile(fpath,fname);
%   configparams.file.featureparamfilename = featureconfigfile;
%   docNode = com.mathworks.xml.XMLUtils.createDocument('params');
%   toc = docNode.getDocumentElement;
%   fnames = fieldnames(configparams);
%   for ndx = 1:numel(fnames)
%     toc.appendChild(createXMLNode(docNode,fnames{ndx},configparams.(fnames{ndx})));
%   end
%   xmlwrite(configfile,docNode);
% end
% 
% featureconfigfile = configparams.file.featureparamfilename;
% params = convertData(handles);
% basicData = get(handles.basicTable,'Data');
% featureWindowSize = round(str2double(get(handles.editSize,'String')));
% docNode = createParamsXML(params,basicData,featureWindowSize);
% xmlwrite(featureconfigfile,docNode);

pushbutton_done_Callback(hObject,eventdata,handles);
push_cancel_Callback(hObject,eventdata,handles);
return
