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

% Last Modified by GUIDE v2.5 01-Apr-2013 01:56:42

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


% -------------------------------------------------------------------------
% --- Executes just before SelectFeatures is made visible.-----------------
function SelectFeatures_OpeningFcn(hObject, eventdata, handles, varargin)  %#ok
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SelectFeatures (see VARARGIN)

% So that we can't be closed before we're done opening
handles.doneWithOpeningFunction=false;
guidata(hObject,handles);

% Process input arguments
figureJLabel = varargin{1};

%
% Populate various aspects that are model-like
%

% Need to keep the parent JLabel's JLabelData object around (a ref to it,
% really), for various annoying reasons
jld = JLabel('getJLabelData',figureJLabel);
handles.jld=jld;

% Get the cached version of maxWindowRadiusCommon out of JLabel
maxWindowRadiusCommonCached= ...
  JLabel('getMaxWindowRadiusCommonCached', ...
         figureJLabel);

% This is the main aspect of the model---the feature vocabulary, in a form
% well-suited to SelectFeatures
featureLexicon = jld.featureLexicon;
scoreFeatures = jld.scoreFeatures;
toBeCalculatedPFNames=jld.allperframefns;
windowFeatureParams = jld.GetPerframeParams();
handles.featureVocabulary= ...
  FeatureVocabularyForSelectFeatures(featureLexicon, ...
                                     scoreFeatures, ...
                                     toBeCalculatedPFNames, ...
                                     windowFeatureParams, ...
                                     maxWindowRadiusCommonCached);

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

% Want to keep this around so we can send a "message" to it when user
% clicks on "Done"
handles.figureJLabel=figureJLabel;



%
% Initialize things that would be part of the View in an MVC arrangement
%

% Set the initial selected per-frame feature and window-feature type
% (indexes to them, actually)
handles.pfNdx = [];
handles.wfTypeNdx = [];

% Compute the figure size in basic and advanced mode, store them
curPos = get(handles.figure_SelectFeatures,'Position');
tablePos = get(handles.basicTable,'Position');
reducedWidth = tablePos(1)+tablePos(3) + 15;
reducedHeight = curPos(4);
handles.advancedSize = curPos(3:4);
handles.basicSize = [reducedWidth reducedHeight];

% Set the initial mode: basic or advanced.  This also sets the figure size
% appropriately.
handles=setViewMode(handles,'basic');

% Center SelectFeatures on the JLabel window
centerOnParentFigure(handles.figure_SelectFeatures,figureJLabel);

% This is the end of initializing the view variables that are not actualy
% graphics objects (although setViewMode() does manipulate some graphics
% objects)

% Commit changes to the figure guidata
guidata(hObject,handles);



%
% Update the graphics objects to match the model and the
% non-graphics-object view instance variables
%

% Change a few things so they still work well on Mac
adjustColorsIfMac(hObject);
% Create the basic table uitable
initializeBasicTable(hObject);
% Put an initial max window radius in the appropriate editbox in the
% "Basic" side of the UI
% And other stuff...
initializeWindowTable(hObject);
initializeCopyFromMenus(hObject);
initializeDescriptionPanels(hObject);
compatibleBasicAdvanced(handles);
%set(hObject,'Visible','on');  % have to do this b/c of Java hacking
initializePfTable(hObject);
%removeRowHeaders(hObject);

% Make sure the graphics items are consistent with the model state
updateMaxWindowRadiusEditBox(handles);
updateWindowTableAndEnablement(handles);
updateWinParamsAndEnablement(handles);

% OK, now it's OK to close the figure
handles=guidata(hObject);
handles.doneWithOpeningFunction=true;
guidata(hObject,handles);

return


% -------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = SelectFeatures_OutputFcn(hObject, eventdata, handles)  %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.figure_SelectFeatures;
return


% % -------------------------------------------------------------------------
% function setJLDobj(hObject,jld)
% % Set the JLabelDataObject.
% handles = guidata(hObject);
% % handles.JLDobj = JLDobj;
% 
% % handles.windowComp = {'default','mean','min','max','hist','prctile',...
% %    'change','std','harmonic','diff_neighbor_mean',...
% %    'diff_neighbor_min','diff_neighbor_max','zscore_neighbors'};
% % 
% % handles.winextraParams = {'','','','','hist_edges','prctiles','change_window_radii',...
% %   '','num_harmonic','','','',''};
% % 
% % handles.winParams = {'max_window_radius','min_window_radius','nwindow_radii',...
% %   'trans_types','window_offsets'};
% % 
% % handles.winextraDefaultParams = {[],[],[],[],[-400000 0 40000],[5 10 30 50 70 90 95],[1 3],...
% %   [],[2],[],[],[],[]};
% % handles.defaultWinParams = {10,1,3,{'none'},0};
% 
% guidata(hObject,handles);
% 
% windowFeatureParams = jld.GetPerframeParams();
% featureLexicon= jld.featureLexicon;
% initData(hObject,windowFeatureParams,featureLexicon,jld);
% 
% return


% % -------------------------------------------------------------------------
% % Initialize the data structure.
% function initData(hObject,featureLexicon,jld)
% readFeatureLexicon(hObject,featureLexicon);
% createBasicTable(hObject,jld);
% 
% handles = guidata(hObject);
% if ~isempty(jld.featureWindowSize)
%   set(handles.editSize,'String',num2str(jld.featureWindowSize));
% end
% 
% % initialize histogramData
% handles.histogramData = struct;
% handles.histogramData.lastPfNdx = nan;
% handles.histogramData.lastType = '';
% handles.histogramData.perframe_idx = [];
% handles.histogramData.hhist = [];
% handles.histogramData.frac = {};
% handles.histogramData.frac_outside = {};
% handles.histogramData.edges = {};
% handles.histogramData.centers_plot = {};
% 
% % handles.data = data;
% handles.pfNdx = [];
% handles.wfTypeNdx = [];
% guidata(hObject,handles);
% createPfTable(hObject);
% createWindowTable(hObject);
% createCopyFromMenus(hObject);
% createDescriptionPanels(hObject);
% compatibleBasicAdvanced(handles);
% return


% -------------------------------------------------------------------------
function initializeWindowTable(hObject)
% Sets values for the window table.

handles = guidata(hObject);

% Deal with windowTable
fv=handles.featureVocabulary;
wfTypes=fv.wfTypes;
set(handles.windowTable, ...
    'ColumnName',{'Computation Type','Select'});

% Initialize the table data  
nWFTypes=length(wfTypes);
nRows=nWFTypes-1;  % don't want to show 'default' pseudo WF type
windowData = cell(nRows,2);
for rowIndex = 1:nRows
  wfTypeIndex=rowIndex+1;
  wfType = fv.wfTypes{wfTypeIndex};
  windowData{rowIndex,1} = wfType;
  windowData{rowIndex,2} = false;
end
set(handles.windowTable,'Data',windowData);

set(handles.windowTable','ColumnWidth',{135 'auto'});
set(handles.windowTable,'ColumnEditable',[false,true]);
set(handles.windowTable,'CellSelectionCallback',@windowSelect);
set(handles.windowTable,'CellEditCallback',@windowEdit);

return


% -------------------------------------------------------------------------
function initializePfTable(hObject)
% Initialize the per-frame feature table.
handles = guidata(hObject);
set(handles.pfTable,'ColumnName',{'Features','Select','Amount','Category'});
set(handles.pfTable,'ColumnEditable',[false,true,true,false]);
set(handles.pfTable, ...
    'ColumnFormat', ...
    {'char',...
     'logical', ...
     [handles.featureVocabulary.wfAmounts' 'custom'], ...
     'char'} );
set(handles.pfTable,'ColumnWidth',{190,50,75,95});
set(handles.pfTable,'CellSelectionCallback',@pfSelect);
set(handles.pfTable,'CellEditCallback',@pfEdit);

% init the static cols of the table
fv=handles.featureVocabulary;  % a ref
pfNameList = fv.subdialectPFNames;
nPFs=length(pfNameList);
tableData = cell(nPFs,4);
for ndx = 1:nPFs
  tableData{ndx,1} = pfNameList{ndx};
  tableData{ndx,2} = false;  % place-holder
  tableData{ndx,3}='normal';  % place-holder
  pfCategories=fv.pfCategoriesFromName.(pfNameList{ndx});
  str = sprintf('%s',pfCategories{1});
  for sndx = 2:numel(pfCategories)
    str = sprintf('%s,%s',str,pfCategories{sndx});
  end
  tableData{ndx,4} = str;
end
set(handles.pfTable,'Data',tableData);

% Now do an update to properly populate the dynamic cols
updatePFTable(handles);
return


% -------------------------------------------------------------------------
function updatePFTable(handles)
% Updates the contents of the per-frame feature table.

% Generate the table data, based on the feature vocabulary
fv=handles.featureVocabulary;  % a ref
pfList = fv.subdialectPFNames;
nPFs=length(pfList);
tableData=get(handles.pfTable,'Data');
for ndx = 1:nPFs
  tableData{ndx,2} = fv.vocabulary{ndx}.enabled;
  tableData{ndx,3}=fv.getWFAmountForPF(ndx);
end

% Update the graphics object
set(handles.pfTable,'Data',tableData);

% % do Java stuff to select the corrent element
% jscrollpane = findjobj(handles.pfTable);
% jtable = jscrollpane.getViewport.getView;
% pfNdx=handles.pfNdx;
% if isempty(pfNdx)
%   jtable.changeSelection(0,0, false, false);
% else
%   jtable.changeSelection(pfNdx-1,0, false, false);
% end  

return


% -------------------------------------------------------------------------
function updatePFTableForCurrentPF(handles)
% Updates the contents of the per-frame feature table, for just the current
% PF.

% Generate the table data, based on the feature vocabulary
fv=handles.featureVocabulary;  % a ref
pfNdx=handles.pfNdx;
% pfList = fv.subdialectPFNames;
% nPFs=length(pfList);
tableData=get(handles.pfTable,'Data');
tableData{pfNdx,2} = fv.vocabulary{pfNdx}.enabled;
tableData{pfNdx,3}=fv.getWFAmountForPF(pfNdx);
set(handles.pfTable,'Data',tableData);

% % do Java stuff to select the current element
% jscrollpane = findjobj(handles.pfTable);
% jtable = jscrollpane.getViewport.getView;
% pfNdx=handles.pfNdx;
% if isempty(pfNdx)
%   jtable.changeSelection(0,0, false, false);
% else
%   jtable.changeSelection(pfNdx-1,0, false, false);
% end  

return


% -------------------------------------------------------------------------
function initializeBasicTable(hObject)
handles = guidata(hObject);

% Set the simple parts of the uitable
set(handles.basicTable,'ColumnName',{'Categories','Select','Amount'});
set(handles.basicTable, ...
    'ColumnFormat', ...
    {'char',...
     {'all' 'none' 'custom'}, ...
     [handles.featureVocabulary.wfAmounts' 'custom']} );
set(handles.basicTable,'ColumnEditable',[false,true,true]);
set(handles.basicTable,'ColumnWidth',{85,65,75});
set(handles.basicTable,'CellSelectionCallback',@basicSelect);
set(handles.basicTable,'CellEditCallback',@basicEdit);

% Get a list of the categories that actually have PFs in them (e.g. scores
% category often has no PFs in it)
fv=handles.featureVocabulary;
pfCategoryNames=fv.pfCategoryNames;
categoryIsNonEmpty= ...
  @(pfCategoryName)(~isempty(fv.getPFNamesInCategory(pfCategoryName)));
pfCategoryNamesNonEmpty= ...
  cellFilter(categoryIsNonEmpty,pfCategoryNames);

% Populate the static parts of the table cells
nPFCategoriesNonEmpty=length(pfCategoryNamesNonEmpty);
tableData = cell(nPFCategoriesNonEmpty,3);
for i = 1:nPFCategoriesNonEmpty
  tableData{i,1} = pfCategoryNamesNonEmpty{i};
  tableData{i,2} = 'none';  % placeholder
  tableData{i,3} = 'normal';  % placeholder
end
set(handles.basicTable,'Data',tableData);

% Call update to populate the dynamic parts of the table
updateBasicTable(handles);
return


% -------------------------------------------------------------------------
function updateBasicTable(handles)
fv=handles.featureVocabulary;
tableData = get(handles.basicTable,'Data');
pfCategoryNames=tableData(:,1);  % these are the nonempty categories
nPFCategories=length(pfCategoryNames);
for i=1:nPFCategories
  pfCategoryName=pfCategoryNames{i};
  tableData{i,2} = fv.getPFCategoryLevel(pfCategoryName);
  tableData{i,3} = fv.getWFAmountForPFCategory(pfCategoryName); 
end
set(handles.basicTable,'Data',tableData);
return


% % -------------------------------------------------------------------------
% function updateBasicTableOnePFCategory(handles,pfCategoryName)
% fv=handles.featureVocabulary;
% pfCategoryList=fv.pfCategoryNames;
% pfCategoryIndex=find(strcmp(pfCategoryName,pfCategoryList));
% tableData = get(handles.basicTable,'Data');
% tableData{pfCategoryIndex,2} = fv.getPFCategoryLevel(pfCategoryName);
% tableData{pfCategoryIndex,3} = fv.getWFAmountForPFCategory(pfCategoryName); 
% set(handles.basicTable,'Data',tableData);
% return


% -------------------------------------------------------------------------
function updateBasicTableAllCategoriesOfCurrentPF(handles)
fv=handles.featureVocabulary;
pfNameList=fv.subdialectPFNames;
pfName=pfNameList{handles.pfNdx};
pfCategoriesThisName=fv.pfCategoriesFromName.(pfName);
allPFCategories=fv.pfCategoryNames;
pfCategoryIndicesThisName = find(ismember(allPFCategories,pfCategoriesThisName));
basicData = get(handles.basicTable,'Data');
pfCategoryNamesInTable=basicData(:,1);
for i = 1:numel(pfCategoryIndicesThisName)
  pfCategoryName=pfCategoriesThisName{i};
  rowIndex=find(strcmp(pfCategoryName,pfCategoryNamesInTable));
  basicData{rowIndex,2} = fv.getPFCategoryLevel(pfCategoryName);
  basicData{rowIndex,3} = fv.getWFAmountForPFCategory(pfCategoryName); 
end
set(handles.basicTable,'Data',basicData);
return


% % -------------------------------------------------------------------------
% function initializeMaxWindowRadiusEditBox(handles,maxWindowRadius)
% if ~isempty(maxWindowRadius)
%   % set the global maxWindowRadius to the given value
%   fv=handles.featureVocabulary;
%   fv.setMaxWindowRadiusForAllWFs(maxWindowRadius);
%   fv.setMaxWindowRadiusForAllWFAmounts(maxWindowRadius);
% end
% % Update the view
% updateMaxWindowRadiusEditBox(handles);
% return


% -------------------------------------------------------------------------
function updateMaxWindowRadiusEditBox(handles)
fv=handles.featureVocabulary;
[thereIsConsensus,maxWindowRadiusConsensus]= ...
  fv.getConsensusMaxWindowRadiusForAllWFs();
if thereIsConsensus ,
  set(handles.editSize,'String',num2str(maxWindowRadiusConsensus), ...
                       'FontAngle','normal');
else
  set(handles.editSize,'String','Custom', ...
                       'FontAngle','italic');  
end
return


% % -------------------------------------------------------------------------
% function removeRowHeaders(hObject)
% % Tweaking the table. Use underlying java objects to do that. Found
% % this at http://undocumentedmatlab.com/blog/uitable-sorting/
% 
% handles = guidata(hObject);
% 
% % Basic Table
% jscrollpane = findjobj(handles.basicTable);
% jtable = jscrollpane.getViewport.getView;
% jtable.setSortable(false);	
% jtable.setAutoResort(false);
% jtable.setMultiColumnSortable(true);
% 
% % Set the size for the row headers.
% rowHeaderViewport=jscrollpane.getComponent(4);
% rowHeader=rowHeaderViewport.getComponent(0);
% newWidth=0; 
% rowHeaderViewport.setPreferredSize(java.awt.Dimension(newWidth,0));
% height=rowHeader.getHeight;
% rowHeader.setPreferredSize(java.awt.Dimension(newWidth,height));
% rowHeader.setSize(newWidth,height); 
% 
% % Pf Table.
% jscrollpane = findjobj(handles.pfTable);
% jtable = jscrollpane.getViewport.getView;
% jtable.setSortable(false);	
% jtable.setAutoResort(false);
% jtable.setMultiColumnSortable(false);
% 
% % Set the size for the row headers.
% rowHeaderViewport=jscrollpane.getComponent(4);
% rowHeader=rowHeaderViewport.getComponent(0);
% newWidth=0; 
% rowHeaderViewport.setPreferredSize(java.awt.Dimension(newWidth,0));
% height=rowHeader.getHeight;
% rowHeader.setPreferredSize(java.awt.Dimension(newWidth,height));
% rowHeader.setSize(newWidth,height); 
% 
% % Window Table.
% jscroll=findjobj(handles.windowTable);
% rowHeaderViewport=jscroll.getComponent(4);
% rowHeader=rowHeaderViewport.getComponent(0);
% rowHeader.setSize(80,360);
% 
% %resize the row header
% newWidth=0; %100 pixels.
% rowHeaderViewport.setPreferredSize(java.awt.Dimension(newWidth,0));
% height=rowHeader.getHeight;
% rowHeader.setPreferredSize(java.awt.Dimension(newWidth,height));
% rowHeader.setSize(newWidth,height); 
% return


% -------------------------------------------------------------------------
function initializeCopyFromMenus(hObject)

handles = guidata(hObject);
% can copy from any window feature type
fv=handles.featureVocabulary;
wfTypes=fv.wfTypes;
pfNameList=fv.subdialectPFNames;
set(handles.popupmenu_copy_windowparams, ...
    'String',[wfTypes,{'---'}],...
    'Value',numel(wfTypes)+1);
% can copy from any per-frame feature
set(handles.popupmenu_copy_windowtypes, ...
    'String',[pfNameList;{'---'}],...
    'Value',numel(pfNameList)+1);
return


% -------------------------------------------------------------------------
function initializeDescriptionPanels(hObject)

handles = guidata(hObject);

% which tab to show
handles.currentTab = 'perframehistogram';
set(handles.togglebutton_tabperframehistogram,'Value',1);

% set visibility
uipanel_tabs_SelectionChangeFcn(handles.uipanel_tabs, struct('NewValue',handles.togglebutton_tabdescription), handles);

guidata(hObject,handles);


% -------------------------------------------------------------------------
function basicSelect(hObject,eventData)  %#ok
return


% -------------------------------------------------------------------------
function basicEdit(hObject,eventData)
% the user selects the category
if isempty(eventData.Indices); return; end
oldPointer=pointerToWatch(gcbf);
if eventData.Indices(2)==2
  % User made a selection in the 2nd column, to select 'all' or 'none' of
  % the PFs in that category.  They can also select 'custom', but that does
  % nothing.
  pfCategoryLevelChanged(hObject,eventData);
elseif eventData.Indices(2) == 3
  pfCategoryWFAmountChanged(hObject,eventData);
end
restorePointer(gcbf,oldPointer);
return


% -------------------------------------------------------------------------
function pfCategoryLevelChanged(hObject,eventData)
  % User selected the level (all, none, custom) for a PF category in the
  % basic table.  Note that custom will be ignored, and the level will be
  % changed back to what it was before.
  handles = guidata(hObject);
  fv=handles.featureVocabulary;  % a ref
  rowIndex=eventData.Indices(1);
  basicData = get(handles.basicTable,'Data');
  pfCategoryName=basicData{rowIndex,1};
  wfAmount=basicData{rowIndex,3};
  newLevel=eventData.NewData;
  switch newLevel
    case 'none'
      h = waitbar(0,'Unselecting the perframe features');
      pfNames=fv.getPFNamesInCategory(pfCategoryName);
      for i = 1:length(pfNames)
        fv.disablePerframeFeature(pfNames{i});
      end            
      close(h);
    case 'all'
      h = waitbar(0,'Selecting the perframe features');
      if ~isequal(wfAmount,'custom')
        fv.setAllPFsInCategoryToWFAmount(pfCategoryName,wfAmount);
      end
      fv.enableAllPFsInCategory(pfCategoryName);
      close(h);
    case 'custom'
      % don't change the model at all---the custom item is only there for
      % when the user adds/deletes individual PFs using the Advanced tab.
  end
  guidata(hObject,handles);
  updateBasicTable(handles);
  updatePFTable(handles);
  updateWindowTableAndEnablement(handles);
  updateWinParamsEnablement(handles);
return


% -------------------------------------------------------------------------
function pfCategoryWFAmountChanged(hObject,eventData)
  % User selected the window-feature
  % amount for a PF category, in the basic table.  Note that selecting
  % 'custom' does nothing.
  handles = guidata(hObject);
  fv=handles.featureVocabulary;  % a ref
  rowIndex=eventData.Indices(1);
  basicData=get(handles.basicTable,'Data');
  pfCategoryName=basicData{rowIndex,1};
  newWFAmount=eventData.NewData;
  if ~isequal(newWFAmount,'custom'), 
    fv.setAllPFsInCategoryToWFAmount(pfCategoryName,newWFAmount);
  end
  updateBasicTable(handles);
  updatePFTable(handles);
  updateWindowTableAndEnablement(handles);
  updateWinParamsAndEnablement(handles);
return


% % -------------------------------------------------------------------------
% function handles = applyCategoryType(handles,iSelectedCategory)
% % get some variables out of handles
% basicData = get(handles.basicTable,'Data');
% wfAmount = basicData{iSelectedCategory,3};  
%   % The amount of window features to use.  Can be 'normal', 'more', or
%   % 'less'
% fv=handles.featureVocabulary;
% selectedCategory = handles.pfCategoryNames{iSelectedCategory};
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
% What does this do?  --ALT, Mar 31, 2013
basicTable = get(handles.basicTable,'Data');
incompatible = '';
fv=handles.featureVocabulary;
pfNameList=fv.subdialectPFNames;
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
% Called when user selects cells in the per-frame feature table.

% fprintf('pfSelect() called.\n');
if isempty(eventData.Indices)
  return;
end

oldPointer=pointerToWatch(gcbf);

handles = guidata(hObject);

%pfData = get(handles.pfTable,'Data');

if size(eventData.Indices,1)>1 ,
  % If multiple elements have been selected, select none of the PFs,
  % because, you know, what the hell are we supposed to do with that?
  handles.pfNdx=[];
  %disableWindowTable(handles);
else
  % the usual case---a single cell is selected
  pfNdx = eventData.Indices(1,1);
  handles.pfNdx = pfNdx;
end
handles = UpdateDescriptionPanels(handles);
guidata(hObject,handles);
updateWindowTableAndEnablement(handles);
updateWinParamsAndEnablement(handles);
restorePointer(gcbf,oldPointer);

return


% -------------------------------------------------------------------------
function pfEdit(hObject,eventData)
% Called when a perframe feature is added or removed in the pfTable.

oldPointer=pointerToWatch(gcbf);

% Update the figure guidata
handles = guidata(hObject);

pfNdx = eventData.Indices(1,1);
handles.pfNdx = pfNdx;
guidata(hObject,handles);

% Call the appropriate sub-handler
tableColumnIndex=eventData.Indices(2);
if tableColumnIndex==2,
  % User edited the enable/disable column
  pfEnablementChanged(handles,eventData);
elseif tableColumnIndex==3,
  % User edited the amount column
  pfAmountChanged(handles,eventData);
end

restorePointer(gcbf,oldPointer);

return


% -------------------------------------------------------------------------
function pfEnablementChanged(handles,eventData)
  % This is a "private" "method" of the "controller".
  fv=handles.featureVocabulary;  % a ref
  pfIndex=eventData.Indices(1);
  enabled=eventData.NewData;  
  fv.setPFEnablement(pfIndex,enabled);
  % Now update the view
  updateBasicTable(handles);
  updateWindowTableAndEnablement(handles);
  updateWinParamsAndEnablement(handles);
return


% -------------------------------------------------------------------------
function pfAmountChanged(handles,eventData)
  % This is a "private" "method" of the "controller".
  fv=handles.featureVocabulary;  % a ref
  pfIndex=eventData.Indices(1);
  newWFAmount=eventData.NewData;  
  fv.setPFToWFAmount(pfIndex,newWFAmount);
  % Now update the view
  updatePFTableForCurrentPF(handles);
  updateBasicTable(handles);
  updateWindowTableAndEnablement(handles);
  updateWinParamsAndEnablement(handles);
return


% -------------------------------------------------------------------------
function updateWindowTableAndEnablement(handles)
updateWindowTable(handles);
updateWindowTableEnablement(handles);
return


% -------------------------------------------------------------------------
function updateWindowTable(handles)
% Sets the data in window table based on selections made in pfTable.
pfNdx=handles.pfNdx;
fv=handles.featureVocabulary;
nWFTypes=length(fv.wfTypes);
%windowData = cell(nWFTypes,2);
windowData=get(handles.windowTable,'Data');
if isempty(pfNdx) ,
  pfParams=struct([]);  % empty struct with no fields
else
  pfParams = fv.vocabulary{pfNdx};
end
nRows=nWFTypes-1; % b/c default not shown
for rowIndex = 1:nRows
  wfTypeIndex=rowIndex+1;
  wfType = fv.wfTypes{wfTypeIndex};
  if isfield(pfParams,wfType),
    %windowData{wfTypeIndex,1} = wfType;
    windowData{rowIndex,2} = pfParams.(wfType).enabled;
  else
    % warning('This error is occurring!');  
      % this isn't an error anymore -- it can happen by design if pfNdx is
      % empty
    %windowData{wfTypeIndex,1} = wfType;
    windowData{rowIndex,2} = false;
  end
end
set(handles.windowTable,'Data',windowData);
% jscrollpane = findjobj(handles.windowTable);
% jtable = jscrollpane.getViewport.getView;
% wfTypeNdx=handles.wfTypeNdx;
% if isempty(wfTypeNdx)
%   jtable.changeSelection(0,0, false, false);
% else
%   jtable.changeSelection(wfTypeNdx-1,0, false, false);
% end  
return


% -------------------------------------------------------------------------
function updateWindowTableRow(handles,rowIndex)
% Updates the checkbox for wfType wfTypeIndex, i.e. one row of the window
% table.
wfTypeIndex=rowIndex+1;  % b/c 'default' not shown
pfNdx=handles.pfNdx;
fv=handles.featureVocabulary;
pfParams = fv.vocabulary{pfNdx};
windowData = get(handles.windowTable,'Data');
wfType = fv.wfTypes{wfTypeIndex};
if ~isfield(pfParams,wfType),
  warning('This error is occurring!');
  windowData{rowIndex,2} = false;
else
  windowData{rowIndex,2} = pfParams.(wfType).enabled;
end
set(handles.windowTable,'Data',windowData);
return


% -------------------------------------------------------------------------
function windowSelect(hObject,eventData)
% Called when user selects cells in windowTable.
% Also seems to get called when the windowTable has been disabled and then
% gets enabled...
%fprintf('In windowSelect()\n');

% When something changes in the pfTable window.
if isempty(eventData.Indices)
  return;
end

oldPointer=pointerToWatch(gcbf);

handles = guidata(hObject);
%winData = get(handles.windowTable,'Data');

if (size(eventData.Indices,1)>1)
  % If more than one cell selected, just don't even try to deal.
  handles.wfTypeNdx = [];
  %disableWinParams(handles);
else
  rowIndex = eventData.Indices(1,1);
  handles.wfTypeNdx = rowIndex+1;  % +1 b/c default is not shown
%   if winData{ndx,2},
%     updateWinParams(handles);
%     enableWinParams(handles);
%   else
%     disableWinParams(handles);
%   end
end
guidata(hObject,handles);
updateWinParamsAndEnablement(handles);
restorePointer(gcbf,oldPointer);
return


% -------------------------------------------------------------------------
function windowEdit(hObject,eventData)
% Called when the window-feature type uitable is edited.

oldPointer=pointerToWatch(gcbf);

handles = guidata(hObject);
fv=handles.featureVocabulary;  % a ref
%updatePFCategoryAmountForCurrentPF(handles);
rowIndex=eventData.Indices(1);
wfTypeNdx = rowIndex+1;  % +1 b/c default is not shown
handles.wfTypeNdx = wfTypeNdx;
wfType = fv.wfTypes{wfTypeNdx};
%handles.data{handles.pfNdx}.(wfType).enabled = eventData.NewData;

% Enable/disable the selected window-feature type in the model
% Unless the the window-feature type is 'default', which cannot be disabled
% except by disabling the per-frame feature
if ~isequal(wfType,'default')
  pfNdx=handles.pfNdx;
  pfName=fv.subdialectPFNames{pfNdx};
  fv.setWFTypeEnablement(pfName,wfType,eventData.NewData)
end

guidata(hObject,handles);
updateWindowTableRow(handles,rowIndex)
updateWinParamsAndEnablement(handles);
updatePFTableForCurrentPF(handles);
updateBasicTableAllCategoriesOfCurrentPF(handles);
% if eventData.NewData
%   updateWinParams(handles);
%   enableWinParams(handles);
% else
%   disableWinParams(handles);
% end
restorePointer(gcbf,oldPointer);
return


% % -------------------------------------------------------------------------
% function updateBasicTablePFCategoryLevel(handles,pfCategoryIndex)
% % Update the amount displayed for a particular per-frame feature category
% % in the basic table.
% fv=handles.featureVocabulary;
% pfCategoryName=fv.pfCategoryNames{pfCategoryIndex};
% basicData = get(handles.basicTable,'Data');
% basicData{pfCategoryIndex,2} = fv.getPFCategoryLevel(pfCategoryName);
% set(handles.basicTable,'Data',basicData);
% return


% -------------------------------------------------------------------------
function updateWinParams(handles)
% Update the GUI controls displaying window-feature parameters to match the current
% model state.

pfNdx=handles.pfNdx;
winNdx=handles.wfTypeNdx;
if isempty(pfNdx) || isempty(winNdx),
  set(handles.MinWindow,'String','');
  set(handles.MaxWindow,'String','');
  set(handles.WindowStep,'String','');
  set(handles.WindowOffsets,'String','');
  set(handles.TransNone,'Value',false);
  set(handles.TransFlip,'Value',false);
  set(handles.TransAbs,'Value',false);
  set(handles.TransRel,'Value',false);
  set(handles.ExtraParams,'String','');
  set(handles.extraParamStatic,'String','No extra params','FontAngle','italic');  
else
  fv=handles.featureVocabulary;
  wfType = fv.wfTypes{winNdx};
  wfParamsThisType = fv.vocabulary{pfNdx}.(wfType);
  set(handles.MinWindow,'String',num2str(wfParamsThisType.values.min_window_radius));
  set(handles.MaxWindow,'String',num2str(wfParamsThisType.values.max_window_radius));
  set(handles.WindowStep,'String',num2str(wfParamsThisType.values.nwindow_radii));
  set(handles.WindowOffsets,'String',num2str(wfParamsThisType.values.window_offsets));
  set(handles.TransNone,'Value',...
      any(strcmp('none',wfParamsThisType.values.trans_types)));
      %bitand(1,curParams.values.trans_types));
  set(handles.TransFlip,'Value',...
      any(strcmp('flip',wfParamsThisType.values.trans_types)));
      %bitand(4,curParams.values.trans_types));
  set(handles.TransAbs,'Value',...
      any(strcmp('abs',wfParamsThisType.values.trans_types)));
      %bitand(2,curParams.values.trans_types));
  set(handles.TransRel,'Value',...
      any(strcmp('relative',wfParamsThisType.values.trans_types)));
      %bitand(8,curParams.values.trans_types));
  if isfield(wfParamsThisType.values,fv.wfExtraParamNames{winNdx})
    extraParam = fv.wfExtraParamNames{winNdx};
    set(handles.ExtraParams, ...
        'String',wfParamsThisType.values.(extraParam));
    set(handles.extraParamStatic, ...
        'String',extraParam, ...
        'FontAngle','normal');
  else
    set(handles.ExtraParams,'String','');
    set(handles.extraParamStatic,'String','No extra params','FontAngle','italic');  
    % set(handles.ExtraParams,'Enable','off','String','--');
    % set(handles.extraParamStatic,'String','Feature Specific');
  end
end
return

% -------------------------------------------------------------------------
function updateWinParamsAndEnablement(handles)
% Update the GUI controls displaying window-feature parameters to match the current
% model state, and also the enablement/disablement of these controls.
updateWinParams(handles);
updateWinParamsEnablement(handles)
return


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
return


% % -------------------------------------------------------------------------
% function disableWindowTable(handles)
% set(handles.windowTable,'enable','off');
% set(handles.pushbutton_copy_windowtypes,'enable','off');
% disableWinParams(handles);
% 
% 
% % -------------------------------------------------------------------------
% function enableWindowTable(handles)
% set(handles.windowTable,'enable','on');
% set(handles.pushbutton_copy_windowtypes,'enable','on');


% -------------------------------------------------------------------------
function updateWindowTableEnablement(handles)
fv=handles.featureVocabulary;  % a ref
pfNdx=handles.pfNdx;
if isempty(pfNdx)
  enabled=false;
else
  enabled=fv.pfIsInVocabulary(pfNdx);
end
set(handles.windowTable,'enable',onIff(enabled));
set(handles.pushbutton_copy_windowtypes,'enable',onIff(enabled));
set(handles.popupmenu_copy_windowtypes,'enable',onIff(enabled));
return


% -------------------------------------------------------------------------
function updateWinParamsEnablement(handles)
pfNdx=handles.pfNdx;
wfTypeNdx=handles.wfTypeNdx;
if isempty(pfNdx) || isempty(wfTypeNdx),
  disableWinParams(handles);
else  
  fv=handles.featureVocabulary;   % a ref
  %pfInVocab=fv.pfIsInVocabulary(pfNdx);
  %wfTypeInVocab=fv.wfTypeIsInVocabulary(pfNdx,wfTypeNdx);
  if fv.pfIsInVocabulary(pfNdx) && ...
     fv.wfTypeIsInVocabulary(pfNdx,wfTypeNdx) ,
    enableWinParams(handles);
  else
    disableWinParams(handles);
  end
end
return


% -------------------------------------------------------------------------
function disableWinParams(handles)
% fprintf('In disableWinParams()\n');
set(handles.MinWindow,'enable','off');
set(handles.MaxWindow,'enable','off');
set(handles.WindowStep,'enable','off');
set(handles.WindowOffsets,'enable','off');
set(handles.TransNone,'enable','off');
set(handles.TransFlip,'enable','off');
set(handles.TransAbs,'enable','off');
set(handles.TransRel,'enable','off');
set(handles.ExtraParams,'enable','off');
set(handles.popupmenu_copy_windowparams,'enable','off');
set(handles.pushbutton_copy_windowparams,'enable','off');
return


% -------------------------------------------------------------------------
function enableWinParams(handles)
%defBack = get(handles.text5,'BackgroundColor');
%defFore = [0 0 0];
set(handles.MinWindow,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.MaxWindow,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.WindowStep,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);
set(handles.WindowOffsets,'enable','on');%,'BackgroundColor',defBack,'ForegroundColor',defFore);

fv=handles.featureVocabulary;   % a ref
if ~isempty(fv.wfExtraParamNames{handles.wfTypeNdx})
  set(handles.ExtraParams,'enable','on');
  %set(handles.extraParamStatic,'String',fv.wfExtraParamNames{handles.wfTypeNdx});
else
  set(handles.ExtraParams,'enable','off');
  %set(handles.extraParamStatic,'String','');
end

set(handles.TransNone,'enable','on');
set(handles.TransFlip,'enable','on');
set(handles.TransAbs,'enable','on');
set(handles.TransRel,'enable','on');
set(handles.popupmenu_copy_windowparams,'enable','on');
set(handles.pushbutton_copy_windowparams,'enable','on');
return


% -------------------------------------------------------------------------
function MinWindow_Callback(hObject, eventdata, handles)  %#ok
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
%wfType = fv.wfTypes{handles.wfTypeNdx};
%handles.data{handles.pfNdx}.(wfType).values.min_window_radius = curVal;

% set in featureVocabulary
pfName=fv.subdialectPFNames{handles.pfNdx};
wfType = FeatureVocabularyForSelectFeatures.wfTypes{handles.wfTypeNdx};
wfParamName='min_window_radius';
fv.setWFParam(pfName,wfType,wfParamName,curVal);

guidata(hObject,handles);
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function MinWindow_CreateFcn(hObject, eventdata, handles)  %#ok
% hObject    handle to MinWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
function MaxWindow_Callback(hObject, eventdata, handles)  %#ok
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
wfType = fv.wfTypes{handles.wfTypeNdx};
%handles.data{handles.pfNdx}.(wfType).values.max_window_radius = curVal;

% set in featureVocabulary
pfName=fv.subdialectPFNames{handles.pfNdx};
%wfType = FeatureVocabularyForSelectFeatures.wfTypes{handles.wfTypeNdx};
wfParamName='max_window_radius';
fv.setWFParam(pfName,wfType,wfParamName,curVal);

guidata(hObject,handles);
updateMaxWindowRadiusEditBox(handles);
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function MaxWindow_CreateFcn(hObject, eventdata, handles)  %#ok
% hObject    handle to MaxWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
function WindowStep_Callback(hObject, eventdata, handles)  %#ok
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
wfType = fv.wfTypes{handles.wfTypeNdx};
%handles.data{handles.pfNdx}.(wfType).values.nwindow_radii = curVal;

% set in featureVocabulary
pfName=fv.subdialectPFNames{handles.pfNdx};
%wfType = FeatureVocabularyForSelectFeatures.wfTypes{handles.wfTypeNdx};
wfParamName='nwindow_radii';
fv.setWFParam(pfName,wfType,wfParamName,curVal);

guidata(hObject,handles);
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function WindowStep_CreateFcn(hObject, eventdata, handles)  %#ok
% hObject    handle to WindowStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return


% -------------------------------------------------------------------------
function WindowOffsets_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to WindowOffsets (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of WindowOffsets as text
%        str2double(get(hObject,'String')) returns contents of WindowOffsets as a double
curVal = str2double(get(hObject,'String'));
if isempty(curVal)
  msgbox('Enter numerical values. eg: "-1 0 1" (without with quotes)');
end
handles = guidata(hObject);
fv=handles.featureVocabulary;
wfType = fv.wfTypes{handles.wfTypeNdx};
%handles.data{handles.pfNdx}.(wfType).values.window_offsets = curVal;

% set in featureVocabulary
pfName=fv.subdialectPFNames{handles.pfNdx};
%wfType = FeatureVocabularyForSelectFeatures.wfTypes{handles.wfTypeNdx};
wfParamName='window_offsets';
fv.setWFParam(pfName,wfType,wfParamName,curVal);

guidata(hObject,handles);
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function WindowOffsets_CreateFcn(hObject, eventdata, handles)  %#ok
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
function TransNone_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to TransNone (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransNone
fv=handles.featureVocabulary;
%curFn = fv.wfTypes{handles.wfTypeNdx};
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
pfName=fv.subdialectPFNames{handles.pfNdx};
wfType = fv.wfTypes{handles.wfTypeNdx};
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
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in TransFlip.
function TransFlip_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to TransFlip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransFlip
fv=handles.featureVocabulary;
%curFn = fv.wfTypes{handles.wfTypeNdx};
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
wfType = fv.wfTypes{handles.wfTypeNdx};
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
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in TransAbs.
function TransAbs_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to TransAbs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransAbs

fv=handles.featureVocabulary;
%curFn = fv.wfTypes{handles.wfTypeNdx};
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
wfType = fv.wfTypes{handles.wfTypeNdx};
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
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in TransRel.
function TransRel_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to TransRel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of TransRel

fv=handles.featureVocabulary;
%curFn = fv.wfTypes{handles.wfTypeNdx};
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
wfType = fv.wfTypes{handles.wfTypeNdx};
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
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in push_cancel.
function push_cancel_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to push_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%uiresume(handles.figure_SelectFeatures);
if handles.doneWithOpeningFunction ,
  delete(handles.figure_SelectFeatures);
end
return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
%basicData = get(handles.basicTable,'Data');
fv=handles.featureVocabulary;
[thereIsConsensus,maxWindowRadiusConsensus]= ...
  fv.getConsensusMaxWindowRadiusForAllWFs();
if ~thereIsConsensus ,
  maxWindowRadiusConsensus=[];
end
%featureWindowSize = str2double(get(handles.editSize,'String'));
%windowFeaturesParams = convertData(handles);
windowFeaturesParams = handles.featureVocabulary.getInJLabelDataFormat();
set(handles.figure_SelectFeatures,'Visible','off');
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
       maxWindowRadiusConsensus);
return


% -------------------------------------------------------------------------
function ExtraParams_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to ExtraParams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ExtraParams as text
%        str2double(get(hObject,'String')) returns contents of ExtraParams as a double

fv=handles.featureVocabulary;
if isempty(fv.wfExtraParamNames{handles.wfTypeNdx}), return; end

str = get(hObject,'String');
str = strtrim(str);
newValue = str2double(str);
%extraParam = fv.wfExtraParamNames{handles.wfTypeNdx};
wfType = fv.wfTypes{handles.wfTypeNdx};
%handles.data{handles.pfNdx}.(wfType).values.(extraParam) = newValue;

% set in featureVocabulary
pfName=fv.subdialectPFNames{handles.pfNdx};
%wfType = handles.featureVocabulary.wfTypes{handles.wfTypeNdx};
wfParamName=fv.wfExtraParamNames{handles.wfTypeNdx};
fv.setWFParam(pfName,wfType,wfParamName,newValue);

guidata(hObject,handles);
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function ExtraParams_CreateFcn(hObject, eventdata, handles)  %#ok
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
function Save_Callback(hObject, eventdata, handles)  %#ok
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
function Load_Callback(hObject, eventdata, handles)  %#ok
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
function popupmenu_copy_windowparams_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to popupmenu_copy_windowparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_copy_windowparams contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_copy_windowparams


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function popupmenu_copy_windowparams_CreateFcn(hObject, eventdata, handles)  %#ok
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
function pushbutton_copy_windowparams_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to pushbutton_copy_windowparams (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(handles.pfNdx) || isempty(handles.wfTypeNdx),
  return;
end
fv=handles.featureVocabulary;
windowComp = fv.wfTypes;

% which perframe fn are we on?
pfNdx = handles.pfNdx;

% which window fn are we copying to?
winNdxTo = handles.wfTypeNdx;

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
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
return


% -------------------------------------------------------------------------
% --- Executes on selection change in popupmenu_copy_windowtypes.
function popupmenu_copy_windowtypes_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to popupmenu_copy_windowtypes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_copy_windowtypes contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_copy_windowtypes


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function popupmenu_copy_windowtypes_CreateFcn(hObject, eventdata, handles)  %#ok
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
pfNameList = fv.subdialectPFNames;
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
updateWindowTableEnablement(handles);
%enableWindowTable(handles);
updateWinParams(handles);
updateWinParamsEnablement(handles);
% curFn = fv.wfTypes{handles.wfTypeNdx};
% if handles.data{handles.pfNdx}.(curFn).valid,
%   enableWinParams(handles);
% end

guidata(hObject,handles);
updateBasicTableAllCategoriesOfCurrentPF(handles);
updatePFTableForCurrentPF(handles);
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

oldPointer=pointerToWatch(gcbf);

if eventdata.NewValue == handles.togglebutton_tabdescription,
  handles.currentTab = 'description';
elseif eventdata.NewValue == handles.togglebutton_tabperframehistogram,
  handles.currentTab = 'perframehistogram';
end

handles = UpdateDescriptionPanels(handles);
guidata(hObject,handles);

restorePointer(gcbf,oldPointer);

return


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
      HistogramPerFrameFeature(handles.jld,fv.subdialectPFNames{handles.pfNdx},...
                               'axes',handles.axes_histogram,...
                               'unknowncolor','w',...
                               'labelcolors',jet(handles.jld.nbehaviors)*.7);
    handles.histogramData.perframe_idx(i) = handles.pfNdx;
  else
    [handles.histogramData.hhist,...
     ~,~,hleg,hxlabel,hylabel,...
     handles.histogramData.frac{i},handles.histogramData.frac_outside{i},...
     handles.histogramData.edges{i},handles.histogramData.centers_plot{i}] = ...
      HistogramPerFrameFeature(handles.jld,fv.subdialectPFNames{handles.pfNdx},...
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
function editSize_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to editSize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSize as text
%        str2double(get(hObject,'String')) returns contents of editSize as a double

oldPointer=pointerToWatch(gcbf);

curVal = str2double(get(hObject,'String'));
if isempty(curVal) || (round(curVal)-curVal)~=0 || curVal<=0 ,
  updateMaxWindowRadiusEditBox(handles);
  restorePointer(gcbf,oldPointer);
  return
end

% Set the max window radius in all the window features, _and_ 
% in all the WF amount presets.  This will mean that any future selection
% of a preset amount with get the new value of max_window_radius.
fv=handles.featureVocabulary;
fv.setMaxWindowRadiusForAllWFs(curVal);
fv.setMaxWindowRadiusForAllWFAmounts(curVal);

% Update the view
updateMaxWindowRadiusEditBox(handles);
updateWinParams(handles);
drawnow('update');  % want to see that number change in window params pronto!
updateBasicTable(handles);
updatePFTable(handles);

restorePointer(gcbf,oldPointer);

return


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function editSize_CreateFcn(hObject, eventdata, handles)  %#ok
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
function togglebuttonMode_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to togglebuttonMode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if get(hObject,'Value')
  handles=setViewMode(handles,'advanced');
else
  handles=setViewMode(handles,'basic');
end
guidata(handles.figure_SelectFeatures,handles);
return


% -------------------------------------------------------------------------
function handles=setViewMode(handles,newMode)
% Set the mode of the view (to 'basic' or 'advanced')

curLoc = get(handles.figure_SelectFeatures,'Position');
if isequal(newMode,'advanced')
  handles.mode = 'advanced';
  set(handles.figure_SelectFeatures,'Position',[curLoc(1:2) handles.advancedSize]);
  set(handles.togglebuttonMode,'String','Basic <');
else
  handles.mode = 'basic';
  set(handles.figure_SelectFeatures,'Position',[curLoc(1:2) handles.basicSize]);  
  set(handles.togglebuttonMode,'String','Advanced >');
end
return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_hist.
function pushbutton_hist_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to pushbutton_hist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prcEdges = [5 15 30 50 70 85 95];

fv=handles.featureVocabulary;
histWFTypeNdx = find(strcmp('hist',fv.wfTypes));
histExtraParamName = fv.wfExtraParamNames{histWFTypeNdx};  %#ok

h = waitbar(0,'Computing hist bins');
for pfIndex = 1:numel(fv.subdialectPFNames)
  pfName = fv.subdialectPFNames{pfIndex};
  waitbar(pfIndex/numel(fv.subdialectPFNames),h);
  
  if ~fv.pfIsInVocabulary(pfIndex) || ~fv.wfTypeIsInVocabulary(pfIndex,'hist')
    continue;
  end
  
  allData = [];
  for expi = 1:handles.jld.nexps,
    
    % load per-frame data for this experiment
    perframeDirName = handles.jld.GetFile('perframedir',expi);
    file = fullfile(perframeDirName,[pfName,'.mat']);
    if ~exist(file,'file'),
      warning('Per-frame data file %s does not exist',file);
      continue;
    end
    
    perframedata = load(file);
    
    for fly = 1:handles.jld.nflies_per_exp(expi),
      
      x = perframedata.data{fly};
      allData = [allData ; x(:)];  %#ok
    end
  end
  bins = prctile(allData,prcEdges);
  minD = min(allData);
  maxD = max(allData);
  binMin = minD - (maxD-minD);
  binMax = maxD + (maxD-minD);
  bins = [binMin bins binMax];  %#ok
  %handles.data{ndx}.hist.values.(histExtraName) = bins;
  fv.setWFParam(pfName,'hist',histExtraParamName,bins)
end
delete(h);
guidata(hObject,handles);
updateWindowTableAndEnablement(handles);
updateWinParamsAndEnablement(handles);
% if ~isempty(handles.pfNdx)
%   updateWindowTable(handles);
%   enableWindowTable(handles);
% end
return


% -------------------------------------------------------------------------
function figure_SelectFeatures_CloseRequestFcn(hObject,eventdata,handles)  %#ok
%fprintf('In CloseRequestFcn()\n');
push_cancel_Callback(hObject,eventdata,handles);
return


% % -------------------------------------------------------------------------
% % --- Executes on button press in pushbutton_done.
% function pushbutton_ok_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton_done (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % configfile = handles.jld.configfilename;
% % [~,~,ext ] = fileparts(configfile);
% % if strcmpi(ext,'.xml'),
% %   configparams = ReadXMLParams(configfile);
% % elseif strcmpi(ext,'.mat'),
% %   configparams = load(configfile);
% % else
% %   errordlg('Project file %s is invalid. Cannot save the window features',configfile);
% % end
% % if ~isfield(configparams.file,'featureparamfilename') || isempty(configparams.file.featureparamfilename)
% %   behaviorname = configparams.behaviors.names;
% %   if iscell(behaviorname),
% %     defaultname = sprintf('WindowFeatures_%s.xml',behaviorname{1});
% %   else
% %     defaultname = sprintf('WindowFeatures_%s.xml',behaviorname);
% %   end
% %   [fname,fpath]= uiputfile(fullfile('params','*.xml'),'Enter a name for feature config file',defaultname);
% %   if isempty(fname),return, end
% %   featureconfigfile = fullfile(fpath,fname);
% %   configparams.file.featureparamfilename = featureconfigfile;
% %   docNode = com.mathworks.xml.XMLUtils.createDocument('params');
% %   toc = docNode.getDocumentElement;
% %   fnames = fieldnames(configparams);
% %   for ndx = 1:numel(fnames)
% %     toc.appendChild(createXMLNode(docNode,fnames{ndx},configparams.(fnames{ndx})));
% %   end
% %   xmlwrite(configfile,docNode);
% % end
% % 
% % featureconfigfile = configparams.file.featureparamfilename;
% % params = convertData(handles);
% % basicData = get(handles.basicTable,'Data');
% % featureWindowSize = round(str2double(get(handles.editSize,'String')));
% % docNode = createParamsXML(params,basicData,featureWindowSize);
% % xmlwrite(featureconfigfile,docNode);
% 
% pushbutton_done_Callback(hObject,eventdata,handles);
% push_cancel_Callback(hObject,eventdata,handles);
% return
