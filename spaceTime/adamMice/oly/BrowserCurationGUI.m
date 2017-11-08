function varargout = BrowserCurationGUI(varargin)
% BROWSERCURATIONGUI MATLAB code for BrowserCurationGUI.fig
%      BROWSERCURATIONGUI, by itself, creates a new BROWSERCURATIONGUI or raises the existing
%      singleton*.
%
%      H = BROWSERCURATIONGUI returns the handle to a new BROWSERCURATIONGUI or the handle to
%      the existing singleton*.
%
%      BROWSERCURATIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in BROWSERCURATIONGUI.M with the given input arguments.
%
%      BROWSERCURATIONGUI('Property','Value',...) creates a new BROWSERCURATIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before BrowserCurationGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to BrowserCurationGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help BrowserCurationGUI

% Last Modified by GUIDE v2.5 20-Jul-2011 14:25:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @BrowserCurationGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @BrowserCurationGUI_OutputFcn, ...
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

% BrowserCurationGUI(browserObj,curationCfg)
% browserObj: handle to OlyDat.Browser object
% curationCfg: struct array providing curation fields. curationCfg
% configures the ith column in the curation table. Note that the 1st column
% is reserved for the "excluded" column. curationCfg fields:
%   * Name: curation field name
%   * TableName: name for column header in table
%   * TableColWidth: width in pixels for col. 0 indicates auto-size.
function BrowserCurationGUI_OpeningFcn(hObject, eventdata, handles, varargin)

handles.output = hObject;
handles.browser = varargin{1};
handles.curationCfg = varargin{2};
handles.selectedCell = [];
 
% table metadata
handles.tableMetadataExcludedColIdx = 1;
handles.tableMetadataCurationField2Col = struct();
handles.tableMetadataCol2Field = {''}; % first col is empty for excluded col

NcurationCfg = numel(handles.curationCfg);
for cfgIdx = 1:NcurationCfg
    name = handles.curationCfg(cfgIdx).Name;
    handles.tableMetadataCurationField2Col.(name) = cfgIdx+1;
    handles.tableMetadataCol2Field{cfgIdx+1,1} = name;
end    

% configure table based on tableMetadata
Ncol = NcurationCfg+1;
tfColEditable = false(1,Ncol);
tfColEditable(handles.tableMetadataExcludedColIdx) = true; % excluded col is editable; rest are not
colFmt = cell(1,Ncol);
colFmt{1} = 'logical'; % unused, see below
colNames = [{'Excluded'} {handles.curationCfg.TableName}];
colWidths = [{60} {handles.curationCfg.TableColWidth}];
for colIdx = 1:Ncol
    if colWidths{colIdx}==0
        colWidths{colIdx} = 'auto';
    end
end

% (*) AL 1/5/2012: including ColumnFormat in the set leads to a SEGV on
% maci64. This is consistent/reproduceable. However simply leaving out the
% ColumnFormat set seems to work fine. The ColumnFormat property will not
% have the right size (eg will not have Ncol number of elements), but
% in this situation MATLAB appears to just use auto-formatting for all
% columns and for now that is working fine.
set(handles.tblCuration,'ColumnEditable',tfColEditable,... % 'ColumnFormat',colFmt,... (*) see above
                        'ColumnName',colNames,...
                        'ColumnWidth',colWidths,...
                        'Data',cell(1,Ncol));
                    
% context menu
uh = uicontextmenu;
uimenu(uh,'label','Enter Curation Info','Callback',@ctxtMnuEnterCurationInfo);
uimenu(uh,'label','Clear Curation Info','Callback',@ctxtMnuClearCurationInfo);
uimenu(uh,'label','Sort Ascending','Callback',@ctxtMnuSortAscending,'Separator','on');
uimenu(uh,'label','Sort Descending','Callback',@ctxtMnuSortDescending);
uimenu(uh,'label','Bring Excluded/Marked to Top','Callback',@ctxtMnuBringExcludedMarkedToTop);
set(handles.tblCuration,'UIContextMenu',uh);

% handles.tableMetadataCurationField2Col = struct('experiment_id',2,'experiment_name',3,'manual_pf',4,
% 'notes_curation',5,'flag_redo',6,'flag_review',7,'manual_curation_date',8,'manual_curator',9);
% handles.tableMetadataCol2Field = {'';'experiment_id';'experiment_name';'manual_pf';'notes_curation';
%     'flag_redo';'flag_review';'manual_curation_date';'manual_curator'};

guidata(hObject,handles);

% This gui is custom-resizeable
set(handles.uipanel1,'Units','Normalized');

% Explicitly call resizefcn b/c the orig layout is not precisely consistent
% with what is done in resize.
uipanel1_ResizeFcn(handles.uipanel1,[],handles);

centerfig(hObject);

function uipanel1_ResizeFcn(hObject,~,handles)
TBLCURATION_Y_TOP_SPACER = 8.0;
TBLCURATION_Y_BOTTOM_SPACER = 10.0;
TBLCURATION_X_SPACER = 8.0;
TBL_CURATION_BETWEEN_BUTTONS_X_SPACER = 3.0;
TBL_CURATION_AFTER_BUTTONS_X_SPACER = 11.0;

% get panel width/height in pixels
origUnits = get(hObject,'Units');
set(hObject,'Units','pixels');
pnlPos = get(hObject,'Position');
set(hObject,'Units',origUnits);
pnlWidth = pnlPos(3);
pnlHeight = pnlPos(4);

% get top posn of unexclude all button
pbUnexcludeAllPos = get(handles.pbUnexcludeAll,'Position');
pbTopPosInPixels = pbUnexcludeAllPos(2)+pbUnexcludeAllPos(4);

% set tbl position
tblPos(1) = TBLCURATION_X_SPACER;
tblPos(2) = pbTopPosInPixels + TBLCURATION_Y_BOTTOM_SPACER;
tblPos(3) = pnlWidth-2.0*TBLCURATION_X_SPACER;
tblPos(4) = pnlHeight-tblPos(2)-TBLCURATION_Y_TOP_SPACER;
set(handles.tblCuration,'Position',tblPos);

% set posn of UnexcludeAll, MarkAllCurated buttons
pbMarkAllUncuratedPos = get(handles.pbMarkAllUncuratedAsPass,'Position');
pbUnexcludeAllXPos = pnlWidth - pbUnexcludeAllPos(3) - pbMarkAllUncuratedPos(3) ...
    - TBL_CURATION_BETWEEN_BUTTONS_X_SPACER - TBL_CURATION_AFTER_BUTTONS_X_SPACER;
pbUnexcludeAllPos(1) = pbUnexcludeAllXPos;
pbMarkAllUncuratedXPos = pnlWidth - pbMarkAllUncuratedPos(3) - TBL_CURATION_AFTER_BUTTONS_X_SPACER;
pbMarkAllUncuratedPos(1) = pbMarkAllUncuratedXPos;
set(handles.pbUnexcludeAll,'Position',pbUnexcludeAllPos);
set(handles.pbMarkAllUncuratedAsPass,'Position',pbMarkAllUncuratedPos);

function figure1_ResizeFcn(~,~,~)
%none 

function varargout = BrowserCurationGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

% A note on exp selection. When rows of the CurationGUI are selected, the
% corresponding experiments are selected in the Browser. However, when
% experiments are selected in the Browser, the corresponding rows of the
% curationGUI are NOT selected. The primary reason for this is that the
% native ML tables do not permit programmatically setting their selection.

function tblCuration_CellSelectionCallback(hObject,eventdata,handles)
if ~isempty(eventdata.Indices)
    handles.selectedCell = eventdata.Indices;
    guidata(hObject,handles);    
    zlclSelectExperimentsInBrowserBasedOnSelectedCellCache(handles);    
else
    handles.selectedCell = [];
    guidata(hObject,handles);
end

function tblCuration_CellEditCallback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
% If there are multiple selections in the table, the upper-left-most will be used
rowIdx = eventdata.Indices(1,1); 
colIdx = eventdata.Indices(1,2);
data = get(hObject,'Data');

switch colIdx
    case handles.tableMetadataExcludedColIdx
        % 'excluded' column
        
        % Make sure this experiment (row) is selected. This may be
        % unnecessary, ie if tblCuration_CellSelectionCallback is always
        % called first.
        expID = data{rowIdx,handles.tableMetadataCurationField2Col.experiment_id};
        handles.browser.selectExperiment(expID);
        
        tf = eventdata.NewData;
        handles.browser.setRemoveExp(tf,false);
        
        parentFigure = ancestor(hObject,'figure');
        uistack(parentFigure,'top');
    otherwise
        assert(false,'This no longer ever occurs');
end

% This is unused
% function s = zlclRow2CurationInfoStruct(row,handles)
% flds = handles.tableMetadataCol2Field;
% Nflds = numel(flds);
% assert(numel(row)==Nflds);
% s = struct();
% for c = 1:Nflds
%     if ~isempty(flds{c}) % excluded col will have an empty field
%         val = row{c};
%         if ischar(val)
%             % convert empty option ' ' to ''
%             val = strtrim(val);
%         end
%         s.(flds{c}) = val;
%     end
% end

%%% ACTIONS/HELPERS
function zlclSelectExperimentsInBrowserBasedOnSelectedCellCache(handles)
rowIdxs = handles.selectedCell(:,1);
data = get(handles.tblCuration,'Data');
eids = data(rowIdxs,handles.tableMetadataCurationField2Col.experiment_id);
handles.browser.selectExperiment(cell2mat(eids));

function zlclDisplaySelectCellMsgBox(msg)
msgbox(msg,'No cell selected');

function zlclActionEnterCurationMarks(handles)
if isempty(handles.selectedCell)
    zlclDisplaySelectCellMsgBox('Please select a cell to specify a row.');
else
    % select experiments in the Browser and open CurationInfoEntry.
    zlclSelectExperimentsInBrowserBasedOnSelectedCellCache(handles);
    handles.browser.curationOpenCurationInfoEntry();
end

function zlclActionClearCurationMarks(handles)
if isempty(handles.selectedCell)
    zlclDisplaySelectCellMsgBox('Please select a cell to specify a row.');
else
    zlclSelectExperimentsInBrowserBasedOnSelectedCellCache(handles);
    handles.browser.curationClearCurationInfoForCurrentExp();
end

function zlclSortByData(tableH,sortData,dir)
data = get(tableH,'Data');
assert(size(data,1)==numel(sortData(:)));

[~,idx] = sort(sortData(:));
switch dir
    case 'ascend'
        %none
    case 'descend'
        idx = idx(end:-1:1);
end

data = data(idx,:);
set(tableH,'Data',data);

function zlclSortHelper(handles,direction)
if isempty(handles.selectedCell)
    zlclDisplaySelectCellMsgBox('Select a cell to specify a column for sorting.');
else
    data = get(handles.tblCuration,'Data');
    colIdx = handles.selectedCell(1,2); % use column of first selected cell
    colData = data(:,colIdx);
    if isnumeric(colData{1}) || islogical(colData{1})
        colData = cell2mat(colData);
    end
    zlclSortByData(handles.tblCuration,colData,direction);
    
    % After setting the table data to a new value, MATLAB de-selects any
    % selected cells. Ideally, would like to have the original cell (in
    % terms of coordinates) selected, so that the user may eg 'sort
    % ascending', then immediately 'sort descending' if they got the wrong
    % direction. Unfortunately, this is not possible with the ML uitable
    % without dropping into Java; this is a potentialy maintenance
    % dependency so for now we don't do it. To sort one way and then the
    % other, the user will have to re-select a cell in the appropriate
    % column.
end

% At the moment, this checks for marking by checking the manual_curator
% field (hardcoded).
function zlclBringExcludedMarkedToTop(handles)
data = get(handles.tblCuration,'Data');
Nrows = size(data,1);

excludedColIdx = handles.tableMetadataExcludedColIdx;
if isfield(handles.tableMetadataCurationField2Col,'manual_curator')
  % todo: this should not be hardcoded against manual_curator
  manualCuratorColIdx = handles.tableMetadataCurationField2Col.manual_curator; 
  tfManualCuration = true;
else
  tfManualCuration = false;
end

indicator = false(Nrows,1);
for c = 1:Nrows
    row = data(c,:);
    
    tfExcluded = row{excludedColIdx};
    if tfExcluded
        indicator(c) = true;
        continue;
    end
    
    if tfManualCuration && ~isempty(row{manualCuratorColIdx})
        indicator(c) = true;
    end
end
zlclSortByData(handles.tblCuration,indicator,'descend');

%%% BUTTONS
function pbEnterCurationMarks_Callback(hObject,eventdata,handles)
zlclActionEnterCurationMarks(handles);

function pbClearMarks_Callback(hObject,eventdata,handles)
zlclActionClearCurationMarks(handles);

function pbUnexcludeAll_Callback(hObject, eventdata, handles)
handles.browser.unexcludeAllExps();

function pbMarkAllUncuratedAsPass_Callback(hObject, eventdata, handles)
handles.browser.curationMarkAllUncuratedAsPass();

%%% CONTEXT MENU
function ctxtMnuEnterCurationInfo(src,~)
zlclActionEnterCurationMarks(guidata(src));

function ctxtMnuClearCurationInfo(src,~)
zlclActionClearCurationMarks(guidata(src));

function ctxtMnuSortAscending(src,~)
zlclSortHelper(guidata(src),'ascend');

function ctxtMnuSortDescending(src,~)
zlclSortHelper(guidata(src),'descend');

function ctxtMnuBringExcludedMarkedToTop(src,~)
zlclBringExcludedMarkedToTop(guidata(src));

%%% REGULAR MENU
function mnuFileLoad_Callback(hObject, eventdata, handles)
handles.browser.curationLoadCurationFile();

function mnuFileAppend_Callback(hObject, eventdata, handles)
handles.browser.curationAppendCurationDataToFile();

function mnuFileWriteNew_Callback(hObject, eventdata, handles)
handles.browser.curationWriteCurationDataToNewFile();

function mnuSortAscending_Callback(hObject,eventdata,handles)
zlclSortHelper(handles,'ascend');

function mnuSortDescending_Callback(hObject,eventdata,handles)
zlclSortHelper(handles,'descend');

function mnuBringExcludedMarkedToTop_Callback(hObject,eventdata,handles)
zlclBringExcludedMarkedToTop(handles);
