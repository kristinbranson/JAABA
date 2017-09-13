function varargout = BrowserGUI(varargin)
%BrowserGUI GUI code for OlyDat Browser
%
%   See also OlyDat.Browser.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
if ismac
    gui_Name = 'BrowserGUI_Mac';
else
    gui_Name = 'BrowserGUI_PC'; % use this for both PC and Unix
end
gui_State = struct('gui_Name',       gui_Name, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @browserGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @browserGUI_OutputFcn, ...
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

% browserGUI = BrowserGUI(browserObj)
function browserGUI_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
handles.output = hObject;
handles.browser = varargin{1};
guidata(hObject, handles);
centerGUI(hObject);

function varargout = browserGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function pumAuxVar_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
handles.browser.auxVarChanged();

function lbBStats_Callback(hObject, eventdata, handles)
handles.browser.statChanged();

function pumBStats_Callback(hObject, eventdata, handles)
handles.browser.statChanged();

function lbGroupVars_Callback(hObject, eventdata, handles)
handles.browser.groupOrGroup2SelectionChanged();

function pumGroup_Callback(hObject, eventdata, handles)
handles.browser.groupChanged();

function lbGroupVars2_Callback(hObject, eventdata, handles)
handles.browser.groupOrGroup2SelectionChanged();

function pumGroup2_Callback(hObject, eventdata, handles)
handles.browser.group2Changed();

function pbUseAllGroupVals_Callback(hObject, eventdata, handles)
strs = get(handles.lbGroupVars,'String');
set(handles.lbGroupVars,'Value',1:numel(strs));
handles.browser.groupOrGroup2SelectionChanged();

function pbUseAllGroupVals2_Callback(hObject, eventdata, handles)
strs = get(handles.lbGroupVars2,'String');
set(handles.lbGroupVars2,'Value',1:numel(strs));
handles.browser.groupOrGroup2SelectionChanged();

function pbPlot_Callback(hObject, eventdata, handles)
handles.browser.plot();

function pumPlotType_Callback(hObject, eventdata, handles)
handles.browser.newPlot();

function cbRemoveExp_Callback(hObject, eventdata, handles)
handles.browser.setRemoveExp(get(hObject,'Value'));

function pbCompare_Callback(hObject, eventdata, handles)
handles.browser.addComparePlot();



function figure1_CloseRequestFcn(hObject, eventdata, handles)
handles.browser.browserGUIWantsToClose();

function pbEnterCurationInfo_Callback(hObject, eventdata, handles)
handles.browser.curationOpenCurationInfoEntry();

function pbNextExpDetail_Callback(hObject, eventdata, handles)
handles.browser.nextExpDetail;

function pbPrevExpDetail_Callback(hObject, eventdata, handles)
handles.browser.prevExpDetail;

function pbCreateCustomGroup_Callback(hObject, eventdata, handles)
grpname = inputdlg('Specify a name for the group');
if ~isempty(grpname)
    grpname = grpname{1};
    handles.browser.createCustomGroup(grpname);
end

function pbLineReport_Callback(hObject, eventdata, handles)
handles.browser.openLineReport();

function pbExpDetail_Callback(hObject, eventdata, handles)
handles.browser.openExperimentDetailPlots();

function pbFileSystem_Callback(hObject, eventdata, handles)
handles.browser.openExperimentFileSystem();

% MENU
function mnuCuration_Load_Callback(hObject, eventdata, handles)
handles.browser.curationLoadCurationFile();

function mnuCuration_Append_Callback(hObject, eventdata, handles)
handles.browser.curationAppendCurationDataToFile();

function mnuCuration_Write_Callback(hObject, eventdata, handles)
handles.browser.curationWriteCurationDataToNewFile();

function mnuCuration_OpenCurationManager_Callback(hObject, eventdata, handles)
handles.browser.curationOpenCurationManager();

function pbExport_Callback(hObject, eventdata, handles)
handles.browser.exportPlot();
