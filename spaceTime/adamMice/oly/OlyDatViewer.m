function varargout = OlyDatViewer(varargin)
% OLYDATVIEWER MATLAB code for OlyDatViewer.fig
%      OLYDATVIEWER, by itself, creates a new OLYDATVIEWER or raises the existing
%      singleton*.
%
%      H = OLYDATVIEWER returns the handle to a new OLYDATVIEWER or the handle to
%      the existing singleton*.
%
%      OLYDATVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in OLYDATVIEWER.M with the given input arguments.
%
%      OLYDATVIEWER('Property','Value',...) creates a new OLYDATVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before OlyDatViewer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to OlyDatViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help OlyDatViewer

% Last Modified by GUIDE v2.5 19-Feb-2011 23:56:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @OlyDatViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @OlyDatViewer_OutputFcn, ...
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


function OlyDatViewer_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.analyzer = OlyDat.Analyzer(varargin{1});

% hardcode box dataset info
dsFam = SAGE.Lab('olympiad').assay('box');
dsNames = {dsFam.dataSets.name}';
Nds = numel(dsNames);
fldObjs = [];
for c = 1:Nds
    ds = dsFam.dataSet(dsNames{c});
    fldObjs = [fldObjs ds.fields]; %#ok<AGROW>
end
fldObjs = fldObjs(:);
fldNames = {fldObjs.name}';
[~,i] = unique(fldNames);
fldObjs = fldObjs(i);
handles.dataSetFields = fldObjs;
for c = 1:numel(fldObjs)
    handles.fieldName2Obj.(fldObjs(c).name) = fldObjs(c);
end

handles.queryComparisons = {'is equal to';'is greater than';'is less than'};

% pull
handles.currentQry = 'New Query';
set(handles.pnlQry,'Title',sprintf('Query: %s','New Query'));
handles.analyzer.DataVarName = 'data';
handles.qcp = QueryComponentPanel(handles.pnlOuter,{fldObjs.name}',handles.queryComparisons);

% anls
set(handles.lbAnalyses,'String',handles.analyzer.AnalysisNames);
zlclRefreshAnlsDesc(handles);
set(handles.edOutputDir,'Enable','off');

% output


% Update handles structure
guidata(hObject, handles);

%% QC
% function zlclRmQC(qc,handles,qcIdx)
% 
% if numel(handles.qc)==1
%     return;
% end
% 
% QCPANELDY = .2;
% 
% delete(handles.qc(qcIdx));
% handles.qc(qcIdx) = [];
% delete(handles.pbKillQC(qcIdx));
% handles.pbKillQC(qcIdx) = [];
% 
% % adjust posn of qcs
% for c = qcIdx:numel(handles.qc)
%     posn = get(handles.qc(c),'Position');
%     posn(2) = posn(2)+handles.qcHeight+QCPANELDY;
%     set(handles.qc(c),'Position',posn);
%     
%     posn = get(handles.pbKillQC(c),'Position');
%     posn(2) = posn(2)+handles.qcHeight+QCPANELDY;
%     set(handles.pbKillQC(c),'Position',posn);
% end
% 
% % add QC button
% posn = get(handles.pbAddQry,'Position');
% posn(2) = posn(2)+handles.qcHeight+QCPANELDY;
% set(handles.pbAddQry,'Position',posn);
% 
% guidata(handles.ax,handles);
% 
% Set compare mode: 'enumerated' or 'unenumerated'
% function zlclQCSetCompareMode(qc,mode,vals)
% rangeLB = findobj(qc,'Tag','lbQCRange');
% comparePUM = findobj(qc,'Tag','pumQryCmp');
% set(qc,'UserData',mode);
% switch mode
%     case 'enumerated'
%         set(comparePUM,'String',{'takes value(s)' 'equals' 'is less than' 'is greater than'});
%         set(comparePUM,'Value',1);
%         set(rangeLB,'String',vals);
%         zlclQCSetRangeLBOrEditB(qc,'range');
%     case 'unenumerated'
%         set(comparePUM,'String',{'equals' 'is less than' 'is greater than'});
%         set(rangeLB,'String',{});
%         zlclQCSetRangeLBOrEditB(qc,'edit');
% end

% % This just toggles the visibility of rangeLB/editLB.
% function zlclQCSetRangeLBOrEditB(qc,mode)
% rangeLB = findobj(qc,'Tag','lbQCRange');
% editB = findobj(qc,'Tag','edQryVal');
% switch mode
%     case 'range'
%         set(rangeLB,'Visible','on');
%         set(editB,'Visible','off');
%     case 'edit'
%         set(rangeLB,'Visible','off');
%         set(editB,'Visible','on');
% end

%%

% function qry = zlclParseQueryTxt(handles)
% qrytxt = get(handles.edExpQry,'String');
% qrytxt = cellstr(qrytxt);
% Nqry = numel(qrytxt);
% qrycell = cell(Nqry,1);
% pat = '(\w+)[=<>]([\*\w]+)';
% for c = 1:Nqry
%     toks = regexp(qrytxt{c},pat,'tokens');
%     toks = toks{1};
%     qrycell{c} = SAGE.Query.Compare(toks{1},'=',toks{2});
% end
% qry = SAGE.Query.All(qrycell{:});

function zlclRefreshAnlsDesc(handles)
name = zlclCurrentAnalysisName(handles);
anls = handles.analyzer.getAnalysisByName(name);
desc = anls.Description;
set(handles.edAnalysisDesc,'String',desc);

function name = zlclCurrentAnalysisName(handles)
v = get(handles.lbAnalyses,'Value');
str = get(handles.lbAnalyses,'String');
name = str{v};


% function zlclRefreshAnalyzerDataVar(handles)
% handles.analyzer.DataVarName = get(handles.edDataVar,'String');
% if ~isvarname(handles.analyzer.DataVarName)
%     error('invalid data var name');
% end

% throws if var not in base WS
function data = zlclGetDataFromBaseWS(handles) %#ok<STOUT>
cmd = sprintf('assignin(''caller'',''data'',%s);',get(handles.edDataVar,'String'));
evalin('base',cmd);
assert(exist('data','var')==1);
        
function varargout = OlyDatViewer_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function pbRunQry_Callback(hObject, eventdata, handles)
qObjs = handles.qcp.queryObjs;
experimentalQry = SAGE.Query.All(qObjs{:});
info = handles.analyzer.pullData(experimentalQry);
msg = sprintf('Data for %d experiments pulled.',info.numExperiments);
msgbox(msg,'Success');

function edDataVar_Callback(hObject, eventdata, handles)
newval = get(hObject,'String');
if ~isvarname(newval)
    set(hObject,'String',handles.analyzer.DataVarName);
else
    handles.analyzer.DataVarName = get(hObject,'String');
end

function pbRunAnalysis_Callback(hObject, eventdata, handles)
% get data from base WS
try
    data = zlclGetDataFromBaseWS(handles);
catch ME %#ok<NASGU>
    msgbox(sprintf('Cannot read variable ''%s'' from base workspace.',...
        get(handles.edDataVar,'String')),'Data not found');
    return;
end

% output dir for publishing
if get(handles.cbPublish,'Value')
    pubDir = get(handles.edOutputDir,'String');
    handles.analyzer.OutputDirectory = pubDir;
end    

% run
name = zlclCurrentAnalysisName(handles);
try
    [anlsStk failures] = handles.analyzer.runAnalysis(name,data);
catch ME
    msgbox(sprintf('Error: %s',ME.message));
    return;
end

% failures
Nanls = numel(anlsStk);
assert(Nanls==numel(failures));
strs = cell(Nanls,1);
for c = 1:Nanls
    strs{c} = sprintf('Anls %s: %d failures.',anlsStk{c}.Name,numel(failures{c}));
end
msgbox(char(strs),'Analysis Complete');

function cbPublish_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
tf = get(hObject,'Value');
if tf
    set(handles.edOutputDir,'Enable','inactive');
    handles.analyzer.PublishToHTML = true;
else
    set(handles.edOutputDir,'Enable','off');
    handles.analyzer.PublishToHTML = false;
end

function edOutputDir_Callback(hObject, eventdata, handles)
set(handles.edOutputDir,'Enable','inactive');

function edOutputDir_ButtonDownFcn(hObject, eventdata, handles)
dir = uigetdir('','Select output directory');
if ~isequal(dir,0)
    set(handles.edOutputDir,'String',dir);
end
set(handles.edOutputDir,'Enable','on');

function lbAnalyses_Callback(hObject, eventdata, handles)
zlclRefreshAnlsDesc(handles);

function edAnalysisDesc_Callback(hObject, eventdata, handles)

function edFailureDesc_Callback(hObject, eventdata, handles)

function pbPrevFailure_Callback(hObject, eventdata, handles)

function pbNextFailure_Callback(hObject, eventdata, handles)
disp('asd');

% function pumQryField_Callback(hObject, eventdata, handles)
% val = get(hObject,'Value');
% strs = get(hObject,'String');
% val = strs{val};
% fObj = handles.fieldName2Obj.(val);
% vv = fObj.validValues;
% if ~isempty(vv)
%     zlclQCSetCompareMode(get(hObject,'Parent'),'enumerated',vv);
% else
%     zlclQCSetCompareMode(get(hObject,'Parent'),'unenumerated');
% end

function zlclSetCurrentQry(obj,handles,qryName)
handles.currentQry = qryName;
set(handles.pnlQry,'Title',sprintf('Query: %s',qryName));
guidata(obj,handles);

function qryName = zlclGetCurrentQry(handles)
qryName = handles.currentQry;

function pbSaveQry_Callback(hObject, eventdata, handles)
[fname pname] = uiputfile('','Save to file',zlclGetCurrentQry(handles));
if isequal(fname,0) || isequal(pname,0)
    return;
end
fullfname = fullfile(pname,fname);
handles.qcp.fileWrite(fullfname);
zlclSetCurrentQry(hObject,handles,fname);

function pbLoadQry_Callback(hObject, eventdata, handles)
[fname pname] = uigetfile('*.qry','Select query file',zlclGetCurrentQry(handles));
if isequal(fname,0) || isequal(pname,0)
    return;
end
fullfname = fullfile(pname,fname);
handles.qcp.fileRead(fullfname);
zlclSetCurrentQry(hObject,handles,fname);
