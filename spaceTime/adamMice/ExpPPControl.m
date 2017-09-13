function varargout = ExpPPControl(varargin)
% EXPPPCONTROL MATLAB code for ExpPPControl.fig
%      EXPPPCONTROL, by itself, creates a new EXPPPCONTROL or raises the existing
%      singleton*.
%
%      H = EXPPPCONTROL returns the handle to a new EXPPPCONTROL or the handle to
%      the existing singleton*.
%
%      EXPPPCONTROL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXPPPCONTROL.M with the given input arguments.
%
%      EXPPPCONTROL('Property','Value',...) creates a new EXPPPCONTROL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ExpPPControl_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ExpPPControl_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ExpPPControl

% Last Modified by GUIDE v2.5 23-Apr-2015 11:01:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ExpPPControl_OpeningFcn, ...
                   'gui_OutputFcn',  @ExpPPControl_OutputFcn, ...
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

function ExpPPControl_OpeningFcn(hObject, eventdata, handles, varargin)
% ExpPPControl(data,callbacks,computeparams,hFig,varargin)
% data: data struct
% callbacks: callbacks
% computeparams: cell of PV args for ExpPP.computeStats
% hFig: optionally, handle to accompanying figure
% optional PVs
% - labelgram, scalar logical
    
data = varargin{1};
callbacks = varargin{2};
computeparams = varargin{3};
hFig = varargin{4};
assert(isstruct(data));
assert(iscell(computeparams));
assert(ishandle(hFig));  

varargin = varargin(5:end);
handles.tfLabelgram = myparse(varargin,'labelgram',false);
set(handles.txFlagLabelgram,'Visible',onIff(handles.tfLabelgram));

data = ExpPP.addNotGroups(data);

handles.data = data;
handles.Nexp = numel(data);
handles.callbacks = callbacks;
handles.computeparams = computeparams;
handles.hFig = hFig;

tbl = ExpPP.PARAMETER_DESCRIPTIONS;
data = [tbl.name num2cell(tbl.default) tbl.desc];
set(handles.tblParam,'Data',data);
handles.table_params = tbl;

handles.pbRecomputeColorActive = [0.85098 0.32549 0.098039];

handles.output = hObject;
guidata(hObject, handles);

mlv = version('-release');
switch mlv
  case '2014b'
    % none
  otherwise
    set(handles.pbExportFig,'Enable','off');
end
initCustomGroups(handles);

% By default, update is called initially
params = ExpPP.loadConfigVal('defaultParameters');
if ~isempty(params)
  loadParameters(handles,params);
  fprintf('Default parameters loaded.\n');
end
pbRecompute_Callback(handles.pbRecompute,[],handles);

set(handles.txExpParamsFile,'Enable','off');
set(handles.pbExpParams,'Enable','off');

% UIWAIT makes ExpPPControl wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function varargout = ExpPPControl_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function loadParameters(handles,parameters)
% parameters: Nx2 cell. Col1: parameter name. Col2: value.

currprms = handles.table_params.name;
data = get(handles.tblParam,'data');
assert(isequal(data(:,1),currprms));

[tf,loc] = ismember(currprms,parameters(:,1));
if any(~tf)
  prmsNoDefaults = currprms(~tf);
  str = civilizedStringFromCellArrayOfStrings(prmsNoDefaults);
  warning('ExpPPControl:noDefaults','No defaults saved for: %s',str);
end
%rows = find(tf);
loc = loc(tf);

data(tf,2) = parameters(loc,2);
set(handles.tblParam,'data',data);
%highlightParameterName(handles,rows);
if any(tf)
  recomputeNecessary(handles);
end

function params = getParamStruct(handles)
% get parameters/values from table

data = get(handles.tblParam,'Data');
names = data(:,1);
values = data(:,2);
assert(isequal(names,handles.table_params.name));
params = ExpPP.formParamStruct(names,values);

function recomputeAndRedraw(handles)
% use parameter values to:
% * update handles.data using ExpPP.computeStats
% * call handles.callback.drawAutoMarks

params = getParamStruct(handles);
params = struct2paramscell(params);

val = get(handles.cbExpParams,'Value');
if val
  expparamfile = get(handles.txExpParamsFile,'String');
  if isempty(expparamfile)
    error('ExpPPControl:fileNotFound','Please specify an experiment-specific parameter file.');
  end
  if exist(expparamfile,'file')==0
    error('ExpPPControl:fileNotFound','Experiment-specific parameter file ''%s'' not found.',expparamfile);
  end
  expparammap = ExpPP.loadExperimentSpecificParams(expparamfile);
else
  expparammap = [];
end

% recompute data
fprintf(1,'Computing...\n');
[handles.data,handles.stats] = ExpPP.computeStats(handles.data,params{:},handles.computeparams{:},'expparammap',expparammap);
guidata(handles.figure1,handles);

handles.callbacks.drawAutoMarks(handles.data);
redrawAutoStats(handles);

function redrawAutoStats(handles)
cgStr = get(handles.lbCustomGroups,'String');
cgIdx = get(handles.lbCustomGroups,'Value');
if ~isempty(cgStr) && ~isempty(cgIdx)
  cgrps = cgStr(cgIdx);
  cgrpflds = cellfun(@(x)sprintf(ExpPP.CUSTOMGROUPNAME_PAT,x),cgrps,'uni',0);
else
  cgrpflds = {};
end
handles.callbacks.drawAutoMarkGroupStats(handles.stats,cgrpflds);

function recomputeNecessary(handles)
set(handles.pbRecompute,...
  'ForegroundColor',handles.pbRecomputeColorActive,...
  'Fontweight','bold');
drawnow;

% function highlightParameterName(handles,rows)
% % This function highlights a parameter name in red. Using this is causing
% % laggy/weird Java behavior; the jtable.setValueAt call triggers the
% % tblParam_CellEdit callback, and callbacks seem to build up, have weird
% % race conditions etc. Short-circuiting recursion in
% % tblParam_CellEditCallback does not seem to help.
% jscroll = findjobj(handles.tblParam); 
% jtable = jscroll.getViewport.getComponent(0); 
% for r = rows(:)'
%   param = handles.table_params.name{r};
%   clrparam = sprintf('<html><font color="red">%s</font></html>',param);
%   jtable.setValueAt(java.lang.String(clrparam),r-1,0); 
% end

function tblParam_CellEditCallback(hObject, eventdata, handles)
% row = eventdata.Indices(1);
% col = eventdata.Indices(2);
%highlightParameterName(handles,row);
recomputeNecessary(handles);

function pbExport_Callback(hObject, eventdata, handles)
uiwait(msgbox('You will be prompted to save two files: one for the raw data, one for statistics.'));
format = 'csv';
if handles.tfLabelgram
  ExpPP.export(handles.data,[],'format',format,'expppstatprefix','labl');
else
  ExpPP.export(handles.data,[],'format',format);
end
ExpPP.exportStats(handles.stats);

function pbExportWS_Callback(hObject, eventdata, handles)
fprintf(1,'Assigning to variables ''data'', ''stats'' in base workspace.\n');
data = ExpPP.cleanNotGroups(handles.data);
assignin('base','data',data);
assignin('base','stats',handles.stats);

function pbDone_Callback(hObject, eventdata, handles)
if ~isempty(handles.hFig) && ishandle(handles.hFig)
  delete(handles.hFig);
end
delete(handles.figure1);
  
function lbCustomGroups_Callback(hObject, eventdata, handles)
customGrps = get(handles.lbCustomGroups,'String');
selIdx = get(handles.lbCustomGroups,'Value');

% update selected exps
if isempty(selIdx)
  tf = false(handles.Nexp,1);
else
  % When more than one group is selected, set selected exps to first
  % selected group only
  selIdx = selIdx(1); 
  selGrp = customGrps{selIdx};
  assert(~isempty(handles.data),'Currently unsupported for empty data.');
  assert(ismember(selGrp,handles.data(1).(ExpPP.FIELDNAME_CUSTOMGROUPLIST)));
  tfFld = sprintf(ExpPP.CUSTOMGROUPNAME_PAT,selGrp);
  tf = [handles.data.(tfFld)]';
  assert(numel(tf)==handles.Nexp);
end
handles.callbacks.setSelectedExps(tf);

% update stats
redrawAutoStats(handles);

function initCustomGroups(handles)
assert(~isempty(handles.data),'Data cannot be empty.');
grps = handles.data(1).(ExpPP.FIELDNAME_CUSTOMGROUPLIST);
set(handles.lbCustomGroups,'String',grps);
for i = 1:numel(grps)
  grpname = grps{i};
  tf = ExpPP.findCustomGroup(handles.data,grpname);
  handles.callbacks.addCustomGroupViz(grpname,tf);
end
if ~isempty(grps)
  set(handles.lbCustomGroups,'Value',1);
  tf = ExpPP.findCustomGroup(handles.data,grps{1});
  handles.callbacks.setSelectedExps(tf);
end

function pbCreateCustomGroup_Callback(hObject, eventdata, handles)
grpname = inputdlg('Name for custom group?','Create custom group');
if isempty(grpname)
  % none; user canceled
else
  grpname = grpname{1};
  tfSelectedExps = handles.callbacks.getSelectedExps();
  handles.data = ExpPP.createCustomGroup(handles.data,tfSelectedExps,grpname);
  handles.callbacks.addCustomGroupViz(grpname,tfSelectedExps);
  customGrps = get(handles.lbCustomGroups,'String');
  customGrps{end+1,1} = grpname;
  selIdx = numel(customGrps);
  set(handles.lbCustomGroups,'String',customGrps,'Value',selIdx);

  guidata(hObject,handles);
  recomputeNecessary(handles);
end

function pbRemoveCustomGroup_Callback(hObject, eventdata, handles)
customGrps = get(handles.lbCustomGroups,'String');
if isempty(customGrps)
  % no-op
  return;
end
idxs = get(handles.lbCustomGroups,'Value');
grpsToRm = customGrps(idxs);
for i = 1:numel(grpsToRm)
  handles.data = ExpPP.removeCustomGroup(handles.data,grpsToRm{i});
  handles.callbacks.rmCustomGroupViz(grpsToRm{i});
end
set(handles.lbCustomGroups,'Value',1);
set(handles.lbCustomGroups,'String',setdiff(customGrps,grpsToRm));
guidata(hObject,handles);
recomputeNecessary(handles);

function pbRecompute_Callback(hObject, eventdata, handles)
recomputeAndRedraw(handles);
fprintf(' Done\n');
% data = get(handles.tblParam,'data');
% data(:,1) = handles.table_params.name;
% set(handles.tblParam,'data',data);
set(hObject,'Fontweight','normal','ForegroundColor',[0 0 0]);

function menuFileSaveDefaults_Callback(hObject, eventdata, handles)
% AL: This parameter stuff (config file, defaults, struct->flattened format 
% etc) could use some thought
data = get(handles.tblParam,'Data');
params = [handles.table_params.name,data(:,2)]; 
ExpPP.saveConfigVal('defaultParameters',params);
fprintf('Default parameters saved.\n');

function menuFileLoadDefaults_Callback(hObject, eventdata, handles)
params = ExpPP.loadConfigVal('defaultParameters');
if isempty(params)
  warningNoTrace('ExpPPControl:noDefaults','No default parameters saved.');
else
  loadParameters(handles,params);
end

function pbExportFig_Callback(hObject, eventdata, handles)
handles.callbacks.exportFig();

function cbExpParams_Callback(hObject, eventdata, handles)
val = get(hObject,'Value');
if val
  set(handles.txExpParamsFile,'Enable','on');
  set(handles.pbExpParams,'Enable','on');
else
  set(handles.txExpParamsFile,'Enable','off');
  set(handles.pbExpParams,'Enable','off');
end
recomputeNecessary(handles);

function pbExpParams_Callback(hObject, eventdata, handles)
expparampath = ExpPP.loadConfigVal('expparampath');
if isempty(expparampath)
  expparampath = pwd;
end
[fname,pname] = uigetfile({'*.csv'},'Select parameter file',expparampath);
if isequal(fname,0)
  return;
end
fname = fullfile(pname,fname);
set(handles.txExpParamsFile,'String',fname);
recomputeNecessary(handles);

