function varargout = CurationInfoEntry(varargin)
% CURATIONINFOENTRY MATLAB code for CurationInfoEntry.fig
%      CURATIONINFOENTRY by itself, creates a new CURATIONINFOENTRY or raises the
%      existing singleton*.
%
%      H = CURATIONINFOENTRY returns the handle to a new CURATIONINFOENTRY or the handle to
%      the existing singleton*.
%
%      CURATIONINFOENTRY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CURATIONINFOENTRY.M with the given input arguments.
%
%      CURATIONINFOENTRY('Property','Value',...) creates a new CURATIONINFOENTRY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CurationInfoEntry_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CurationInfoEntry_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CurationInfoEntry

% Last Modified by GUIDE v2.5 03-Jun-2011 12:30:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CurationInfoEntry_OpeningFcn, ...
                   'gui_OutputFcn',  @CurationInfoEntry_OutputFcn, ...
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

% infoStruct: curation info struct. Can be nonscalar for batch edit.
% username: Browser user (manual curator)
% cbkFcn: function handle to be called on a successful curation edit. Sig
% is cbkFcn(s), where s is a scalar struct with the curation fields.
function CurationInfoEntry_OpeningFcn(hObject, eventdata, handles, varargin)

% TODO: Note similarity between this and AggCurationInfoEntry_OpeningFcn

info = varargin{1};
username = varargin{2};
handles.cbkFcn = varargin{3};

assert(numel(info)>=1);

% Based on how this dialog works when closed, info should correspond
% exactly to the currently selected experiments in browser.

% Pre-populate fields as appropriate
tfBatchEntry = numel(info)>1;
if tfBatchEntry
    set(handles.pnlCurationInfo,'Title','Multiple experiments selected');
    zlclSetPUMToValueWithDefault(handles.pumManualPf,[]);
    zlclSetPUMToValueWithDefault(handles.pumFlagRedo,[]);
    zlclSetPUMToValueWithDefault(handles.pumFlagReview,[]);
    set(handles.etNotes,'String',[]);   
else
    set(handles.pnlCurationInfo,'Title',sprintf('Experiment %d:  %s',info.experiment_id,info.experiment_name));
    zlclSetPUMToValueWithDefault(handles.pumManualPf,info.manual_pf);
    zlclSetPUMToValueWithDefault(handles.pumFlagRedo,info.flag_redo);
    zlclSetPUMToValueWithDefault(handles.pumFlagReview,info.flag_review);
    set(handles.etNotes,'String',info.notes_curation);
end

% Auto-populate curator, curationDate
pnlBgColor = get(handles.pnlCurationInfo,'BackgroundColor');
set(handles.etManualCurator,'String',username,'BackGroundColor',pnlBgColor);
set(handles.etCurationDate,'String',datestr(now,30),'BackGroundColor',pnlBgColor);

% Note: Logic for "existing marks" warning differs from AggCurationInfoEntry
allCurators = {info.manual_curator}';
allCurationDates = {info.manual_curation_date}';
tfNoExistingCurators = all(cellfun(@isempty,allCurators));
tfNoExistingCurationDates = all(cellfun(@isempty,allCurationDates));
tfNoExistingMarks = tfNoExistingCurators && tfNoExistingCurationDates;

if tfNoExistingMarks
    set(handles.txWarningExistingMarks,'Visible','off');
else
    set(handles.txWarningExistingMarks,'Visible','on');
end
guidata(hObject,handles);

tweakGUILayout(handles);

% Determine the position of the dialog - centered on the callback figure
% if available, else, centered on the screen
FigPos=get(0,'DefaultFigurePosition');
OldUnits = get(hObject, 'Units');
set(hObject, 'Units', 'pixels');
OldPos = get(hObject,'Position');
FigWidth = OldPos(3);
FigHeight = OldPos(4);
if isempty(gcbf)
    ScreenUnits=get(0,'Units');
    set(0,'Units','pixels');
    ScreenSize=get(0,'ScreenSize');
    set(0,'Units',ScreenUnits);

    FigPos(1)=1/2*(ScreenSize(3)-FigWidth);
    FigPos(2)=2/3*(ScreenSize(4)-FigHeight);
else
    GCBFOldUnits = get(gcbf,'Units');
    set(gcbf,'Units','pixels');
    GCBFPos = get(gcbf,'Position');
    set(gcbf,'Units',GCBFOldUnits);
    FigPos(1:2) = [(GCBFPos(1) + GCBFPos(3) / 2) - FigWidth / 2, ...
                   (GCBFPos(2) + GCBFPos(4) / 2) - FigHeight / 2];
end
FigPos(3:4)=[FigWidth FigHeight];
set(hObject, 'Position', FigPos);
set(hObject, 'Units', OldUnits);

% modal
set(handles.figure1,'WindowStyle','modal')
uiwait(handles.figure1);

function varargout = CurationInfoEntry_OutputFcn(hObject, eventdata, handles)
varargout{1} = [];
delete(handles.figure1);

function pbOk_Callback(hObject, eventdata, handles) %#ok<*INUSL>
s = struct();
s.manual_pf = zlclGetPUMValueString(handles.pumManualPf);
s.flag_redo = zlclGetPUMValueString(handles.pumFlagRedo);
s.flag_review = zlclGetPUMValueString(handles.pumFlagReview);
s.notes_curation = get(handles.etNotes,'String');
s.manual_curator = get(handles.etManualCurator,'String');
s.manual_curation_date = get(handles.etCurationDate,'String');
handles.cbkFcn(s);
uiresume(handles.figure1);

function pbCancel_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
uiresume(handles.figure1);

function figure1_CloseRequestFcn(hObject, eventdata, handles) %#ok<*INUSD>
if isequal(get(hObject, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, us UIRESUME
    uiresume(hObject);
else
    % The GUI is no longer waiting, just close it
    delete(hObject);
end

function figure1_KeyPressFcn(hObject, eventdata, handles)
if isequal(get(hObject,'CurrentKey'),'escape')
    uiresume(handles.figure1);
end        

function s = zlclGetPUMValueString(hObject)
str = get(hObject,'String');
val = get(hObject,'Value');
s = str{val};

% pumH: handle to popup-menu uicontrol
% val: str enum: {'', <opt1>, <opt2>, ...}
% If val is one of the nonempty options, set pumH to that value (pumH is
% expected to have that value as a valid option). If val is the empty str,
% set pumH to its first value.
function zlclSetPUMToValueWithDefault(pumH,val)
if isempty(val)
    set(pumH,'Value',1);    
else
    str = get(pumH,'String');
    tf = strcmp(val,str);
    assert(nnz(tf)==1);
    set(pumH,'Value',find(tf));
end

function pumManualPf_Callback(hObject, eventdata, handles)
%uicontrol(handles.pumFlagRedo);

function pumFlagRedo_Callback(hObject, eventdata, handles)
%uicontrol(handles.pumFlagReview);

function pumFlagReview_Callback(hObject, eventdata, handles)
%uicontrol(handles.etManualCurator);

function etNotes_Callback(hObject, eventdata, handles)
c = get(handles.figure1,'CurrentCharacter');
if strcmp(c,sprintf('\r'))
    uicontrol(handles.pbOk);
end

function pbOk_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'return'
        pbOk_Callback(hObject,[],handles);
    case 'escape'
        % this is not working
        set(handles.figure1,'CurrentObject',handles.figure1);
end

function pbCancel_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'return'
        pbCancel_Callback(hObject,[],handles);
    case 'escape'
        % this is not working
        set(handles.figure1,'CurrentObject',handles.figure1);
end

function pumManualPf_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'return'
        uicontrol(handles.pumFlagRedo);
end

function pumFlagRedo_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'return'
        uicontrol(handles.pumFlagReview);
end

function pumFlagReview_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'return'
        uicontrol(handles.etNotes);
end

function etCurationDate_Callback(hObject, eventdata, handles)
