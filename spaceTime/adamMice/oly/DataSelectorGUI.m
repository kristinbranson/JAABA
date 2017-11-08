function varargout = DataSelectorGUI(varargin)
%DataSelectorGUI MATLAB code for DataSelectorGUI.fig
%
%   See also OlyDat.DataSelector.

% Edit the above text to modify the response to help DataSelectorGUI

% Last Modified by GUIDE v2.5 31-May-2013 23:29:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DataSelectorGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @DataSelectorGUI_OutputFcn, ...
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

% dsGui = DataSelectorGUI(DataSelectorObj)
% DataSelectorObj needs to have its fDataPuller field set.
function DataSelectorGUI_OpeningFcn(hObject, eventdata, handles, varargin)
handles.output = hObject;
handles.qb = varargin{1};

pbWidth = 55;
dx = 5;
x = 10;
y = 8;
handles.pbEdit = uicontrol('tag','pbEdit','callback',@pbEditCbk,'keyPressFcn',@zlclPbKeyPressFcn,...
    'Units','pixels','position', [x y pbWidth 20],'string','Replace',...
    'FontName','default','FontUnits','pixels','FontSize',13,'parent',handles.pnlOuter);
handles.pbAppend = uicontrol('tag','pbAppend','callback',@pbAppendCbk,'keyPressFcn',@zlclPbKeyPressFcn,...
    'Units','pixels','position',[x+pbWidth+dx y pbWidth 20],'string','Add',...
    'FontName','default','FontUnits','pixels','FontSize',13,'parent',handles.pnlOuter);
handles.pbInsert = uicontrol('tag','pbInsert','callback',@pbInsertCbk,'keyPressFcn',@zlclPbKeyPressFcn,...
    'Units','pixels','position',[x+2*(pbWidth+dx) y pbWidth 20],'string','Insert',...
    'FontName','default','FontUnits','pixels','FontSize',13,'parent',handles.pnlOuter);
handles.pbDelete = uicontrol('tag','pbDelete','callback',@pbDeleteCbk,'keyPressFcn',@zlclPbKeyPressFcn,...
    'Units','pixels','position',[x+3*(pbWidth+dx) y pbWidth 20],'string','Delete',...
    'FontName','default','FontUnits','pixels','FontSize',13,'parent',handles.pnlOuter);
handles.pbDeleteAll = uicontrol('tag','pbDeleteAll','callback',@pbDeleteAllCbk,'keyPressFcn',@zlclPbKeyPressFcn,...
    'Units','pixels','position',[x+4*(pbWidth+dx) y pbWidth 20],'string','Reset',...
    'FontName','default','FontUnits','pixels','FontSize',13,'parent',handles.pnlOuter);

guidata(hObject, handles);
set(hObject,'Name',sprintf('OlyDat DataSelector: %s',handles.qb.AssayName));

if handles.qb.fDataPullerHasOptions
    set(handles.pumPullerOptions,'String',handles.qb.fDataPuller.options,'Visible','on');
else
    set(handles.pumPullerOptions,'Visible','off');
end
if handles.qb.fDataPuller.tfPullFSDefault
    set(handles.cbPullFileSys,'Value',1);    
else
    set(handles.cbPullFileSys,'Value',0);
end

% tab order
uistack(handles.pnlQuery,'top');
uistack(handles.pumField,'top');
uistack(handles.pumCompare,'top');
uistack(handles.etValue,'top');
uistack(handles.lbValidValues,'top');
uistack(handles.pbEdit,'top');
uistack(handles.pbAppend,'top');
uistack(handles.pbInsert,'top');
uistack(handles.pbDelete,'top');
uistack(handles.pbDeleteAll,'top');
uistack(handles.pbLoad,'top');
uistack(handles.pbSave,'top');

zlclUpdateFileSysControls(handles);

centerfig(hObject);

function varargout = DataSelectorGUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function pbSave_Callback(hObject,~,handles) %#ok<*INUSL,*DEFNU>
[fname pname] = uiputfile('*.qry','Save to file',handles.qb.getCurrentQuerySet);
if isequal(fname,0) || isequal(pname,0)
    return;
end
fullfname = fullfile(pname,fname);
handles.qb.setCurrentQuerySet(fullfname);
handles.qb.fileWrite();

function pbLoad_Callback(hObject,~,handles)
[fname pname] = uigetfile('*.qry','Select query file',handles.qb.getCurrentQuerySet);
if isequal(fname,0) || isequal(pname,0)
    return;
end
fullfname = fullfile(pname,fname);
handles.qb.fileRead(fullfname);
handles.qb.setCurrentQuerySet(fullfname);

function pumField_Callback(hObject, eventdata, handles)
handles.qb.newFieldSelection;

function pumCompare_Callback(hObject, eventdata, handles)
handles.qb.newCompareSelection;

function pbEditCbk(src,~,handles) %#ok<*INUSD>
if nargin < 3
    handles = guidata(src);
end
handles.qb.replaceQuery;

function pbAppendCbk(src,~,handles)
if nargin < 3
    handles = guidata(src);
end
handles.qb.appendQuery;

function pbInsertCbk(src,~,handles)
if nargin < 3
    handles = guidata(src);
end
handles.qb.insertQuery;

function pbDeleteCbk(src,~,handles)
if nargin < 3
    handles = guidata(src);
end
handles.qb.deleteQuery;

function pbDeleteAllCbk(src,~,handles)
if nargin < 3
    handles = guidata(src);
end
handles.qb.deleteAllQueries;
set(handles.pumField,'Value',1);
set(handles.pumCompare,'Value',1);
set(handles.etValue,'String','');
set(handles.lbValidValues,'Value',1);
handles.qb.newFieldSelection();

function pbRunQuery_Callback(hObject,~,handles)
handles.qb.runQuery();

function etValue_KeyPressFcn(hObject, eventdata, handles)
switch eventdata.Key
    case 'return'
        drawnow;
        handles.qb.appendQuery;
end

function lbValidValues_Callback(hObject, eventdata, handles)
switch get(gcf,'SelectionType')
    case 'open' % captures double-click and return
        handles.qb.appendQuery;
end

% Generic key-press function for pushbuttons that translates 'return' into
% 'press'
function zlclPbKeyPressFcn(hObject,eventdata)
handles = guidata(hObject);
c = get(handles.figure1,'CurrentCharacter');
if strcmp(c,sprintf('\r'))
    cbk = get(hObject,'Callback');
    cbk(hObject,[]);
end

function pbRunQuery_KeyPressFcn(hObject, eventdata, handles)
zlclPbKeyPressFcn(hObject,eventdata);

function pbLoad_KeyPressFcn(hObject, eventdata, handles)
zlclPbKeyPressFcn(hObject,eventdata);

function pbSave_KeyPressFcn(hObject, eventdata, handles)
zlclPbKeyPressFcn(hObject,eventdata);

function zlclUpdateFileSysControls(handles)
tf = logical(get(handles.cbPullFileSys,'Value'));
if tf
    enab = 'on';
else
    enab = 'off';
end
set(handles.etDataDir,'enable',enab);
set(handles.pbBrowseFS,'enable',enab);

function cbPullFileSys_Callback(hObject, eventdata, handles)
zlclUpdateFileSysControls(handles);

function pbBrowseFS_Callback(hObject, eventdata, handles)
d = uigetdir([],'Select data directory');
if isequal(d,0)
    % user canceled
elseif exist(d,'dir')~=7
    warning('DataSelectorGUI:dirNotFound','Directory ''%s'' not found.',d);
else
    set(handles.etDataDir,'String',d);
end
