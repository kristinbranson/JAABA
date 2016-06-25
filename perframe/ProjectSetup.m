function varargout = ProjectSetup(varargin)
% PROJECTSETUP MATLAB code for ProjectSetup.fig
%      PROJECTSETUP, by itself, creates a new PROJECTSETUP or raises the existing
%      singleton*.
%
%      H = PROJECTSETUP returns the handle to a new PROJECTSETUP or the handle to
%      the existing singleton*.
%
%      PROJECTSETUP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PROJECTSETUP.M with the given input arguments.
%
%      PROJECTSETUP('Property','Value',...) creates a new PROJECTSETUP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ProjectSetup_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ProjectSetup_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ProjectSetup

% Last Modified by GUIDE v2.5 22-Apr-2013 21:12:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ProjectSetup_OpeningFcn, ...
                   'gui_OutputFcn',  @ProjectSetup_OutputFcn, ...
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
% --- Executes just before ProjectSetup is made visible.
function ProjectSetup_OpeningFcn(hObject, ~, handles, varargin)
% ProjectSetup(p1,v1,...)
%
% ProjectSetup configures a JAABA project. It can optionally accept an
% existing project configuration (basicParamsStruct). It returns a
% configured project in the supplied 'handleobj'. A project configuration
% currently looks like a stripped-down Macguffin and consists of 
%   - Names of behaviors/classifiers
%   - Name of experimental movie files
%   - etc.
%
% Call ProjectSetup as follows:
% hobj = HandleObj;
% uiwait(ProjectSetup('handleobj',hobj,...<other args>...));
% processResultsOrDoWhatever(hobj.data); % hobj.data is a Macguffin
%
% P-V arguments:
% * basicParamsStruct: Either [], or Macguffin (or Macguffin-like struct).
% Supplies initial values for project parameters. If basicParamStruct==[],
% then we are creating a new project.
% * handleobj: A HandleObj. 

handles.output = hObject;

[figureJLabel, ...
 basicParamsStruct,...
 defaultmoviefilename,...
 defaultmovieindexfilename,...
 defaulttrxfilename,...
 handles.handleobj] = ...
   myparse(varargin,...
           'figureJLabel',[],...
           'basicParamsStruct',[],...
           'defaultmoviefilename',0,...
           'defaultmovieindexfilename',0,...
           'defaulttrxfilename',0,...
           'handleobj',[]);

handles.new = isempty(basicParamsStruct);
handles.figureJLabel = figureJLabel;
handles.basicParamsStruct = basicParamsStruct;

adjustColorsIfMac(hObject);

% Need to derive the basic and advanced sizes (part of the model) from
% the current figure dimensions
handles = setBasicAndAdvancedSizesToMatchFigure(handles);

handles.mode = 'basic';
handles = updateFigurePosition(handles);

centerOnParentFigure(hObject,figureJLabel);

% Initialize the list of possible feature lexicon names
[handles.featureLexiconNameList,xmlList] = getFeatureLexiconListsFromXML();
handles.featureLexiconNameListIsST = cellfun(@xmlParamsIsST,xmlList);
assert(numel(handles.featureLexiconNameList)==numel(handles.featureLexiconNameListIsST));

% Populate the featureLexiconName popuplist with options
set(handles.featureconfigpopup,'String',handles.featureLexiconNameList);

if isempty(basicParamsStruct)
  % New project
  warnst = warning('off','MATLAB:structOnObject');
  handles.basicParamsStruct = struct(Macguffin('flies')); % the default featureLexiconName
  warning(warnst);
  if ischar(defaultmoviefilename)
    handles.basicParamsStruct.file.moviefilename = defaultmoviefilename;
  end
  if ischar(defaultmovieindexfilename)
    handles.basicParamsStruct.file.movieindexfilename = defaultmovieindexfilename;
  end
  if ischar(defaulttrxfilename)
    handles.basicParamsStruct.file.trxfilename = defaulttrxfilename;
  end  
else
  warnst = warning('off','MATLAB:structOnObject');
  handles.basicParamsStruct = struct(basicParamsStruct);
  warning(warnst);
  
  % Quirk: basicParamsStruct.behaviors
  % For existing projects basicParamsStruct.behaviors.names includes the
  % no-behavior names (eg 'None', 'No_Chase', etc). In ProjectSetup, the
  % user currently only specifies the positive names.
  % So, within ProjectSetup, we truncate .behaviors.names to just the real
  % names; as a quirk, we leave .behaviors.labelcolors unmodified, so
  % .labelcolors DOES include no-behaviors, as the user may want to change
  % those colors.
  behnames = Labels.verifyBehaviorNames(handles.basicParamsStruct.behaviors.names);
  handles.basicParamsStruct.behaviors.names = behnames;
end

if ~isfield(handles.basicParamsStruct,'extra'),
  handles.basicParamsStruct.extra = struct;
end
if ~isfield(handles.basicParamsStruct.extra,'usePastOnly'),
  handles.basicParamsStruct.extra.usePastOnly = false;
end

% Update the configuration table in the GUI
updateConfigTable(handles);

set(hObject,'name',fif(handles.new,'New...','Edit...'));
onoff = fif(handles.new,'on','off');
set(handles.featureconfigpopup,'enable',onoff);
set(handles.pushbutton_copy,'visible',onoff);
set(handles.editName,'enable',onoff);
set(handles.edittrxfilename,'enable',onoff);
set(handles.pushbutton_perframe,'visible',onoff);

% Set the current score-as-feature file
fileNameList = {handles.basicParamsStruct.scoreFeatures(:).classifierfile};
handles.indexOfScoreFeaturesFile = fif(isempty(fileNameList),[],1);
  
guidata(hObject, handles);

% Update all the text fields in the figure
updateEditsListboxesAndPopupmenus(handles);

% Make the figure visible
set(hObject,'Visible','on');

return


% -------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function varargout = ProjectSetup_OutputFcn(hObject, ~, handles)  %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
%varargout{2} = handles.outfile;
%delete(handles.figureProjectSetup);
return


% -------------------------------------------------------------------------
function handles=setBasicAndAdvancedSizesToMatchFigure(handles)
% Derive the basic and advanced sizes (part of the model) from
% the current figure dimensions
curPos = get(handles.figureProjectSetup,'Position');
tablePos = get(handles.config_table,'Position');
reducedWidth = tablePos(1)-15;
reducedHeight = curPos(4);
handles.advancedSize = curPos(3:4);
handles.basicSize = [reducedWidth reducedHeight];
return


% -------------------------------------------------------------------------
function handles = updateFigurePosition(handles)
% Update the figure to match the current mode (which is in the model)
if strcmp(handles.mode,'basic')
  curPos = get(handles.figureProjectSetup,'Position');
  set(handles.figureProjectSetup,'Position',[curPos(1:2) handles.basicSize]);
  set(handles.togglebutton_advanced,'Value',0);
else
  curPos = get(handles.figureProjectSetup,'Position');
  set(handles.figureProjectSetup,'Position',[curPos(1:2) handles.advancedSize]);
  set(handles.togglebutton_advanced,'Value',1);  
end
return


% -------------------------------------------------------------------------
function updateEditsListboxesAndPopupmenus(handles)

% update the behavior name editbox
if isfield(handles.basicParamsStruct.behaviors,'names');
  names = handles.basicParamsStruct.behaviors.names;
  if numel(names)>0    
    namestr = String.cellstr2CommaSepList(names);
  else
    namestr = '';
  end
  set(handles.editName,'String',namestr);
else
  set(handles.editName,'String','');
end

fnames = {'scorefilename','moviefilename','trxfilename'};
boxnames = {'editscorefilename','editmoviefilename','edittrxfilename'};
for ndx = 1:numel(fnames)
  curf = fnames{ndx};
  curbox = boxnames{ndx};
  if isfield(handles.basicParamsStruct.file,curf);
    val = handles.basicParamsStruct.file.(curf);
    if iscellstr(val)
      val = String.cellstr2CommaSepList(val);
    end
    set(handles.(curbox),'String',val);
  else
    set(handles.(curbox),'String','');
  end
end

% Update the select feature dictionary name
indexOfFeatureLexicon = find(strcmp(handles.basicParamsStruct.featureLexiconName,handles.featureLexiconNameList));
set(handles.featureconfigpopup,'Value',indexOfFeatureLexicon);
tfST = handles.featureLexiconNameListIsST(indexOfFeatureLexicon);
if tfST
  visVal = 'on';
else
  visVal = 'off';
end
set(handles.txST,'Visible',visVal);

% % Update the list of scores-an-inputs
% fileNameList = {handles.basicParamsStruct.scoreFeatures(:).classifierfile};
% set(handles.listbox_inputscores,'String',fileNameList);
% set(handles.listbox_inputscores,'Value',handles.indexOfScoreFeaturesFile);

% % Disble the Remove button iff the list of scores-as-inputs is empty
% set(handles.pushbutton_removelist,'enable',offIff(isempty(fileNameList)));

return


% -------------------------------------------------------------------------
function handles = fileNameEditTwiddled(handles,editName)
% Deal with one of the file name edits being modified.
% This is like a controller method.
varName=editName(5:end);  % trim 'edit' off of start
newFileName=strtrim(get(handles.(editName),'String'));
allowEmpty = ismember(varName,{'moviefilename'});
allowMulti = ismember(varName,{'scorefilename'});
[handles,whatHappened]=setFileName(handles,varName,newFileName,allowEmpty,allowMulti);
warnIfInvalid(varName,whatHappened);
updateFileNameEdit(handles,editName);
updateConfigTable(handles);
return


% -------------------------------------------------------------------------
function [handles,whatHappened] = setFileName(handles,varName,newFileName,allowEmpty,allowMulti)
% Set the given file name in the model, if it's a valid name.  whatHappened
% is a string that can be 'notChanged', 'emptyEntry','changed', or
% 'invalidEntry', depending.

fileName=handles.basicParamsStruct.file.(varName);

if allowEmpty && isempty(newFileName),
  handles.basicParamsStruct.file.(varName)=newFileName;
  if strcmp(newFileName,fileName),
    whatHappened='notChanged';
  else
    whatHappened='changed';
  end
  return;
end

if isequal(newFileName,fileName)
  whatHappened='notChanged';
elseif isempty(newFileName)
  whatHappened='emptyEntry';
elseif ~allowMulti && IsNiceFileName(newFileName)
  whatHappened='changed';
  handles.basicParamsStruct.file.(varName)=newFileName;
elseif allowMulti
  newFileNameCellStr = String.commaSepList2CellStr(newFileName);
  if all(cellfun(@IsNiceFileName,newFileNameCellStr))
    whatHappened = 'changed';
    handles.basicParamsStruct.file.(varName)=newFileNameCellStr;
  else
    whatHappened = 'invalidEntry';
  end
else
  whatHappened='invalidEntry';
end

return


% -------------------------------------------------------------------------
function warnIfInvalid(varName,whatHappened)
% Throw up a warning dialog if whatHappened equals 'invalidEntry'
if isequal(whatHappened,'invalidEntry')
  message=...
    sprintf(['The name specified for %s cannot have special characters.'...
             'Please use only alphanumeric characters and underscore.'],varName);
  title='Invalid File Name';
  uiwait(warndlg(message,title,'modal'));
end
return


% -------------------------------------------------------------------------
function updateFileNameEdit(handles,editName)
  % Update the named file name edit to match the model.
  varName=editName(5:end);  % trim 'edit' off of start
  val = handles.basicParamsStruct.file.(varName);
  if iscellstr(val)
    val = String.cellstr2CommaSepList(val);
  end
  set(handles.(editName),'String',val);
return


% -------------------------------------------------------------------------
function updateConfigTable(handles)
% Update the config table (a GUI element) to match the current "model"
% state
basicParamsStruct = handles.basicParamsStruct;
fields2remove = {'featureLexicon','windowFeaturesParams','scoreFeatures', ...
                 'sublexiconPFNames','labels','gtLabels','expDirNames', ...
                 'gtExpDirNames', 'classifierStuff', 'version','classifierfile',...
                 'windowFeaturesParams'};
if ~handles.new
  editModeRemove = {'featureLexiconName' 'behaviors.names' 'behaviors.type' ...
    'behaviors.nbeh' 'file.trxfilename' 'extra.perframe'};
  fields2remove = [fields2remove editModeRemove];
end
basicParamsStruct = structrmfield(basicParamsStruct,fields2remove);
data = GetParamsAsTable(basicParamsStruct);
set(handles.config_table,'Data',data);
return


% -------------------------------------------------------------------------
function data = GetParamsAsTable(basicParamsStruct)
data = addToList(basicParamsStruct,{},'');
idx = cellfun(@iscell,data(:,2));
if any(idx),
  for i = find(idx(:)'),
    if all(cellfun(@ischar,data{i,2})),
      data{i,2} = sprintf('%s,',data{i,2}{:});
      if numel(data{i,2})>0 && data{i,2}(end) == ',',
        data{i,2} = data{i,2}(1:end-1);
      end
    end
  end
end

if any(cellfun(@iscell,data(:,2))),
  data = {}; 
  return;
end


% -------------------------------------------------------------------------
function list = addToList(curStruct,list,pathTillNow)
if isempty(fieldnames(curStruct)), return; end
fnames = fieldnames(curStruct);
for ndx = 1:numel(fnames)
  if isstruct(curStruct.(fnames{ndx})),
    list = addToList(curStruct.(fnames{ndx}),list,[pathTillNow fnames{ndx} '.']);
  else
    list{end+1,1} = [pathTillNow fnames{ndx}]; %#ok<AGROW>
    param = curStruct.(fnames{ndx});
    if isnumeric(param)
      if isempty(param)
        q='';
      else
        q = num2str(param(1));
        for jj = 2:numel(param)
          q = [q ',' num2str(param(jj))]; %#ok<AGROW>
        end
      end
      list{end,2} = q;
    else
      list{end,2} = param;
    end
  end
end

% -------------------------------------------------------------------------
function editName_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editName as text
%        str2double(get(hObject,'String')) returns contents of editName as a double
name = String.commaSepList2CellStr(get(hObject,'String'));
if any(cellfun(@(x)isempty(regexp(x,'^[a-zA-Z][\w_,]*$','once','start')),name));
   uiwait(warndlg(['The behavior name cannot have special characters.'...
                   'Please use only alphanumeric characters and _']));
   updateEditsListboxesAndPopupmenus(handles);              
   return
end
    
handles.basicParamsStruct.behaviors.names = name;
if numel(name)>1
  % multiple classifiers; set labelcolors
  nlabels = 2*numel(name);
  clrs = Labels.cropOrAugmentLabelColors(zeros(1,0),nlabels,'lines');
  handles.basicParamsStruct.behaviors.labelcolors = clrs;
end
handles.basicParamsStruct.file.scorefilename = cellfun(@(x)sprintf('scores_%s.mat',x),name,'uni',0);
guidata(hObject,handles);
updateEditsListboxesAndPopupmenus(handles); 
updateConfigTable(handles);
return

function behaviorsStruct = enforceBehaviorParamConstraints(behaviorsStruct,dowarn)

if ~exist('dowarn','var')
  dowarn = false;
end

% Add no-behavior names
behaviorsStruct.names = Labels.behnames2labelnames(behaviorsStruct.names);

nlabels = numel(behaviorsStruct.names);
if numel(behaviorsStruct.labelcolors)~=3*nlabels
  if dowarn
    warning('ProjectSetup:behaviorParams',...
      'Specified label colors not consistent with number of behaviors. Updating.');
  end
  behaviorsStruct.labelcolors = Labels.cropOrAugmentLabelColors(behaviorsStruct.labelcolors,nlabels,'lines');
end


% -------------------------------------------------------------------------
function editName_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return


% -------------------------------------------------------------------------
% --- Executes on selection change in featureconfigpopup.
function featureconfigpopup_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to featureconfigpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
iFeatureLexicon=get(hObject,'Value');
featureLexiconNameNew = handles.featureLexiconNameList{iFeatureLexicon};
% Save the current values of the main four fields
% (Note that any changes in the advanced part of the window will be lost.)
behaviorName=handles.basicParamsStruct.behaviors.names;
movieFileName=handles.basicParamsStruct.file.moviefilename;
movieIndexFileName=handles.basicParamsStruct.file.movieindexfilename;
trackFileName=handles.basicParamsStruct.file.trxfilename;
scoreFileName=handles.basicParamsStruct.file.scorefilename;
% Replace basicParamsStruct with one appropriate to the new feature lexicon
% name
old=warning('query','MATLAB:structOnObject');
warning('off','MATLAB:structOnObject');  % turn off annoying warning
handles.basicParamsStruct = struct(Macguffin(featureLexiconNameNew));
warning(old);  % restore annoying warning  
% Now restore fields from the old basicParamsStruct
handles.basicParamsStruct.behaviors.names=behaviorName;
handles.basicParamsStruct.file.moviefilename=movieFileName;
handles.basicParamsStruct.file.movieindexfilename=movieIndexFileName;
handles.basicParamsStruct.file.trxfilename=trackFileName;
handles.basicParamsStruct.file.scorefilename=scoreFileName;
guidata(hObject,handles);
updateConfigTable(handles);
updateEditsListboxesAndPopupmenus(handles);
return


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function featureconfigpopup_CreateFcn(hObject, ~, ~)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
function editscorefilename_Callback(hObject, eventdata, handles)  %#ok
editName=get(hObject,'tag');
handles = fileNameEditTwiddled(handles,editName);
guidata(hObject,handles);


% -------------------------------------------------------------------------
function editscorefilename_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
function editmoviefilename_Callback(hObject, eventdata, handles)  %#ok
editName=get(hObject,'tag');
handles = fileNameEditTwiddled(handles,editName);
guidata(hObject,handles);


% -------------------------------------------------------------------------
function editmoviefilename_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
function edittrxfilename_Callback(hObject, eventdata, handles)  %#ok
editName=get(hObject,'tag');
handles = fileNameEditTwiddled(handles,editName);
guidata(hObject,handles);


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function edittrxfilename_CreateFcn(hObject, eventdata, handles)  %#ok
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% -------------------------------------------------------------------------
function pushbutton_cancel_Callback(hObject, eventdata, handles)  %#ok
handles.outfile = 0;
guidata(hObject,handles);
delete(findAncestorFigure(hObject));              
%uiresume(handles.figureProjectSetup);
return


% -------------------------------------------------------------------------
function pushbutton_done_Callback(hObject, ~, handles)  %#ok
behaviorName = handles.basicParamsStruct.behaviors.names;
if isempty(behaviorName)
  uiwait(errordlg('You must enter a behavior name.', ...
                  'No behavior name', ...
                  'modal'));
  return
end

% Check for an empty track file name
trackFileName = handles.basicParamsStruct.file.trxfilename;
if isempty(trackFileName)
  uiwait(errordlg('You must enter a track file name.', ...
                  'No track file name', ...
                  'modal'));
  return
end

% Check for an empty score file name
scoreFileName = handles.basicParamsStruct.file.scorefilename;
if isempty(scoreFileName)
  uiwait(errordlg('You must enter a score file name.', ...
                  'No score file name', ...
                  'modal'));
  return
end

nbehavior = numel(handles.basicParamsStruct.behaviors.names);
assert(iscellstr(handles.basicParamsStruct.file.scorefilename));
nscorefile = numel(handles.basicParamsStruct.file.scorefilename);
if nbehavior~=nscorefile
  uiwait(errordlg('The number of score files must match the number of behaviors.', ...
                  'Invalid score file names', ...
                  'modal'));
  return
end

% if movie is an mjpg, ask for an index file
[~,~,ext] = fileparts(handles.basicParamsStruct.file.moviefilename);
if ismember(lower(ext),{'.mjpg','.mjpeg'}),
  
  res = questdlg('Is this an indexed MJPG file? If so, you will be prompted to select the corresponding index file location.');
  if strcmpi(res,'Yes'),    
    if isfield(handles.basicParamsStruct.file,'movieindexfilename') && ...
        ischar(handles.basicParamsStruct.file.movieindexfilename),
      defaultindexfile = handles.basicParamsStruct.file.movieindexfilename;
    else
      defaultindexfile = 'index.txt';
    end
    res = inputdlg({'MJPG index file name: '},'MJPG index file name',1,{defaultindexfile});
    % cancel
    if isempty(res),
      return;
    end
    handles.basicParamsStruct.file.movieindexfilename = res{1};
  else
    handles.basicParamsStruct.file.movieindexfilename = 0;
  end
  
else
  handles.basicParamsStruct.file.movieindexfilename = 0;
end

% Finish initialization of basicParamsStruct and return as Macguffin in handleobj
handles.basicParamsStruct.behaviors = ...
  enforceBehaviorParamConstraints(handles.basicParamsStruct.behaviors,true);

if handles.new  
  for i = nbehavior:-1:1
    cs(i,1) = ClassifierStuff();
    wfp{i,1} = handles.basicParamsStruct.windowFeaturesParams{1};
  end
  handles.basicParamsStruct.classifierStuff = cs;
  handles.basicParamsStruct.windowFeaturesParams = wfp;
end

basicParamsStruct = handles.basicParamsStruct;
% AL: Not sure why we return a Macguffin vs just the struct, prob doesn't
% make a difference
everythingParams = Macguffin(basicParamsStruct); 
assert(~isempty(handles.handleobj));
handles.handleobj.data = everythingParams;

% Delete the ProjectSetup window
delete(get(hObject,'parent'));


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_copy.
function pushbutton_copy_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to pushbutton_copy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fnames,pname] = uigetfile('*.jab','Select a project to copy settings from');
if fnames == 0; return; end;

list = {'Target Type','Behavior Name and File Names','Window Features',...
  'List of perframe features','Classifier files used as input',...
  'Classifier Settings',...
  'Advanced Parameters'};

dstr = {'Select the parameters to copy.','Use Control/Command click to select multiple entries'};
[sel,ok] = listdlg('PromptString',dstr,'ListSize',[350 120],'Name','Parameters to copy',...
  'ListString',list);
sellist = list(sel);

if ok == 0, return; end

try

  origparams = load(fullfile(pname,fnames),'-mat');

if ismember('Target Type',sellist)
    handles.basicParamsStruct.featureLexiconName = origparams.x.featureLexiconName;
    handles.basicParamsStruct.featureLexicon = origparams.x.featureLexicon;
end

if ismember('Behavior Name and File Names',sellist)
  behnames = origparams.x.behaviors.names;
  isNone=strcmpi('none',behnames);
  behnames=behnames(~isNone);
  handles.basicParamsStruct.behaviors.names = behnames;
  handles.basicParamsStruct.file = origparams.x.file;
end


if ismember('List of perframe features',sellist)
  if ~isequal(handles.basicParamsStruct.featureLexicon,origparams.x.featureLexicon),
    res = questdlg(['Target type (or one of the target settings) are not the same for the current project '...
        'and the original project. Are you sure you want to import the list of perframe features?'],...
        'Import list of perframe features','Yes','No','No');
      if strcmp(res,'Yes')
        handles.basicParamsStruct.sublexiconPFNames= origparams.x.sublexiconPFNames;
      end
  else
    handles.basicParamsStruct.sublexiconPFNames= origparams.x.sublexiconPFNames;
    
  end
end

if ismember('Window Features',sellist),
  if ~isequal(handles.basicParamsStruct.featureLexicon,origparams.x.featureLexicon) || ...
    ~isequal(handles.basicParamsStruct.sublexiconPFNames,origparams.x.sublexiconPFNames)
    res = questdlg(['Target type or list of perframe features (or one of the traget settings) are not the same for the current project '...
      'and the original project. Are you sure you want to import the window features?'],...
      'Import Window features','Yes','No','No');
    if strcmp(res,'Yes')
      handles.basicParamsStruct.windowFeaturesParams = origparams.x.windowFeaturesParams;
    end
  else
      handles.basicParamsStruct.windowFeaturesParams = origparams.x.windowFeaturesParams;
  end
end



if ismember('Classifier files used as input', sellist);
  handles.basicParamsStruct.scoreFeatures = origparams.x.scoreFeatures;
end

if ismember('Classifier Settings',sellist),
  clf_params_specific = {'classifierStuff.type'
    'classifierStuff.postProcessParams'
    'classifierStuff.trainingParams'
    'classifierStuff.savewindowdata'
    'classifierStuff.selFeatures'
    'classifierStuff.predictOnlyCurrentFly'
    };
  
  for str = clf_params_specific(:)',
    try %#ok<TRYNC>
      eval(sprintf('handles.basicParamsStruct.%s = origparams.x.%s;',str{1},str{1}));
    end
  end
  
end


if ismember('Advanced Parameters',sellist),
  adv_params_specific = {'behaviors.labelcolors'
    'behaviors.unknowncolor' 
    'trxGraphicParams'
    'extra'};
  
  for str = adv_params_specific(:)',
    try %#ok<TRYNC>
      eval(sprintf('handles.basicParamsStruct.%s = origparams.x.%s;',str{1},str{1}));
    end
  end
  
end

guidata(hObject,handles);
updateEditsListboxesAndPopupmenus(handles);
updateConfigTable(handles);

catch ME,
  uiwait(warndlg(sprintf('Could not copy settings from previous project:%s',ME.message)));
end
return

% -------------------------------------------------------------------------
% --- Executes on button press in togglebutton_advanced.
function togglebutton_advanced_Callback(hObject, eventdata, handles)  %#ok
handles.mode=fif(get(hObject,'Value'),'advanced','basic');
handles = updateFigurePosition(handles);
guidata(hObject,handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_perframe.
function pushbutton_perframe_Callback(hObject, eventdata, handles)  %#ok
featureLexiconName=handles.basicParamsStruct.featureLexiconName;
sublexiconPFNames = handles.basicParamsStruct.sublexiconPFNames;
[featureLexiconPFNames,isSelected,missingPFNames] = ...
  GetAllPerframeList(featureLexiconName,sublexiconPFNames);

% If some to-be-calculated PFs are not in the lexicon, warn the user.
if ~isempty(missingPFNames),
  list = naturalLanguageListFromStringList(missingPFNames);
  if length(missingPFNames)==1
    wstr = ...
      sprintf(['Perframe feature %s is defined in the project file ' ...
               'but is not defined for the current target type.\n' ...
               'Ignoring it.'],...
              list);
  else
    wstr = ...
      sprintf(['Perframe features %s are defined in the project file ' ...
               'but are not defined for the current target type.\n' ...
               'Ignoring them.'],...
              list);
  end
  uiwait(warndlg(wstr)); 
end

% Put up the dialog that allows users to pick which PFs will be computed

[isSelected,ok] = listdlg('ListString',featureLexiconPFNames, ...
                          'InitialValue',find(isSelected), ...
                          'Name','Selecte perframe features',...
                          'PromptString',{sprintf('%s-click to',fif(ismac(),'Command','Control')), ...
                                          'select/deselect perframe features'}, ...
                                          'ListSize',[250,700]);

% If the user clicked OK, update the list of to-be-calculated PFs                 
if ok,
  sublexiconPFNames=featureLexiconPFNames(isSelected);
  handles.basicParamsStruct.sublexiconPFNames = sublexiconPFNames;
  guidata(hObject,handles);
end

return


% -------------------------------------------------------------------------
function [featureLexiconPFNames,isSelected,missingPFNames]= ...
  GetAllPerframeList(featureLexiconName,sublexiconPFNames)
% Returns the list of all per-frame feature names in the lexicon indicated
% by featureLexiconName.  Also returns isSelected, a boolean array of the
% same size as featureLexiconPFNames, which is true iff that PF is in the
% sublexiconPFNames list.  Finally, missingPFNames is a list of all PF
% names in sublexiconPFNames that are not present in
% featureLexiconPFNames.

%featureLexiconName=handles.basicParamsStruct.featureLexiconName;
%featureconfigfile = basicParamsStruct.file.featureconfigfile;
featureLexicon=featureLexiconFromFeatureLexiconName(featureLexiconName);
%basicParamsStruct = ReadXMLParams(featureconfigfile);
%allPfList = fieldnames(basicParamsStruct.perframe);
featureLexiconPFNames = fieldnames(featureLexicon.perframe);
nFeatureLexiconPFNames=length(featureLexiconPFNames);
isSelected = false(nFeatureLexiconPFNames,1);
missingPFNames = {};
%toBeCalculatedPFNames = basicParamsStruct.toBeCalculatedPFNames;
if isempty(sublexiconPFNames),
  isSelected = true(nFeatureLexiconPFNames,1);
  missingPFNames = {};
else
  for ndx = 1:length(sublexiconPFNames),
    sublexiconPFName=sublexiconPFNames{ndx};
    allndx = find(strcmp(sublexiconPFName,featureLexiconPFNames));
    if isempty(allndx)
      missingPFNames{end+1} = sublexiconPFName;  %#ok
    else
      isSelected(allndx) = true;
    end
  end
end
return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_addadvanced.
function pushbutton_addadvanced_Callback(hObject, eventdata, handles)  %#ok
in = inputdlg({'Configuration Parameter Name','Configuration Parameter Value'});
if isempty(in) || length(in) < 2
  return;
end

handles = AddConfig(handles,in{1},in{2});
updateConfigTable(handles);
updateEditsListboxesAndPopupmenus(handles);
guidata(hObject,handles);


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_removeadvanced.
function pushbutton_removeadvanced_Callback(hObject, eventdata, handles)  %#ok
data = get(handles.config_table,'Data');

jUIScrollPane = findjobj(handles.config_table);
jUITable = jUIScrollPane.getViewport.getView;
allndx = jUITable.getSelectedRows + 1;
if numel(allndx)==1 && allndx <1, return, end

for ndx = allndx(:)'
  handles = RemoveConfig(handles,data{ndx,1});
end
updateConfigTable(handles);
updateEditsListboxesAndPopupmenus(handles);
guidata(hObject,handles);
return


% -------------------------------------------------------------------------
% --- Executes when entered data in editable cell(s) in config_table.
function config_table_CellEditCallback(hObject, eventdata, handles)  %#ok
% hObject    handle to config_table (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
if isempty(eventdata.Indices); return; end
data = get(handles.config_table,'Data');
ndx = eventdata.Indices(1);
if eventdata.Indices(2) == 2
  handles = EditConfigValue(handles,data{ndx,1},data{ndx,2});
else
  handles = EditConfigName(handles,eventdata.PreviousData,eventdata.NewData);
end
updateEditsListboxesAndPopupmenus(handles);
guidata(hObject,handles);


% -------------------------------------------------------------------------
function handles = RemoveConfig(handles,name)

[fieldNamePath,dotLastField] = splitext(name);
if isempty(dotLastField)
  handles.basicParamsStruct = rmfield(handles.basicParamsStruct,fieldNamePath);
else
  lastField=dotLastField(2:end);  %#ok
  evalStr = sprintf(...
    'handles.basicParamsStruct.%s = rmfield(handles.basicParamsStruct.%s,lastField);',...
    fieldNamePath,fieldNamePath);
  eval(evalStr);
end


% -------------------------------------------------------------------------
function handles = EditConfigName(handles,oldName,newName)
eval_str = sprintf(...
  'value = handles.basicParamsStruct.%s;',...
  oldName);
eval(eval_str);
handles = AddConfig(handles,newName,value);
handles = RemoveConfig(handles,oldName);
return


% -------------------------------------------------------------------------
function handles = EditConfigValue(handles,name,value) 

switch name
  case 'behaviors.names'
    value = regexp(value,',','split'); %#ok<NASGU>
  otherwise
end

% try to convert it to a number if possible
valuen = regexp(value,',','split');
valuen = str2double(valuen);
if all(~isnan(valuen)),
  value = valuen; %#ok<NASGU>
end
eval_str = sprintf('handles.basicParamsStruct.%s = value;',name);
eval(eval_str);
return


% -------------------------------------------------------------------------
function handles = AddConfig(handles,name,value)

% If value is a string, convert to a double matrix, if that works
if ischar(value) && ~isempty(str2num(value)),  %#ok<ST2NM>
  value = str2num(value);  %#ok
end
evalStr=sprintf('handles.basicParamsStruct.%s = value;',name);
eval(evalStr);

return


% -------------------------------------------------------------------------
function figureProjectSetup_CloseRequestFcn(hObject, eventdata, handles)  %#ok
% Hint: delete(hObject) closes the figure
%uiresume(hObject);
pushbutton_cancel_Callback(hObject, eventdata, handles)
return


% -------------------------------------------------------------------------
% --- Executes on key press with focus on figureProjectSetup and none of its controls.
function figureProjectSetup_WindowKeyPressFcn(hObject, eventdata, handles)  %#ok
% hObject    handle to figureProjectSetup (see GCBO)
% eventdata  structure with the following fields (see FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

if isequal(eventdata.Key,'f') && isControlLike(eventdata.Modifier)
  handles.basicParamsStruct.file.moviefilename = 'movie.ufmf';
  handles.basicParamsStruct.file.trxfilename = 'registered_trx.mat';
  guidata(hObject,handles);  % write the changes to the figure
  updateEditsListboxesAndPopupmenus(handles);
end

return
