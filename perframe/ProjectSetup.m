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
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ProjectSetup (see VARARGIN)

% Choose default command line output for ProjectSetup
handles.output = hObject;
handles.defpath = pwd;

% Parse the arguments
[figureJLabel, ...
 basicParamsStruct] = ...
   myparse(varargin,...
           'figureJLabel',[],...
           'basicParamsStruct',[]);

% If we got called via New..., basicParamsStruct should be empty
% If via Basic Settings..., basicParamsStruct should be nonempty
handles.new=isempty(basicParamsStruct);
handles.figureJLabel=figureJLabel;
handles.basicParamsStruct = basicParamsStruct;

% Change a few things so they still work well on Mac
adjustColorsIfMac(hObject);

% Need to derive the basic and advanced sizes (part of the model) from
% the current figure dimensions
handles=setBasicAndAdvancedSizesToMatchFigure(handles);

% Set to basic mode, update the figure
handles.mode = 'basic';
handles = updateFigurePosition(handles);

% Center the window on the JLabel figure
centerOnParentFigure(hObject,figureJLabel);

% Initialize the list of possible feature lexicon names
handles.featureLexiconNameList = getFeatureLexiconListsFromXML();

% Populate the featureLexiconName popuplist with options
set(handles.featureconfigpopup,'String',handles.featureLexiconNameList);

% Either copy or 
if isempty(basicParamsStruct)
  old=warning('query','MATLAB:structOnObject');
  warning('off','MATLAB:structOnObject');  % turn off annoying warning
  handles.basicParamsStruct = struct(Macguffin('flies'));  % the default featureLexiconName
  warning(old);  % restore annoying warning  
else  
  handles.basicParamsStruct = basicParamsStruct;
end

% Update the configuration table in the GUI
updateConfigTable(handles);

% Set the window title
set(hObject,'name',fif(handles.new,'New...','Basic Settings...'));

% Set the current score-as-feature file
fileNameList = {handles.basicParamsStruct.scoreFeatures(:).classifierfile};
handles.indexOfScoreFeaturesFile=fif(isempty(fileNameList),[],1);
  
% Update handles structure
guidata(hObject, handles);

% Update all the text fields in the figure
updateEditsListboxesAndPopupmenus(handles);

% Make the figure visible
set(hObject,'Visible','on');

return


% -------------------------------------------------------------------------
% --- Outputs from this function are returned to the command line.
function ProjectSetup_OutputFcn(hObject, ~, handles)  %#ok
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
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


% % -------------------------------------------------------------------------
% function handles = initFeatureLexiconLists(handles)
% [featureLexiconNameList, ...
%  featureLexiconFileNameList, ...
%  featureLexiconAnimalTypeList] = ...
%   getFeatureLexiconListsFromXML();
% % store the three cell arrays in the handles
% handles.featureLexiconNameList=featureLexiconNameList;
% handles.featureLexiconFileNameList=featureLexiconFileNameList;
% handles.featureLexiconAnimalTypeList=featureLexiconAnimalTypeList;
% return


% % -------------------------------------------------------------------------
% function basicParams = newBasicParams(featureLexiconName)
% basicParams.featureLexiconName=featureLexiconName;
% basicParams.behaviors.names = {};
% basicParams.file.perframedir = 'perframe';
% basicParams.file.clipsdir = 'clips';
% basicParams.scoreFeatures = struct('classifierfile',{},'ts',{},'scorefilename',{});
% %basicParams.windowfeatures = struct;
% basicParams.behaviors.labelcolors = [0.7,0,0,0,0,0.7];
% basicParams.behaviors.unknowncolor = [0,0,0];
% basicParams.trxGraphicParams.colormap = 'jet';
% basicParams.trxGraphicParams.colormap_multiplier = 0.7;
% basicParams.trxGraphicParams.extra_linestyle = '-';
% basicParams.trxGraphicParams.extra_marker = '.';
% basicParams.trxGraphicParams.extra_markersize = 12;
% basicParams.labelGraphicParams.colormap = 'line';
% basicParams.labelGraphicParams.linewidth = 3;
% %basicParams.file.labelfilename = '';
% %basicParams.file.gtlabelfilename = '';
% basicParams.file.scorefilename = '';
% basicParams.file.trxfilename = '';
% basicParams.file.moviefilename = '';
% %handles = addversion(handles);
% %basicParams.scoresinput = struct('classifierfile',{},'ts',{},'scorefilename',{});
% featureLexicon=featureLexiconFromFeatureLexiconName(featureLexiconName);
% featureLexiconPFNames = fieldnames(featureLexicon.perframe);
% basicParams.sublexiconPFNames = featureLexiconPFNames;
% 
% return


% -------------------------------------------------------------------------
function updateEditsListboxesAndPopupmenus(handles)

% update the behavior name editbox
if isfield(handles.basicParamsStruct.behaviors,'names');
  names = handles.basicParamsStruct.behaviors.names;
  if numel(names)>0
    namestr = [sprintf('%s,',names{1:end-1}),names{end}];
  else
    namestr = '';
  end
  set(handles.editName,'String',namestr);
else
  set(handles.editName,'String','');
end

% Update all the editboxes
% fnames = {'labelfilename','gtlabelfilename','scorefilename',...
%           'moviefilename','trxfilename'};
% boxnames = {'editlabelfilename','editgtlabelfilename','editscorefilename',...
%             'editmoviefilename','edittrxfilename'};
fnames = {'scorefilename',...
          'moviefilename','trxfilename'};
boxnames = {'editscorefilename',...
            'editmoviefilename','edittrxfilename'};
for ndx = 1:numel(fnames)
  curf = fnames{ndx};
  curbox = boxnames{ndx};
  if isfield(handles.basicParamsStruct.file,curf);
    set(handles.(curbox),'String',handles.basicParamsStruct.file.(curf));
  else
    set(handles.(curbox),'String','');
  end
end

% Update the select feature dictionary name
indexOfFeatureLexicon=find(strcmp(handles.basicParamsStruct.featureLexiconName,handles.featureLexiconNameList));
set(handles.featureconfigpopup,'Value',indexOfFeatureLexicon);

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
[handles,whatHappened]=setFileName(handles,varName,newFileName);
warnIfInvalid(varName,whatHappened);
updateFileNameEdit(handles,editName);
updateConfigTable(handles);
return


% -------------------------------------------------------------------------
function [handles,whatHappened] = setFileName(handles,varName,newFileName)
% Set the given file name in the model, if it's a valid name.  whatHappened
% is a string that can be 'notChanged', 'emptyEntry','changed', or
% 'invalidEntry', depending.

fileName=handles.basicParamsStruct.file.(varName);
if isequal(newFileName,fileName)
  whatHappened='notChanged';
elseif isempty(newFileName)
  whatHappened='emptyEntry';
elseif IsNiceFileName(newFileName)
  whatHappened='changed';
  handles.basicParamsStruct.file.(varName)=newFileName;
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
  set(handles.(editName),'String',handles.basicParamsStruct.file.(varName));
return


% -------------------------------------------------------------------------
function updateConfigTable(handles)
% Update the config table (a GUI element) to match the current "model"
% state
basicParamsStruct = handles.basicParamsStruct;
fields2remove = {'featureLexicon','windowfeatures','scoreFeatures', ...
                 'sublexiconPFNames','labels','gtLabels','expDirNames', ...
                 'gtExpDirNames', 'classifierStuff', 'version'};
for ndx = 1:numel(fields2remove)
  if isfield(basicParamsStruct,fields2remove{ndx}),
    basicParamsStruct = rmfield(basicParamsStruct,fields2remove{ndx});
  end  
end
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
      q = num2str(param(1));
      for jj = 2:numel(param)
        q = [q ',' num2str(param(jj))]; %#ok<AGROW>
      end
      list{end,2} = q;
    else
      list{end,2} = param;
    end
  end
end


% % -------------------------------------------------------------------------
% function handles = addversion(handles)
% if ~isfield(handles.basicParams,'ver')
%   vid = fopen('version.txt','r');
%   vv = textscan(vid,'%s');
%   fclose(vid);
%   handles.basicParams.ver = vv{1};
% end


% -------------------------------------------------------------------------
function editName_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editName as text
%        str2double(get(hObject,'String')) returns contents of editName as a double
name = get(hObject,'String');
if isempty(regexp(name,'^[a-zA-Z][\w_,]*$','once','start'));
   uiwait(warndlg(['The behavior name cannot have special characters.'...
                   'Please use only alphanumeric characters and _']));
   return;
end
    
name = regexp(name,',','split');
name_str = [sprintf('%s_',name{1:end-1}),name{end}];
handles.basicParamsStruct.behaviors.names = name;
%handles.basicParamsStruct.file.moviefilename = 'movie.ufmf';
%handles.basicParamsStruct.file.trxfilename = 'registered_trx.mat';
%handles.basicParamsStruct.file.labelfilename = sprintf('label_%s.mat',name_str);
%handles.basicParamsStruct.file.gtlabelfilename = sprintf('gt_label_%s.mat',name_str);
handles.basicParamsStruct.file.scorefilename = sprintf('scores_%s.mat',name_str);
guidata(hObject,handles);
updateEditsListboxesAndPopupmenus(handles);
updateConfigTable(handles);
return


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function editName_CreateFcn(hObject, eventdata, handles)  %#ok
% hObject    handle to editName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
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
handles.basicParamsStruct.featureLexiconName=featureLexiconNameNew;
%animalTypeNew=handles.featureLexiconAnimalTypeList{iFeatureLexicon};
%handles.basicParamsStruct.behaviors.type=animalTypeNew;
guidata(hObject,handles);
updateConfigTable(handles);
return


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function featureconfigpopup_CreateFcn(hObject, ~, ~)  %#ok
% hObject    handle to featureconfigpopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
function editlabelfilename_Callback(hObject, eventdata, handles) %#ok
% hObject    handle to editlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editlabelfilename as text
%        str2double(get(hObject,'String')) returns contents of editlabelfilename as a double
editName=get(hObject,'tag');
handles = fileNameEditTwiddled(handles,editName);
guidata(hObject,handles);


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function editlabelfilename_CreateFcn(hObject, eventdata, handles) %#ok
% hObject    handle to editlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
function editgtlabelfilename_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to editgtlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editgtlabelfilename as text
%        str2double(get(hObject,'String')) returns contents of editgtlabelfilename as a double
editName=get(hObject,'tag');
handles = fileNameEditTwiddled(handles,editName);
guidata(hObject,handles);


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function editgtlabelfilename_CreateFcn(hObject, eventdata, handles)  %#ok
% hObject    handle to editgtlabelfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
function editscorefilename_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to editscorefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editscorefilename as text
%        str2double(get(hObject,'String')) returns contents of editscorefilename as a double
editName=get(hObject,'tag');
handles = fileNameEditTwiddled(handles,editName);
guidata(hObject,handles);


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function editscorefilename_CreateFcn(hObject, eventdata, handles)  %#ok
% hObject    handle to editscorefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
function editmoviefilename_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to editmoviefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editmoviefilename as text
%        str2double(get(hObject,'String')) returns contents of editmoviefilename as a double
editName=get(hObject,'tag');
handles = fileNameEditTwiddled(handles,editName);
guidata(hObject,handles);


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function editmoviefilename_CreateFcn(hObject, eventdata, handles)  %#ok
% hObject    handle to editmoviefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
function edittrxfilename_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to edittrxfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edittrxfilename as text
%        str2double(get(hObject,'String')) returns contents of edittrxfilename as a double
editName=get(hObject,'tag');
handles = fileNameEditTwiddled(handles,editName);
guidata(hObject,handles);


% -------------------------------------------------------------------------
% --- Executes during object creation, after setting all properties.
function edittrxfilename_CreateFcn(hObject, eventdata, handles)  %#ok
% hObject    handle to edittrxfilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% -------------------------------------------------------------------------
% --- Executes on butuiresume(handles.figureProjectSetup);ton press in pushbutton_setfeatures.
function pushbutton_setfeatures_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to pushbutton_setfeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_copyfeatures.
function pushbutton_copyfeatures_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to pushbutton_copyfeatures (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% % -------------------------------------------------------------------------
% % --- Executes on selection change in listbox_inputscores.
% function listbox_inputscores_Callback(hObject, eventdata, handles)
% % hObject    handle to listbox_inputscores (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% handles.indexOfScoreFeaturesFile=get(hObject,'Value');
% guidata(hObject,handles);
% return

% % -------------------------------------------------------------------------
% % --- Executes during object creation, after setting all properties.
% function listbox_inputscores_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to listbox_inputscores (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: listbox controls usually have a white background on Windows.
% %       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% % -------------------------------------------------------------------------
% % --- Executes on button press in pushbutton_addlist.
% function pushbutton_addlist_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton_addlist (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% [fname,pname] = ...
%   uigetfile({'*.jab','JAABA Files (*.jab)'}, ...
%             'Add .jab file containing classifier to be used as input');
% if fname == 0; return; end;
% 
% fileNameAbs = fullfile(pname,fname);
% everythingParams = load(fileNameAbs,'-mat');
% if isempty(everythingParams.classifier)
%   uiwait(errordlg(sprintf('%s does not contain a classifier.',fname), ...
%                   'Error', ...
%                   'modal'));
%   return
% end               
% 
% % [fnames,pname] = uigetfile('*.mat','Classifier whose scores should be used as input');
% % if fnames == 0; return; end;
% % 
% % cfile = fullfile(pname,fnames);
% % classifier = load(cfile);
% 
% curs = struct;  % this is the thing we'll return
% 
% % Add the .jab file name, with maybe not the best field name
% curs.classifierfile = fileNameAbs;
% 
% % Add the classifier time stamp
% classifier=everythingParams.classifier;
% if isfield(classifier,'timeStamp');
%   curs.ts = classifier.timeStamp;
% else
%   uiwait(errordlg('The classifier in the selected field lacks a timestamp.  Aborting.', ...
%                   'Error', ...
%                   'modal'));
%   return
% end
% 
% % Add the name of the score file (without the .mat extension)
% if isfield(classifier,'file') && isfield(classifier.file,'scorefilename')
%   scorefilename = classifier.scorefilename;
%   [~,scoreBaseName] = fileparts(scorefilename);
%   curs.scorefilename = scoreBaseName;
% elseif isfield(everythingParams,'behaviors') && ...
%        isfield(everythingParams.behaviors,'names') && ...
%        ~isempty(everythingParams.behaviors.names)
%   behaviorName=everythingParams.behaviors.names{1};
%   curs.scorefilename = sprintf('scores_%s',behaviorName);
% else
%   uiwait(errordlg('Unable to determine score file name.  Aborting.', ...
%                   'Error', ...
%                   'modal'));
%   return
% end
% 
% handles.basicParams.scoreFeatures(end+1) = curs;
% handles.indexOfScoreFeaturesFile = length(handles.basicParams.scoreFeatures);
% updateEditsListboxesAndPopupmenus(handles);
% guidata(hObject,handles);
% return


% % -------------------------------------------------------------------------
% % --- Executes on button press in pushbutton_removelist.
% function pushbutton_removelist_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton_removelist (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% i = handles.indexOfScoreFeaturesFile;
% if isempty(i), return; end
% nScoreFeatures=length(handles.basicParams.scoreFeatures);
% handles.basicParams.scoreFeatures(i) = [];
% if (i==nScoreFeatures)
%   i=i-1;
%   if i==0
%     i=[];
%   end
% end
% handles.indexOfScoreFeaturesFile=i;
% guidata(hObject,handles);
% updateEditsListboxesAndPopupmenus(handles);
% return

% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to pushbutton_cancel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.outfile = 0;
guidata(hObject,handles);
delete(findAncestorFigure(hObject));              
%uiresume(handles.figureProjectSetup);
return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_done.
function pushbutton_done_Callback(hObject, ~, handles)  %#ok
% hObject    handle to pushbutton_done (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check for an empty behavior name
behaviorName=handles.basicParamsStruct.behaviors.names;
if isempty(behaviorName) ,
  uiwait(errordlg('You must enter a behavior name.', ...
                  'No behavior name', ...
                  'modal'));
  return
end

% Check for an empty movie file name
movieFileName=handles.basicParamsStruct.file.moviefilename;
if isempty(movieFileName) ,
  uiwait(errordlg('You must enter a movie file name.', ...
                  'No movie file name', ...
                  'modal'));
  return
end

% Check for an empty track file name
trackFileName=handles.basicParamsStruct.file.trxfilename;
if isempty(trackFileName) ,
  uiwait(errordlg('You must enter a track file name.', ...
                  'No track file name', ...
                  'modal'));
  return
end

% Check for an empty score file name
scoreFileName=handles.basicParamsStruct.file.scorefilename;
if isempty(scoreFileName) ,
  uiwait(errordlg('You must enter a score file name.', ...
                  'No score file name', ...
                  'modal'));
  return
end

% Get the info we need out of the handles
basicParamsStruct=handles.basicParamsStruct;
figureJLabel=handles.figureJLabel;
everythingParams=Macguffin(basicParamsStruct);

% Call the appropriate function to notify the JLabel "object" that 
% we're done.
if (handles.new)
  % we were called via New..., so act accordingly
  JLabel('newFileSetupDone', ...
         figureJLabel, ...
         everythingParams);
else
  % we were called via Basic Settings..., so act accordingly
  JLabel('basicSettingsChanged', ...
         figureJLabel, ...
         everythingParams); 
end

% Delete the ProjectSetup window
delete(get(hObject,'parent'));

return


% % -------------------------------------------------------------------------
% % --- Executes on button press in pushbutton_copy.
% function pushbutton_copy_Callback(hObject, eventdata, handles)
% % hObject    handle to pushbutton_copy (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% [fnames,pname] = uigetfile('*.mat','Select a project to copy settings from');
% if fnames == 0; return; end;
% 
% list = {'Target Type','Behavior Name and File Names','Window Features',...
%   'List of perframe features','Classifier files used as input',...
%   'Advanced Parameters'};
% 
% dstr = {'Select the parameters to copy.','Use Control/Command click to select multiple entries'};
% [sel,ok] = listdlg('PromptString',dstr,'ListSize',[350 120],'Name','Parameters to copy',...
%   'ListString',list);
% sellist = list(sel);
% 
% if ok == 0, return; end
% [~,~,ext] = fileparts(fnames);
% if strcmp(ext,'.xml')
%   origparams = ReadXMLParams(fullfile(pname,fnames));
% else
%   origparams = load(fullfile(pname,fnames));
% end
% 
% if ismember('Target Type',sellist)
%   alltypes = fieldnames(handles.featureLexiconFileLookup);
%   listndx = [];
%   for andx = 1:numel(alltypes)
%     curfname = handles.featureLexiconFileLookup.(alltypes{andx}).file;
%     if strcmp(origparams.file.featureconfigfile,curfname); 
%       listndx = andx;
%     end
%   end
%   if ~isempty(listndx),
%     set(handles.featureconfigpopup,'Value',listndx);
%   else
%     uiwait(warndlg(['The feature config file used (', origparams.file.featureconfigfile,...
%       ') does not exist. Cannot import from the project']));
%     return;
%   end
%   handles.basicParams.file.featureconfigfile = origparams.file.featureconfigfile;
% end
% 
% 
% if ismember('Behavior Name and File Names',sellist)
%   origfeatureconfigfile = handles.basicParams.file.featureconfigfile;
%   handles.basicParams.file = origparams.file;
%   handles.basicParams.file.featureconfigfile = origfeatureconfigfile;
%   if iscell(origparams.behaviors.names),
%     handles.basicParams.behaviors.names = origparams.behaviors.names;
%   else
%     handles.basicParams.behaviors.names = {origparams.behaviors.names};
%   end
%   if ~isfield(origparams.file,'scorefilename')
%      name = handles.basicParams.behaviors.names;
%      name_str = [sprintf('%s_',name{1:end-1}),name{end}];
%      handles.basicParams.file.scorefilename =  sprintf('scores_%s.mat',name_str);
%   end
% end
% 
% if ismember('Window Features',sellist),
%   if isfield(origparams,'windowfeatures')
%     
%     if ~strcmp(handles.basicParams.file.featureconfigfile,origparams.file.featureconfigfile),
%       res = questdlg(['Target type are not the same for the current project '...
%         'and the original project. Are you sure you want to import the window features?'],...
%         'Import Window features','Yes','No','No');
%       if strcmp(res,'Yes')
%         handles.basicParams.windowfeatures = origparams.windowfeatures;
%       end
%     else
%       handles.basicParams.windowfeatures = origparams.windowfeatures;
%     end
%   elseif isfield(origparams.file,'featureparamfilename')
%     uiwait(warndlg(['The selected configuration file does not have any '...
%       'window features, but points to a file that does have the window features. '...
%       'Loading the window features from the referenced file:' origparams.file.featureparamfilename]));
%     [windowfeaturesparams,windowfeaturescellparams,basicFeatureTable,featureWindowSize] = ...
%       ReadPerFrameParams(origparams.file.featureparamfilename,handles.basicParams.file.featureconfigfile); 
% 
%     handles.basicParams.windowfeatures.windowfeaturesparams = windowfeaturesparams;
%     handles.basicParams.windowfeatures.windowfeaturescellparams = windowfeaturescellparams;
%     handles.basicParams.windowfeatures.basicFeatureTable = basicFeatureTable;
%     handles.basicParams.windowfeatures.featureWindowSize = featureWindowSize;
%     if isfield(handles.basicParams.file,'featureparamfilename'),
%         handles.basicParams.file = rmfield(handles.basicParams.file,'featureparamfilename');
%     end
%   else
%     uiwait(warndlg(['The selected configuration file does not have any '...
%       'window features. Not copying the window features']));
%   end
% end
% 
% if ismember('List of perframe features',sellist)
%   
%   if isfield(origparams,'featureparamlist')
%     
%     if ~strcmp(handles.basicParams.file.featureconfigfile,origparams.file.featureconfigfile),
%       res = questdlg(['Target type are not the same for the current project '...
%         'and the original project. Are you sure you want to import the list of perframe features?'],...
%         'Import list of perframe features','Yes','No','No');
%       if strcmp(res,'Yes')
%         handles.basicParams.featureparamlist = origparams.featureparamlist;
%       end
%     else
%       handles.basicParams.featureparamlist= origparams.featureparamlist;
%       
%     end
%   end
% end
% 
% if ismember('Classifier files used as input', sellist);
%   if isfield(origparams,'scoreFeatures')
%     
%     handles.basicParams.scoreFeatures = origparams.scoreFeatures;
%   end
% end
% 
% if ismember('Advanced Parameters',sellist),
%   adv_params = {'behaviors.labelcolors',...
%     'behaviors.unknowncolor',...
%     'plot.trx.colormap',...
%     'plot.trx.colormap_multiplier',...
%     'plot.trx.nextra_markers',...
%     'plot.trx.extra_marker',...
%     'plot.trx.extra_markersize',...
%     'plot.trx.extra_linestyle',...
%     'plot.labels.colormap',...
%     'plot.labels.linewidth',...
%     'perframe.basicParams.fov',...
%     'perframe.basicParams.nbodylengths_near',...
%     'perframe.basicParams.thetafil',...
%     'perframe.basicParams.max_dnose2ell_anglerange',...
%     'perframe.landmark_params.arena_center_mm_x',...
%     'perframe.landmark_params.arena_center_mm_y',...
%     'perframe.landmark_params.arena_radius_mm',...
%     'perframe.landmark_params.arena_type'};
%   
%   for str = adv_params(:)',
%     try %#ok<TRYNC>
%       eval(sprintf('handles.basicParams.%s = origparams.%s;',str{1},str{1}));
%     end
%   end
% 
% end
% 
% handles = addversion(handles);
% guidata(hObject,handles);
% updateEditsListboxesAndPopupmenus(handles);
% updateConfigTable(handles);
% return

% -------------------------------------------------------------------------
% --- Executes on button press in togglebutton_advanced.
function togglebutton_advanced_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to togglebutton_advanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebutton_advanced
handles.mode=fif(get(hObject,'Value'),'advanced','basic');
handles = updateFigurePosition(handles);
guidata(hObject,handles);
return


% -------------------------------------------------------------------------
% --- Executes on button press in pushbutton_perframe.
function pushbutton_perframe_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to pushbutton_perframe (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the boolean array that tells which PFs in the feature lexicon are in 
% the to-be-calculated feature vocabulary
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
                          'PromptString',{'Control/Command click to', ...
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
% hObject    handle to pushbutton_addadvanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
% hObject    handle to pushbutton_removeadvanced (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
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
function handles = EditConfigValue(handles,name,value)  %#ok
eval_str = sprintf('handles.basicParamsStruct.%s = value;',name);
eval(eval_str);
return


% -------------------------------------------------------------------------
function handles = AddConfig(handles,name,value)

% If value is a string, convert to a double matrix, if that works
if ischar(value) && ~isempty(str2num(value)),  %#ok<ST2NM>
  value = str2num(value);  %#ok
end

% The code below doesn't seem to matter, since curstruct is never passed
% out of the function, nor is anything that depends upon it.
% iname = fliplr(name);
% curstruct = handles.basicParamsStruct;
% while true,
%   [iname,lastfield] = splitext(iname);
%   if isempty(lastfield)
%     fexist = isfield(curstruct,iname);
%     break;
%   else
%     fexist = isfield(curstruct,fliplr(iname(2:end)));
%     if ~fexist, break;    end
%     curstruct = curstruct.(fliplir(iname(2:end)));
%   end
% end

evalStr=sprintf('handles.basicParamsStruct.%s = value;',name);
eval(evalStr);

return


% -------------------------------------------------------------------------
% --- Executes when user attempts to close figureProjectSetup.
function figureProjectSetup_CloseRequestFcn(hObject, eventdata, handles)  %#ok
% hObject    handle to figureProjectSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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
