function varargout = APTProject(varargin)
% APTPROJECT MATLAB code for APTProject.fig
%      APTPROJECT, by itself, creates a new APTPROJECT or raises the existing
%      singleton*.
%
%      H = APTPROJECT returns the handle to a new APTPROJECT or the handle to
%      the existing singleton*.
%
%      APTPROJECT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APTPROJECT.M with the given input arguments.
%
%      APTPROJECT('Property','Value',...) creates a new APTPROJECT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before APTProject_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to APTProject_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help APTProject

% Last Modified by GUIDE v2.5 07-Feb-2020 07:35:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @APTProject_OpeningFcn, ...
                   'gui_OutputFcn',  @APTProject_OutputFcn, ...
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


% --- Executes just before APTProject is made visible.
function APTProject_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to APTProject (see VARARGIN)

% Choose default command line output for APTProject
handles.output = hObject;

[lbl_file,aptStruct] = myparse(varargin,'lbl_file','','aptStruct',struct());

aptStruct.lbl_file = lbl_file;
if ~isfield(aptStruct,'featureLexicon')
  aptStruct.featureLexicon = [];
  aptStruct.featureLexiconName = '';
end

if ~isfield(aptStruct,'origFeatureLexiconName')
  aptStruct.origFeatureLexiconName = '';
end
if isfield(aptStruct,'animalType'),
else
  aptStruct.animalType = '';
end
if isfield(aptStruct,'n_pts'),
else
  aptStruct.n_pts = 0;
end
if isfield(aptStruct,'n_view'),
else
  aptStruct.n_view = 0;
end
if isfield(aptStruct,'projname'),
else
  aptStruct.projname = '';
end
if isfield(aptStruct,'trkfilename'),
else
  aptStruct.trkfilename = '';
end

handles.aptStruct = aptStruct;
% Update handles structure
guidata(hObject, handles);

% Checks
if isempty(lbl_file)
  % close the gui if there is no label file.
  % the gui is closed in outputfcn.
  return
end

h = waitbar(0,'Loading label file');
handles.waitbar = h;
guidata(hObject,handles);
try
  apt = loadLbl(lbl_file);
catch ME
  errmsg = sprintf('Could not load label file %s',ME.message);
  uiwait(errordlg(errmsg, 'Unable to load label file'));
  handles.aptStruct.lbl_file = '';
  guidata(hObject,handles);
  return;
end

if isempty(apt.movieFilesAll)
  uiwait(errordlg('APT Label file should have at least one movie','Empty APT project'));
  return;
end

%
handles.apt = apt;
handles.has_trx = apt.projectHasTrx;
if isfield(apt,'maIsMA')
  handles.is_ma = apt.maIsMA;
else
  handles.is_ma = false;
end
handles.has_crops = ~isempty(apt.movieFilesAllCropInfo{1});
aptStruct.n_pts = apt.cfg.NumLabelPoints;
aptStruct.skeletonEdges = apt.skeletonEdges;
if handles.is_ma
  if ~isempty(apt.skelHead)
    handles.head_tail = [apt.skelHead apt.skelTail];
    handles.use_theta = true;
  else
    handles.head_tail = [];
    handles.use_theta = false;
  end
  sz = apt.trackParams.ROOT.MultiAnimal.TargetCrop.ManualRadius;
  
%  assert 'APT multi-animal projects not yet implemented';
elseif handles.has_trx
  if isfield(apt.trackParams.ROOT.ImageProcessing.MultiTarget,'TargetCrop')
    handles.use_theta = apt.trackParams.ROOT.ImageProcessing.MultiTarget.TargetCrop.AlignUsingTrxTheta;
    sz = apt.trackParams.ROOT.ImageProcessing.MultiTarget.TargetCrop.Radius;
  else
    handle.use_theta = apt.trackParams.ROOT.MultiAnimal.TargetCrop.AlignUsingTrxTheta;
    sz = apt.trackParams.ROOT.MultiAnimal.TargetCrop.ManualRadius;
  end
  handles.sz = [sz,sz];
elseif handles.has_crops
  % dont use theta if there is no trx file. 
  % might need updating later.
  handles.use_theta = true;
  ff = apt.movieFilesAllCropInfo{1}(1).roi;
  ht = ff(4)-ff(3)+1; wd = ff(2)-ff(1)+1;
  handles.sz = [ht,wd];
else
  handles.use_theta = true;
  ht = apt.movieInfoAll{1}.info.nr;
  wd = apt.movieInfoAll{1}.info.nc;
  handles.sz = [ht,wd];
end

[featureList,xmlList] = getFeatureLexiconListsFromXML();
stFeatures = cellfun(@xmlParamsIsST,xmlList);
featureList = featureList(~stFeatures);
xmlList = xmlList(~stFeatures);
handles.trxList = {featureList, xmlList};

% handles.global_ftrs = {'velmag'};
% handles.body_ftrs = {'velmag',...
%           'distcenter', 'ddistcenter'};
% handles.theta_ftrs = {'x','y','dx','dy',...
%           'dtheta', 'sin', 'cos'};
handles.common_ftrs = {'global_velmag','body_velmag',...
          'body_distcenter', 'body_ddistcenter',...
          'body_x','body_y','body_dx','body_dy',...
          'body_dtheta', 'body_sin', 'body_cos'};
handles.pair_ftrs = {'velmag',...
          'ddist','dist', ... % These are same body_ftrs
          'x','y','dx','dy',...
          'dtheta','cos','sin',... % these are same theta_ftrs
          'areaswept'};
handles.triad_ftrs = {'cos','sin','dangle',...
          'area','darea','dlen'};

waitbar(.5,handles.waitbar,'Reading example frame from movie file');        
handles.aptStruct = aptStruct;
handles = initialize(handles);
delete(handles.waitbar);
guidata(hObject,handles);

% UIWAIT makes APTProject wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = APTProject_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles.aptStruct.lbl_file)
  uiwait(handles.figure1);
end
if ~isgraphics(hObject)
  aptStruct = struct();
  aptStruct.featureLexicon = [];
  varargout{1} = aptStruct;
  return;
end
handles = guidata(hObject);
varargout{1} = handles.aptStruct;
if isgraphics(handles.waitbar)
  delete(handles.waitbar);
end
figure1_CloseRequestFcn(hObject,eventdata,handles);


function handles = initialize(handles)
is_multi = handles.has_trx|handles.is_ma;
set(handles.trx_pop,'Enable',onoffif(is_multi));
set(handles.trx_pop,'Visible',onoffif(is_multi));
set(handles.trx_checkbox,'Visible',onoffif(is_multi));
set(handles.trx_checkbox,'Value',is_multi);
set(handles.group_notrx,'Visible',onoffif(~is_multi));
set(handles.rb_centroid,'Visible',onoffif(~is_multi));
set(handles.rb_frame,'Visible',onoffif(~is_multi));
set(handles.rb_crop,'Visible',onoffif(~is_multi));
set(handles.rb_custom,'Visible',onoffif(~is_multi));
set(handles.text_centroid,'Visible',onoffif(~is_multi));

if is_multi
  set(handles.trx_pop,'String',handles.trxList{1})
  trx_type_ndx = find(strcmp(handles.aptStruct.origFeatureLexiconName,handles.trxList{1}),1);
  if ~isempty(trx_type_ndx),
    set(handles.trx_pop,'Value',trx_type_ndx);
  else
    set(handles.trx_pop,'Value',1);
  end
  set(handles.trx_checkbox,'Value',true);
else
  set(handles.trx_pop,'String',{''})
  set(handles.trx_pop,'Value',0)
  set(handles.trx_checkbox,'Value',true);
  if handles.has_crops
    set(handles.rb_crop,'Value',true);
    set(handles.rb_frame,'Value',false);
    set(handles.rb_crop,'Enable','on');
  else
    set(handles.rb_crop,'Value',false);
    set(handles.rb_crop,'Enable','off');
    set(handles.rb_frame,'Value',true);
  end
  set(handles.rb_centroid,'Value',false);
end

% load a reference image and point locations.
apt = handles.apt;
has_ref_img = false;
has_mov = true;
for ndx = 1:size(apt.movieFilesAll,1)
  has_label = false;
  tndx_to_use = 1;
  locs = [];
  if isfield(apt,'labeledpos')
    s = apt.labeledpos{ndx};    
    cur_pts = nan(s.size);
    cur_pts(s.idx) = s.val; % scalar expansion for 'log'
    for tndx = 1:size(cur_pts,4)
      frm = find(~isnan(cur_pts(1,1,:,tndx)),1);
      if isempty(frm)
        continue;
      else
        has_label = true;
        tndx_to_use = tndx;
        locs = cur_pts(:,:,frm,tndx);
        break;
      end
    end
  else
    s = apt.labels{ndx};
    assert(isstruct(s));
    cur_pts = s.p;
    if numel(s.frm)>0
      has_label = true;
      frm = s.frm(1);
      tndx_to_use = s.tgt(1);
      locs = reshape(s.p(:,1),[s.npts,2]);
    else
      has_label = false;
    end
  end
  
  if ~has_label
    continue;
  end
  
  mov_files = {};
  trx_files = {};
  for m_ndx = 1:size(apt.movieFilesAll,2)
    mov_file = unMacroise(apt.movieFilesAll{ndx,m_ndx},apt);
    if ~exist(mov_file,'file')      
      [bdir,mname,mext] =fileparts(mov_file);
      qstr = sprintf('JAABA needs movie file %s%s to show a reference labeled example, but it does not exist at location %s. Browse to movie? If you select no, no reference image will be shown.',mname,mext,bdir);
      res = questdlg(qstr,'Missing movie file...','Yes','No','Cancel','Yes');
      if strcmp(res,'Yes')
        [fname,fpath] = uigetfile('',['Missing ' mov_file]);
        if ~isequal(fname,0)
          mov_files{end+1} = fullfile(fpath,fname);
          if handles.has_trx
            trx_file = getTrxFile(apt.trxFilesAll{ndx},mov_file);
            [fname,fpath] = uigetfile('',['Missing ' trx_file]);
            trx_files{end+1} = fullfile(fpath,fname);
          end
        else
          has_mov = false;
          break;
        end
      else
        has_mov = false;
        break;
      end
    else
      mov_files{end+1} = mov_file;
      trx_files{end+1} = unMacroise(apt.trxFilesAll{ndx,m_ndx},apt);
    end
  end
  
  if ~has_mov
    break;
  end
  
  has_ref_img = true;
  img = []; all_locs = []; prev_width = 0;
  for m_ndx = 1:size(apt.movieFilesAll,2)
    mov_file = mov_files{m_ndx};
    if has_mov
      rfn = get_readframe_fcn(mov_file);
      cur_img = rfn(frm);
    else
      nr = apt.movieInfoAll{ndx,m_ndx}.info.nr;
      nc = apt.movieInfoAll{ndx,m_ndx}.info.nc;
      cur_img = zeros(nr,nc);
    end
    if handles.has_trx
      trx = load(trx_files{m_ndx},'-mat');
      trx = trx.trx(tndx_to_use);
      off = trx.off;
      x = trx.x(frm+off);
      y = trx.y(frm+off);
      theta = trx(tndx_to_use).theta(frm+off);
      [cur_img, aff_mat] = CropImAroundTrx(cur_img,x,y,theta,handles.sz(2),handles.sz(1));
      ll = locs*aff_mat(1:2,1:2);
      ll(:,2) = ll(:,2)+aff_mat(3,2) + handles.sz(1) + 0.5;
      ll(:,1) = ll(:,1)+aff_mat(3,1) + handles.sz(2) + 0.5;
      ll(:,1) = ll(:,1) + prev_width;
    elseif handles.is_ma
      if ~isempty(handles.head_tail)
        h_pt = locs(handles.head_tail(1),:);
        t_pt = locs(handles.head_tail(2),:);
        x = (h_pt(1)+t_pt(1))/2;
        y = (h_pt(2)+t_pt(2))/2;
        dy = (t_pt(2)-h_pt(2));
        dx = (t_pt(1)-h_pt(1));
        theta = atan2( dy,dx );
      else
        x = mean(locs(:,1));
        y = mean(locs(:,2));
        theta = 0;
      end
      [cur_img, aff_mat] = CropImAroundTrx(cur_img,x,y,theta,handles.sz(2),handles.sz(1));
      ll = locs*aff_mat(1:2,1:2);
      ll(:,2) = ll(:,2)+aff_mat(3,2) + handles.sz(1) + 0.5;
      ll(:,1) = ll(:,1)+aff_mat(3,1) + handles.sz(2) + 0.5;
      ll(:,1) = ll(:,1) + prev_width;      
    else
      n_pts = handles.aptStruct.n_pts;
      st = (m_ndx-1)*n_pts + 1; en = m_ndx*n_pts;
      ll = locs(st:en,:);
      ll(:,1) = ll(:,1) + prev_width;
    end
    img = cat(2,img,cur_img);
    prev_width = size(img,2);
    all_locs = cat(1,all_locs,ll);
  end
  
  break
end

if has_ref_img
  handles.im = img;
  handles.locs = all_locs;
else
  handles.im = [];
  handles.locs = [];
end

set(handles.text_centroid,'String','');

npts = apt.cfg.NumLabelPoints;
if isfield(handles.aptStruct,'pairs') && ~isempty(handles.aptStruct.pairs),
  isvalid = all(cellfun(@(x) all(x==round(x)) && numel(x) == 2 && x(1) ~= x(2) && ...
    all(x>=1) && all(x<=npts),handles.aptStruct.pairs));
  if isvalid,
    pairstr = sprintf('%d %d, ',handles.aptStruct.pairs{:});
    pairstr = pairstr(1:end-2);
    set(handles.edit_pair,'String',pairstr);
  end
end

if isfield(handles.aptStruct,'triads') && ~isempty(handles.aptStruct.triads),
  isvalid = all(cellfun(@(x) all(x==round(x)) && numel(x) == 3 && numel(unique(x))==3 && ...
    all(x>=1) && all(x<=npts),handles.aptStruct.triads));
  if isvalid,
    str = sprintf('%d %d %d, ',handles.aptStruct.triads{:});
    str = str(1:end-2);
    set(handles.edit_triad,'String',str);
  end
end

if isfield(handles.aptStruct,'custom_list') && ~isempty(handles.aptStruct.custom_list),
  str = sprintf('%s, ',handles.aptStruct.custom_list{:});
  str = str(1:end-2);
  set(handles.edit_custom,'String',str);
end

isvalid = ~isempty(handles.aptStruct.trkfilename) && ...
  isvalid_default_trkfilename(handles,handles.aptStruct.trkfilename);
if isvalid,
  str = sprintf('%s, ',handles.aptStruct.trkfilename{:});
  str = str(1:end-2);
  set(handles.edit_trk,'String',str);
else
  trkfilename = reset_trk_edit(handles);
  handles.aptStruct.trkfilename = trkfilename;
end

handles.aptStruct.n_pts = apt.cfg.NumLabelPoints;
handles.aptStruct.n_view = apt.cfg.NumViews;
handles.aptStruct.projname = apt.projname;
handles.aptStruct.head_tail = handles.head_tail;
% handles.aptStruct.featureLexicon = [];
% handles.aptStruct.featureLexiconName = '';
% handles.aptStruct.animalType = '';
handles.aptStruct.view = 1;
handles = UpdateHandles(handles);
UpdatePlots(handles)



function handles = UpdateHandles(handles)
handles.use_trx = get(handles.trx_checkbox,'Value');
if handles.has_trx
  trx_type_ndx = get(handles.trx_pop,'Value');
  featureLexiconName = handles.trxList{1}{trx_type_ndx};
  [featureLexicon,animalType]=featureLexiconFromFeatureLexiconName(featureLexiconName);
else
  featureLexiconName = 'apt';
  animalType = 'apt';
end
handles.aptStruct.origFeatureLexiconName = featureLexiconName;
if handles.has_trx
  handles.aptStruct.featureLexiconName = [featureLexiconName '_' handles.aptStruct.projname];
else
  handles.aptStruct.featureLexiconName = featureLexiconName;
end

handles.aptStruct.animalType = animalType;
% pairs
str = get(handles.edit_pair,'String');
pair_str = strsplit(str,',');
pairs = {};
if ~isempty(pair_str{1}),
  for ndx = 1:numel(pair_str)
    cur_str = strtrim(pair_str{ndx});
    cur_pair = strsplit(cur_str);
    if numel(cur_pair)~=2
      uiwait(errordlg('Incorrect format for specifying pairs.'));
      pairs = {}; break;
    end
    pts = [];
    is_error = false;
    for pndx = 1:numel(cur_pair)
      cur_p = str2double(cur_pair{pndx});
      if round(cur_p) ~= cur_p || isnan(cur_p)
        uiwait(errordlg('Please specify points as integer'));
        pairs = {}; is_error = true; break;
      end
      if cur_p > handles.aptStruct.n_pts
        uiwait(errordlg(sprintf('Points should be between 1 and %d',handles.n_pts)));
        pairs = {};  is_error = true;break;
      end
      pts(end+1) = cur_p;
    end
    if is_error
      pairs = {}; break;
    end
    if pts(1) == pts(2)
      uiwait(errordlg('Cannot use same point in pairs'));
      pairs = {}; break;
    end
    pairs{end+1} = pts;
  end
end
handles.pairs = pairs;

% triads
str = get(handles.edit_triad,'String');
pair_str = strsplit(str,',');
triads = {};
if ~isempty(pair_str{1}),
  for ndx = 1:numel(pair_str)
    cur_str = strtrim(pair_str{ndx});
    cur_pair = strsplit(cur_str);
    if numel(cur_pair)~=3
      uiwait(errordlg('Incorrect format for specifying triads.'));
      triads = {}; break;
      return;
    end
    pts = [];
    is_error = false;
    for pndx = 1:numel(cur_pair)
      cur_p = str2double(cur_pair{pndx});
      if round(cur_p) ~= cur_p || isnan(cur_p)
        uiwait(errordlg('Please specify points as integer'));
        triads = {}; is_error = true; break;
      end
      if cur_p > handles.aptStruct.n_pts
        uiwait(errordlg(sprintf('Points should be between 1 and %d',handles.aptStruct.n_pts)));
        triads = {};  is_error = true;break;
      end
      pts(end+1) = cur_p;
    end
    if is_error, triads = {}; break; end
    if numel(unique(pts))~=3
      uiwait(errordlg('Cannot use same points twice in the triads'));
      triads = {}; break;
    end
    triads{end+1} = pts;
  end
end
handles.triads = triads;

% custom
str = get(handles.edit_custom,'String');
pair_str = strsplit(str,',');
custom_list = {};
if ~isempty(pair_str{1}),
  for ndx = 1:numel(pair_str)
    custom_list{ndx} = strtrim(pair_str{ndx});
  end
end
handles.custom_list = custom_list;
trk_file_str = strtrim(get(handles.edit_trk,'String'));
trkfiles = strsplit(trk_file_str,',');
handles.aptStruct.trkfilename = trkfiles;
UpdatePlots(handles);  


function filename = unMacroise(filename,apt)
macros = fieldnames(apt.projMacros); 
for fndx = 1:numel(macros)
  f = macros{fndx};
  filename = strrep(filename,sprintf('$%s',f),apt.projMacros.(f));
end


function trxfile = getTrxFile(trxfile,movfile)
movdir = fileparts(movfile);
trxfile = strrep(trxfile,'$movdir',movdir);


function UpdatePlots(handles)
if isempty(handles.im)
  axis(handles.axes1,'off');
  return
end
axes(handles.axes1);
img = handles.im;
ll = double(handles.locs);

imshow(img);
hold on;
scatter(ll(:,1),ll(:,2),'.');
for ndx = 1:size(ll,1)
  cur_n = modrange(ndx,1,handles.aptStruct.n_pts+1);
  text(ll(ndx,1)+1,ll(ndx,2)+1,num2str(cur_n),'Color',[1 0 0]);
end
for pndx = 1:numel(handles.pairs)
  for view = 1:handles.aptStruct.n_view
    n_pts = handles.aptStruct.n_pts;
    p = ll(handles.pairs{pndx}+(view-1)*n_pts,:);
    plot(p(:,1),p(:,2));
  end
end
hold off;
handles.im = img;
handles.locs = ll;
axis(handles.axes1,'off');


function all_trkfilenames = reset_trk_edit(handles)
apt = handles.apt;
trk_file_str = '';
all_trkfilenames = cell(1,size(apt.movieFilesAll,2));
for ndx = 1:size(apt.movieFilesAll,2)
  %mov_file = apt.movieFilesAll{1,ndx};
  [~,movfilename,~] = fileparts(apt.movieFilesAll{1,ndx});
  [~,projfile,~] = fileparts(handles.aptStruct.lbl_file);
  trkfilename = [movfilename '_' projfile '.trk'];
  if isempty(trk_file_str)
    trk_file_str = trkfilename;
  else
    trk_file_str = sprintf('%s,%s',trk_file_str,trkfilename);
  end
  all_trkfilenames{ndx} = trkfilename;
end
set(handles.edit_trk,'String',trk_file_str)

function isvalid = isvalid_default_trkfilename(handles,trkfilename)
apt = handles.apt;
isvalid = false;
% check number of views
if numel(trkfilename) ~= size(apt.movieFilesAll,2),
  return;
end

[~,projfile,~] = fileparts(handles.aptStruct.lbl_file);
for ndx = 1:size(apt.movieFilesAll,2),
  [~,movfilename,~] = fileparts(apt.movieFilesAll{1,ndx});
  basename = [movfilename '_' projfile];
  if ~startsWith(trkfilename{ndx},basename) || ~endsWith(trkfilename{ndx},'.trk'),
    return;
  end
  
end
isvalid = true;

function onoff = onoffif(in)
if in
  onoff = 'on';
else
  onoff = 'off';
end

% --- Executes on button press in trx_checkbox.
function trx_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to trx_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trx_checkbox


% --- Executes on button press in body_checkbox.
function body_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to body_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of body_checkbox


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox5



function edit_pair_Callback(hObject, eventdata, handles)
% hObject    handle to edit_pair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_pair as text
%        str2double(get(hObject,'String')) returns contents of edit_pair as a double
handles = UpdateHandles(handles);
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function edit_pair_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_pair (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_triad_Callback(hObject, eventdata, handles)
% hObject    handle to edit_triad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_triad as text
%        str2double(get(hObject,'String')) returns contents of edit_triad as a double
handles = UpdateHandles(handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function edit_triad_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_triad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_custom_Callback(hObject, eventdata, handles)
% hObject    handle to edit_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_custom as text
%        str2double(get(hObject,'String')) returns contents of edit_custom as a double
handles = UpdateHandles(handles);
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function edit_custom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in trx_pop.
function trx_pop_Callback(hObject, eventdata, handles)
% hObject    handle to trx_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns trx_pop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from trx_pop
handles = UpdateHandles(handles);
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function trx_pop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to trx_pop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in global_checkbox.
function global_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to global_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of global_checkbox
handles = UpdateHandles(handles);
guidata(hObject,handles)


% --- Executes on button press in theta_checkbox.
function theta_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to theta_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of theta_checkbox
handles = UpdateHandles(handles);
guidata(hObject,handles)

% --- Executes on button press in done_pb.
function done_pb_Callback(hObject, eventdata, handles)
% hObject    handle to done_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = UpdateHandles(handles);
handles.perframe = struct;
if handles.has_trx
  featureLexiconName = handles.aptStruct.origFeatureLexiconName;
  [featureLexicon,animalType]=featureLexiconFromFeatureLexiconName(featureLexiconName);
  if ~handles.use_trx
    f = fieldnames(featureLexicon.perframe);
    for fndx = 1:numel(f)
      featureLexicon.perframe = rmfield(featureLexicon.perframe,f{fndx});
    end
  end
else
  featureLexiconName = 'apt';
  animalType = 'apt';
  [featureLexicon,~]=featureLexiconFromFeatureLexiconName('flies');
  f = fieldnames(featureLexicon.perframe);
  for fndx = 1:numel(f)
    featureLexicon.perframe = rmfield(featureLexicon.perframe,f{fndx});
  end
end

handles.aptStruct.animalType = animalType;
trans_type = {'none' 'relative'}; 
common_struct = struct();
common_struct.trans_types=trans_type;
common_struct.type='apt_locomotion';
dphi_struct = common_struct;
dphi_struct.trans_types = {'none' 'abs'};
dtheta_struct = common_struct;
dtheta_struct.trans_types = {'none'};

n_pts = handles.aptStruct.n_pts;
n_view = handles.aptStruct.n_view;
cur_set = handles.common_ftrs;
for fndx = 1:numel(cur_set)
  for ndx = 1:n_pts
    for view = 1:n_view
      cur_ftr = sprintf('apt_view%d_%s_%d',view,cur_set{fndx},ndx);
      if contains(cur_ftr,'dphi')
        featureLexicon.perframe.(cur_ftr) = dphi_struct;        
      elseif contains(cur_ftr,'dtheta')
        featureLexicon.perframe.(cur_ftr) = dtheta_struct;
      else
        featureLexicon.perframe.(cur_ftr) = common_struct;
      end
    end
  end
end

% if handles.use_global
%   if handles.has_trx 
%     cur_set = handles.global_ftrs;
%   else
%     cur_set = cat(2,handles.global_ftrs,handles.theta_ftrs);
%   end
%   for fndx = 1:numel(cur_set)
%     for ndx = 1:n_pts
%       cur_ftr = sprintf('apt_global_%s_%d',cur_set{fndx},ndx);
%       if contains(cur_ftr,'dphi')
%         featureLexicon.perframe.(cur_ftr) = dphi_struct;        
%       elseif contains(cur_ftr,'dtheta')
%         featureLexicon.perframe.(cur_ftr) = dtheta_struct;
%       else
%         featureLexicon.perframe.(cur_ftr) = common_struct;
%       end
%     end
%   end
% end
% if handles.use_body
%   cur_set = handles.body_ftrs;
%   for fndx = 1:numel(cur_set)
%     for ndx = 1:n_pts
%       cur_ftr = sprintf('apt_body_%s_%d',cur_set{fndx},ndx);
%       if contains(cur_ftr,'dphi')
%         featureLexicon.perframe.(cur_ftr) = dphi_struct;        
%       elseif contains(cur_ftr,'dtheta')
%         featureLexicon.perframe.(cur_ftr) = dtheta_struct;
%       else
%         featureLexicon.perframe.(cur_ftr) = common_struct;
%       end
%     end
%   end
% end
% if handles.use_theta
%   cur_set = handles.theta_ftrs;
%   for fndx = 1:numel(cur_set)
%     for ndx = 1:n_pts
%       cur_ftr = sprintf('apt_body_%s_%d',cur_set{fndx},ndx);
%       if contains(cur_ftr,'dphi')
%         featureLexicon.perframe.(cur_ftr) = dphi_struct;        
%       elseif contains(cur_ftr,'dtheta')
%         featureLexicon.perframe.(cur_ftr) = dtheta_struct;
%       else
%         featureLexicon.perframe.(cur_ftr) = common_struct;
%       end
%     end
%   end
% end

% pair features
for pndx = 1:numel(handles.pairs)
  cur_set = handles.pair_ftrs;
  cur_p = handles.pairs{pndx};
  for fndx = 1:numel(cur_set)
    for view = 1:n_view

      cur_ftr = sprintf('apt_view%d_pair_%s_%d_%d',view,cur_set{fndx},cur_p(1),cur_p(2));
      if contains(cur_ftr,'dphi')
        featureLexicon.perframe.(cur_ftr) = dphi_struct;        
      elseif contains(cur_ftr,'dtheta')
        featureLexicon.perframe.(cur_ftr) = dtheta_struct;
      else
        featureLexicon.perframe.(cur_ftr) = common_struct;
      end
    end
  end
end
% triad features
for pndx = 1:numel(handles.triads)
  cur_set = handles.triad_ftrs;
  cur_p = handles.triads{pndx};
  for fndx = 1:numel(cur_set)
    for view = 1:n_view
      cur_ftr = sprintf('apt_view%d_triad_%s_%d_%d_%d',view,cur_set{fndx},cur_p(1),cur_p(2),cur_p(3));
      if contains(cur_ftr,'dphi')
        featureLexicon.perframe.(cur_ftr) = dphi_struct;        
      elseif contains(cur_ftr,'dtheta')
        featureLexicon.perframe.(cur_ftr) = dtheta_struct;
      else
        featureLexicon.perframe.(cur_ftr) = common_struct;
      end
    end
  end
end

for pndx = 1:numel(handles.custom_list)
  trans_type = {'none' 'relative'}; 
  cur_ftr = handles.custom_list{pndx};
  featureLexicon.perframe.(cur_ftr) = common_struct;
end
handles.aptStruct.featureLexicon = featureLexicon;
handles.aptStruct.has_trx = handles.has_trx;
handles.aptStruct.has_crops = handles.has_crops;
handles.aptStruct.pairs = handles.pairs;
handles.aptStruct.triads = handles.triads;
handles.aptStruct.custom_list = handles.custom_list;

if handles.has_trx
  handles.aptStruct.apt_trx_type = '';
elseif handles.is_ma
  handles.aptStruct.apt_trx_type = 'ma';
  handles.aptStruct.ma_head_tail = handles.head_tail;
  handles.aptStruct.ma_use_theta = handles.use_theta;
else
  if get(handles.rb_frame,'Value')
    handles.aptStruct.apt_trx_type = 'frame';
  elseif  get(handles.rb_crop,'Value')
    handles.aptStruct.apt_trx_type = 'crop';
  elseif  get(handles.rb_custom,'Value')
    handles.aptStruct.apt_trx_type = 'custom';
  else
    handles.aptStruct.apt_trx_type = 'trk';    
    handles.aptStruct.apt_trx_centers = handles.apt_trx_centers;
    handles.aptStruct.apt_trx_orient = handles.apt_trx_orient;    
  end
end
guidata(hObject,handles);
uiresume(handles.figure1);

  
% --- Executes on button press in cancel_pb.
function cancel_pb_Callback(hObject, eventdata, handles)
% hObject    handle to cancel_pb (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.aptStruct.featureLexicon = [];
handles.aptStruct.featureLexiconName = '';
handles.aptStruct.animalType = '';
guidata(hObject,handles);
uiresume(handles.figure1);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
delete(hObject);



function edit_trk_Callback(hObject, eventdata, handles)
% hObject    handle to edit_trk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_trk as text
%        str2double(get(hObject,'String')) returns contents of edit_trk as a double
trk_file_str = strtrim(get(hObject,'String'));
trkfiles = strsplit(trk_file_str,',');
if numel(trkfiles)~=handles.aptStruct.n_view
  uiwait(errormsg('Number of trk files should be same as number of views'));
  handles.aptStruct.trkfilename = reset_trk_edit(handles);
  guidata(hObject,handels);
  return;
end
handles.aptStruct.trkfilename = trkfiles;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function edit_trk_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_trk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rb_crop.
function rb_crop_Callback(hObject, eventdata, handles)
% hObject    handle to rb_crop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_crop


% --- Executes on button press in rb_frame.
function rb_frame_Callback(hObject, eventdata, handles)
% hObject    handle to rb_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_frame


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over rb_centroid.
function rb_centroid_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to rb_centroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in rb_centroid.
function rb_centroid_Callback(hObject, eventdata, handles)
% hObject    handle to rb_centroid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_centroid
in_str = {'Centroid of landmarks to be used as center for JAABA Trx. Eg ''1 2 3''',...
  'Use direction from centroid to this landmark as orientation. Eg ''3''. For no orientation leave it emtpy'};
gg = inputdlg(in_str,'Select landmarks for JAABA Trx');

handles = guidata(hObject);
if isempty(gg)
  reset_rb(handles);
  % reset to default if user cancels.
  return;
end

bb = strsplit(strtrim(gg{1}));
centers = [];
n_pts = handles.aptStruct.n_pts;
for ndx = 1:numel(bb)
  curo = str2double(bb{ndx});
  if isnan(curo) || round(curo) - curo ~= 0 || curo < 1 || curo > n_pts
    uiwait(errordlg(sprintf('Please input only integers or integers less than %d',n_pts)));
    reset_rb(handles)
    return;
  end
  centers(end+1) = curo;
end

if isempty(strtrim(gg{2})),
  orient = 0;
else
  bb = strsplit(strtrim(gg{2}));
  if numel(bb)>1
    uiwait(errordlg('Please input only 1 integer for orientation'));
    reset_rb(handles)
    return;
  end
  
  curo = str2double(bb{1});
  if isnan(curo) || round(curo) - curo ~= 0 || curo < 1 || curo > handles.aptStruct.n_pts
    uiwait(errordlg(sprintf('Please input only integers for orientation or integers less than %d',n_pts)));
    reset_rb(handles)
    return;
  end
  orient = curo;
end

handles.apt_trx_centers = centers;
handles.apt_trx_orient = orient;
handles.apt_trx_type = 'trk';
tt = sprintf('%d ', centers);
cstr = sprintf('Centroid: %s',tt);
tstr = sprintf('Orientation: %d',orient);
set(handles.text_centroid,'String',{cstr,tstr});
guidata(hObject,handles);

function reset_rb(handles)
  if handles.has_crops
    set(handles.rb_crop,'Value',true);
    set(handles.rb_frame,'Value',false);
  else
    set(handles.rb_crop,'Value',false);
    set(handles.rb_frame,'Value',true);
  end
  set(handles.rb_custom,'Value',false);
  set(handles.rb_centroid,'Value',false);
  set(handles.text_centroid,'String','');

% --- Executes on button press in rb_custom.
function rb_custom_Callback(hObject, eventdata, handles)
% hObject    handle to rb_custom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_custom


% --- Executes on button press in pushbutton_help.
function pushbutton_help_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_help (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
web('https://github.com/kristinbranson/APT/wiki/Using-APT-tracking-data-in-JAABA','-browser');
