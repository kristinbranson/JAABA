function [success,msg] = ConvertMotrCtc2JAABA(varargin)

% For converting Motr/Catalytic .ctc track files to JAABA format.

success = false;
msg = {};

%% parse inputs

[inmoviefile,seqindexfile,...
  intrxfile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  pxpermm,...
  dofliplr,doflipud,dotransposeimage,...
  dosoftlink,frameinterval,...
  sex] = myparse(varargin,...
  'inmoviefile','','seqindexfile','',...
  'intrxfile','',...
  'expdir','','moviefilestr','movie.ufmf','trxfilestr','trx.mat','perframedirstr','perframe',...
  'arenatype','None','arenacenterx',0,'arenacentery',0,...
  'arenaradius',123,'arenawidth',123,'arenaheight',123,...
  'pxpermm',1,...
  'fliplr',false,'flipud',false,'dotransposeimage',false,...
  'dosoftlink',false,...
  'frameinterval',[],...
  'sex',{});

%% check that required inputs are given
if isempty(inmoviefile),
  msg = 'Input movie file is empty';
  return;
end
[~,~,ext] = fileparts(inmoviefile);
% if ~strcmpi(ext,'.seq'),
%   msg = 'Input movie file must be a .seq file';
%   return;
% end
isseqmovie = strcmp(ext,'.seq');
if isseqmovie && isempty(seqindexfile),
  msg = 'Input seq index file is empty';
  return;
end
if isempty(intrxfile),
  msg = 'Input trx mat file is empty';
  return;
end

%% output file locations
moviefile = fullfile(expdir,moviefilestr);
[~,name] = fileparts(moviefilestr);
seqindexfilestr = [name,'.mat'];
outseqindexfile = fullfile(expdir,seqindexfilestr);
trxfile = fullfile(expdir,trxfilestr);
perframedir = fullfile(expdir,perframedirstr);

%% load in the data
td = load('-mat',intrxfile);  % .ctc file
trx=td.trx;
clear td;  % no further need for this

% Make sure trx is a struct
if ~isstruct(trx) ,
    msg='Tracks file is not a valid .ctc file';
    return
end

% If trx has a pxpermm field, use that instead of the pixpermm from the
% PrepareJAABAData UI, so user doesn't have to re-enter it.
if ~isempty(trx) && isfield(trx,'pxpermm') ,
    pxpermm=trx(1).pxpermm;
end

% remove all but the 'core' fields
coreFieldNames={'x' 'y' 'a' 'b' 'theta'};
allFieldNames=fieldnames(trx);
fieldNamesToRemove=setdiff(allFieldNames,coreFieldNames);
trx=rmfield(trx,fieldNamesToRemove);

% % orientation is backwards
% for i = 1:numel(td.astrctTrackers),
%   td.astrctTrackers(i).m_afTheta = -td.astrctTrackers(i).m_afTheta;
% end
n_mice = numel(trx);
n_frames = numel(trx(1).x);
       
if isseqmovie
  headerinfo = r_readseqinfo(inmoviefile);
  vidmetadata.imwidth = headerinfo.m_iWidth;
  vidmetadata.imheight = headerinfo.m_iHeight;
  vidmetadata.timestamps = headerinfo.m_afTimestamp;
else
  try
    vrdr = VideoReader(inmoviefile);
    vidmetadata.imwidth = vrdr.Width;
    vidmetadata.imheight = vrdr.Height;
    vidmetadata.timestamps = [];
    vidmetadata.framerate = vrdr.FrameRate;
  catch ME
    msg = sprintf('Cannot read metadata from video file: %s',ME.message);
    return;
  end
end

%% if sex not entered for each mice, ask for it

if numel(sex) <= n_mice,
  prompt = cell(1,n_mice);
  for i = 1:n_mice,
    prompt{i} = sprintf('Mouse %d',i);
  end
  defans = sex;
  for i = numel(sex)+1:n_mice,
    if i > n_mice/2,
      defans{i} = 'M';
    else
      defans{i} = 'F';
    end
  end
  while true,
    res = inputdlg(prompt,'Enter sex for each mouse. M = male, F = female.',1,defans,'on');
    if isempty(res),
      msg = 'Sex of each mouse not entered.';
      return;
    end
    res = upper(res);
    isok = ~cellfun(@isempty,regexp(res,'^(M|F)$','once'));
    if all(isok),
      sex = res;
      break;
    end
    uiwait(warndlg('Sex must be M or F for each mouse'));
    if any(isok),
      [defans{isok}] = deal(res{isok});
    end
  end
end

%% flip if necessary
if dofliplr,
  for i = 1:n_mice,
    trx(i).x = vidmetadata.imwidth-trx(i).x+1;
    trx(i).theta = modrange(pi - trx(i).theta,-pi,pi);
  end
  msg{end+1} = 'Flipped the trajectories left-right.';
end
if doflipud,
  for i = 1:n_mice,
    trx(i).y = vidmetadata.imheight-trx(i).y+1;
    trx(i).theta = modrange(-trx(i).theta,-pi,pi);
  end
  msg{end+1} = 'Flipped the trajectories up-down.';
end
if dotransposeimage,
  for i = 1:numel(tx),
    tmp = trx(i).x;
    trx(i).x = trx(i).y;
    trx(i).y = tmp;
    c = cos(trx(i).theta);
    s = sin(trx(i).theta);
    trx(i).theta = atan2(c,s);
  end
  msg{end+1} = 'Switched x and y in trajectories.';
end

%%

% crop if desired
if ~isempty(frameinterval),
  frameinterval(1) = max(frameinterval(1),1);
  frameinterval(2) = min(frameinterval(2),numel(trx(1).x));
  fns = fieldnames(trx);
  for i = 1:numel(fns),
    fn = fns{i};
    for j = 1:numel(trx),
      trx(j).(fn) = trx(j).(fn)(frameinterval(1):frameinterval(2));
    end
  end
  if ~isempty(vidmetadata.timestamps)
    timestamps = vidmetadata.timestamps(frameinterval(1):frameinterval(2));
  else
    timestamps = vidmetadata.timestamps;
  end
else
  timestamps = vidmetadata.timestamps;
end

% fps
if isempty(timestamps)
  fps = vidmetadata.framerate;
  dt = 1/fps;
else
  dt = diff(timestamps);
  fps = 1/median(dt); % computed from *cropped* timestamps
  msg{end+1} = sprintf('Computed frame rate = %f fps from timestamps.',fps);
end

if strcmpi(arenatype,'None'),
  arenacenterx = 0;
  arenacentery = 0;
end

%% create new trx

trx = struct('x',{trx.x},...
             'y',{trx.y},...
             'theta',{trx.theta},...
             'a',{trx.a},...
             'b',{trx.b},...
             'firstframe',num2cell(ones(1,n_mice)),...
             'arena',cell(1,n_mice),...
             'off',num2cell(zeros(1,n_mice)),...
             'nframes',num2cell(n_frames(ones(1,n_mice))),...
             'endframe',num2cell(n_frames(ones(1,n_mice))),...
             'timestamps',repmat({timestamps},[1,n_mice]),...
             'moviename',repmat({inmoviefile},[1,n_mice]),...
             'annname',repmat({intrxfile},[1,n_mice]),...
             'matname',repmat({trxfile},[1,n_mice]),...
             'x_mm',cellfun(@(x) (x-arenacenterx)/pxpermm,{trx.x},'UniformOutput',false),...
             'y_mm',cellfun(@(x) (x-arenacentery)/pxpermm,{trx.y},'UniformOutput',false),...
             'a_mm',cellfun(@(x) x/pxpermm,{trx.a},'UniformOutput',false),...
             'b_mm',cellfun(@(x) x/pxpermm,{trx.b},'UniformOutput',false),...
             'theta_mm',{trx.theta},...
             'dt',repmat({dt},[1,n_mice]),...
             'fps',repmat({fps},[1,n_mice]),...
             'pxpermm',repmat({pxpermm},[1,n_mice]));
for i = 1:numel(sex),
  trx(i).sex = sex{i};
end

msg{end+1} = sprintf('Read trx for %d mice, frame range [%d,%d]',n_mice,min([trx.firstframe]),max([trx.endframe]));
msg{end+1} = sprintf('Using input pixels-per-mm scaling = %f.',pxpermm);

arenaradius_mm = arenaradius / pxpermm;
arenawidth_mm = arenawidth / pxpermm;
arenaheight_mm = arenaheight / pxpermm;
arenacenterx_mm = 0;
arenacentery_mm = 0;
trx = SetLandmarkParameters(trx,arenatype,arenacenterx_mm,arenacentery_mm,...
  arenaradius_mm,arenawidth_mm,arenaheight_mm);

if strcmpi(arenatype,'Circle')
  msg{end+1} = sprintf('Set circular arena parameters. Center = (0,0) mm, radius = %.1f mm.',arenaradius_mm);
elseif strcmpi(arenatype,'Rectangle'),
  msg{end+1} = sprintf('Set rectangular arena parameters. Center = (0,0) mm, width= %.1f mm, height = %.1f mm',arenawidth_mm,arenaheight_mm);  
end
msg{end+1} = sprintf('x_mm ranges within [%.1f,%.1f], y_mm ranges within [%.1f,%.1f]',...
  min([trx.x_mm]),max([trx.x_mm]),min([trx.y_mm]),max([trx.y_mm]));

%% create the experiment directory
if ~exist(expdir,'dir'),
  [success1,msg1] = mkdir(expdir);
  if ~success1,
    msg = msg1;
    return;
  end
  msg{end+1} = sprintf('Created experiment directory %s.',expdir);
end

%% save the trx file
try
  save(trxfile,'trx','timestamps');
catch ME,
  msg = sprintf('Could not save to file %s: %s',trxfile,getReport(ME));
  return;
end
if ~exist(trxfile,'file'),
  msg = sprintf('Failed to save trx to file %s',trxfile);
  return;
end
msg{end+1} = sprintf('Saved trx to file %s.',trxfile);

%% copy/soft-link movie

if strcmp(fullfile(inmoviefile),fullfile(moviefile)),
  msg{end+1} = 'Input and out movie files are the same, not copying/linking.';
else
  
  
  if dosoftlink,
    if exist(moviefile,'file'),
      delete(moviefile);
    end
    if isunix,
      cmd = sprintf('ln -s ''%s'' ''%s''',inmoviefile,moviefile);
      unix(cmd);
      % test to make sure it worked
      [status,result] = unix(sprintf('readlink %s',moviefile));
      result = strtrim(result);
      if status ~= 0 || ~strcmp(result,inmoviefile),
        res = questdlg(sprintf('Failed to make soft link. Copy %s to %s instead?',inmoviefile,moviefile));
        if ~strcmpi(res,'Yes'),
          msg = sprintf('Failed to make soft link from %s to %s.',inmoviefile,moviefile);
          return;
        end
        dosoftlink = false;
      end
    elseif ispc,
      if exist([moviefile,'.lnk'],'file'),
        delete([moviefile,'.lnk']);
      end
      cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',inmoviefile,moviefile);
      fprintf('Making a Windows shortcut file at "%s" with target "%s"\n',inmoviefile,moviefile);
      system(cmd);
      % test to make sure that worked
      [equalmoviefile,didfind] = GetPCShortcutFileActualPath(moviefile);
      if ~didfind || ~strcmp(equalmoviefile,inmoviefile),
        res = questdlg(sprintf('Failed to make shortcut. Copy %s to %s instead?',inmoviefile,moviefile));
        if ~strcmpi(res,'Yes'),
          msg = sprintf('Failed to make shortcut from %s to %s.',inmoviefile,moviefile);
          return;
        end
        dosoftlink = false;
      end
    else
      res = questdlg(sprintf('Unknown OS, cannot soft-link movie file %s. Copy instead?',inmoviefile));
      if ~strcmpi(res,'Yes'),
        msg = sprintf('Failed to make softlink from %s to %s.',inmoviefile,moviefile);
        return;
      end
      dosoftlink = false;
    end
    if dosoftlink,
      msg{end+1} = sprintf('Made a link to movie file %s at %s',inmoviefile,moviefile);
    end
  end
  
  if ~dosoftlink,
    if ispc,
      if exist([moviefile,'.lnk'],'file'),
        delete([moviefile,'.lnk']);
      end
    end
    if exist(moviefile,'file'),
      delete(moviefile);
    end
    [success1,msg1] = copyfile(inmoviefile,moviefile);
    if ~success1,
      msg = msg1;
      success = false;
      return;
    end
    msg{end+1} = sprintf('Copied movie file %s to %s',inmoviefile,moviefile);
  end
  
end

%% copy/soft-link index file

if isseqmovie
  if isempty(frameinterval),
    
    if strcmp(fullfile(seqindexfile),fullfile(outseqindexfile)),
      msg{end+1} = 'Input and out seq index files are the same, not copying/linking.';
    else
      
      if dosoftlink,
        if ispc,
          cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',seqindexfile,outseqindexfile);
          fprintf('Making a Windows shortcut file at "%s" with target "%s"\n',outseqindexfile,seqindexfile);
        else
          cmd = sprintf('ln -s "%s" "%s"',seqindexfile,outseqindexfile);
          fprintf('Soft-linking from "%s" with target "%s"\n',outseqindexfile,seqindexfile);
        end
        system(cmd);
        
      else
        
        fprintf('Copying "%s" to "%s"\n',seqindexfile,outseqindexfile);
        copyfile(seqindexfile,outseqindexfile);
        
      end
      
    end
    
  else
    
    if strcmp(fullfile(seqindexfile),fullfile(outseqindexfile)),
      msg = 'Input and output seq index file names are the same, but trying to crop seq file. Data will be lost, so this is not allowed.';
      return;
    end
    
    % load in index file
    indexdata = load(seqindexfile);
    
    % crop
    indexdata.aiSeekPos = indexdata.aiSeekPos(frameinterval(1):frameinterval(2));
    indexdata.afTimestamp = indexdata.afTimestamp(frameinterval(1):frameinterval(2));
    indexdata.frameinterval = frameinterval;
    
    % save
    save(outseqindexfile,'-struct','indexdata');
  end
end

%% make per-frame directory
if ~exist(perframedir,'dir'),
  [success1,msg1] = mkdir(perframedir);
  if ~success1,
    msg = msg1;
    return;
  end
end

%% done

success = true;

end  % function
