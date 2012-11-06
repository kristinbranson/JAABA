function [success,msg] = ConvertMoTr2JAABA(varargin)

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
if ~strcmpi(ext,'.seq'),
  msg = 'Input movie file must be a .seq file';
  return;
end
if isempty(seqindexfile),
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

headerinfo = r_readseqinfo(inmoviefile);
td = load(intrxfile);
% orientation is backwards
for i = 1:numel(td.astrctTrackers),
  td.astrctTrackers(i).m_afTheta = -td.astrctTrackers(i).m_afTheta;
end
% get image size if flipping
imheight = headerinfo.m_iHeight;
imwidth = headerinfo.m_iWidth;
nmice = numel(td.astrctTrackers);

%% if sex not entered for each mice, ask for it

if numel(sex) <= nmice,
  prompt = cell(1,nmice);
  for i = 1:nmice,
    prompt{i} = sprintf('Mouse %d',i);
  end
  defans = sex;
  for i = numel(sex)+1:nmice,
    if i > nmice/2,
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
  for i = 1:numel(td.astrctTrackers),
    td.astrctTrackers(i).m_afX = imwidth-td.astrctTrackers(i).m_afX+1;
    td.astrctTrackers(i).m_afTheta = modrange(pi - td.astrctTrackers(i).m_afTheta,-pi,pi);
  end
  msg{end+1} = 'Flipped the trajectories left-right.';
end
if doflipud,
  for i = 1:numel(td.astrctTrackers),
    td.astrctTrackers(i).m_afY = imheight-td.astrctTrackers(i).m_afY+1;
    td.astrctTrackers(i).m_afTheta = modrange(-td.astrctTrackers(i).m_afTheta,-pi,pi);
  end
  msg{end+1} = 'Flipped the trajectories up-down.';
end
if dotransposeimage,
  for i = 1:numel(tx),
    tmp = td.astrctTrackers(i).m_afX;
    td.astrctTrackers(i).m_afX = td.astrctTrackers(i).m_afY;
    td.astrctTrackers(i).m_afY = tmp;
    c = cos(td.astrctTrackers(i).m_afTheta);
    s = sin(td.astrctTrackers(i).m_afTheta);
    td.astrctTrackers(i).m_afTheta = atan2(c,s);
  end
  msg{end+1} = 'Switched x and y in trajectories.';
end

%%

% crop if desired
if ~isempty(frameinterval),
  frameinterval(1) = max(frameinterval(1),1);
  frameinterval(2) = min(frameinterval(2),numel(td.astrctTrackers(1).m_afX));
  fns = fieldnames(td.astrctTrackers);
  for i = 1:numel(fns),
    fn = fns{i};
    for j = 1:numel(td.astrctTrackers),
      td.astrctTrackers(j).(fn) = td.astrctTrackers(j).(fn)(frameinterval(1):frameinterval(2));
    end
  end
  timestamps = headerinfo.m_afTimestamp(frameinterval(1):frameinterval(2));
else
  timestamps = headerinfo.m_afTimestamp;
end

% count
nframes = numel(td.astrctTrackers(1).m_afX);

% fps
dt = diff(timestamps);
fps = 1/median(dt);
msg{end+1} = sprintf('Computed frame rate = %f fps from timestamps.',fps);

if strcmpi(arenatype,'None'),
  arenacenterx = 0;
  arenacentery = 0;
end

%% create new trx

trx = struct('x',{td.astrctTrackers.m_afX},...
  'y',{td.astrctTrackers.m_afY},...
  'theta',{td.astrctTrackers.m_afTheta},...
  'a',cellfun(@(x) x/2,{td.astrctTrackers.m_afA},'UniformOutput',false),...
  'b',cellfun(@(x) x/2,{td.astrctTrackers.m_afB},'UniformOutput',false),...
  'firstframe',num2cell(ones(1,nmice)),...
  'arena',cell(1,nmice),...
  'off',num2cell(zeros(1,nmice)),...
  'nframes',num2cell(nframes(ones(1,nmice))),...
  'endframe',num2cell(nframes(ones(1,nmice))),...
  'timestamps',repmat({timestamps},[1,nmice]),...
  'moviename',repmat({inmoviefile},[1,nmice]),...
  'annname',repmat({intrxfile},[1,nmice]),...
  'matname',repmat({trxfile},[1,nmice]),...
  'x_mm',cellfun(@(x) (x-arenacenterx)/pxpermm,{td.astrctTrackers.m_afX},'UniformOutput',false),...
  'y_mm',cellfun(@(x) (x-arenacentery)/pxpermm,{td.astrctTrackers.m_afY},'UniformOutput',false),...
  'a_mm',cellfun(@(x) x/pxpermm/2,{td.astrctTrackers.m_afA},'UniformOutput',false),...
  'b_mm',cellfun(@(x) x/pxpermm/2,{td.astrctTrackers.m_afB},'UniformOutput',false),...
  'theta_mm',cellfun(@(x) -x, {td.astrctTrackers.m_afTheta},'UniformOutput',false),...
  'dt',repmat({dt},[1,nmice]),...
  'fps',repmat({fps},[1,nmice]),...
  'pxpermm',repmat({pxpermm},[1,nmice]));
for i = 1:numel(sex),
  trx(i).sex = sex{i};
end

msg{end+1} = sprintf('Read trx for %d mice, frame range [%d,%d]',nmice,min([trx.firstframe]),max([trx.endframe]));
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

if dosoftlink,
  if exist(moviefile,'file'),
    delete(moviefile);
  end
  if isunix,
    cmd = sprintf('ln -s %s %s',inmoviefile,moviefile);
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

%% copy/soft-link index file

if isempty(frameinterval),
  
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
  
else
  
  % load in index file
  indexdata = load(seqindexfile);
  
  % crop
  indexdata.aiSeekPos = indexdata.aiSeekPos(frameinterval(1):frameinterval(2));
  indexdata.afTimestamp = indexdata.afTimestamp(frameinterval(1):frameinterval(2));
  indexdata.frameinterval = frameinterval;
  
  % save
  save(outseqindexfile,'-struct','indexdata');
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