function [success,msg] = ConvertMouseHouse2JAABA(varargin)

success = false;
msg = '';

%% parse inputs

[inmoviefile,seqindexfile,...
  intrxfile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  pxpermm,fps,overridefps,...
  dosoftlink,frameinterval,...
  sex] = myparse(varargin,...
  'inmoviefile','','seqindexfile','',...
  'intrxfile','',...
  'expdir','','moviefilestr','movie.ufmf','trxfilestr','trx.mat','perframedirstr','perframe',...
  'arenatype','None','arenacenterx',0,'arenacentery',0,...
  'arenaradius',123,'arenawidth',123,'arenaheight',123,...
  'pxpermm',1,'fps',30,...
  'overridefps',false,...
  'dosoftlink',false,...
  'frameinterval',[],...
  'sex',{});

%% check that required inputs are given
if isempty(inmoviefile),
  msg = 'Input movie file is empty';
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
tmp = load(intrxfile);

% crop if desired
if ~isempty(frameinterval),
  fns = fieldnames(tmp.astrctTrackers);
  for i = 1:numel(fns),
    fn = fns{i};
    for j = 1:numel(tmp.astrctTrackers),
      tmp.astrctTrackers(j).(fn) = tmp.astrctTrackers(j).(fn)(frameinterval(1):frameinterval(2));
    end
  end
  timestamps = headerinfo.m_afTimestamp(frameinterval(1):frameinterval(2));
else
  timestamps = headerinfo.m_afTimestamp;
end

% count
nmice = numel(tmp.astrctTrackers);
nframes = numel(tmp.astrctTrackers(1).m_afX);

% fps
if overridefps,
  timestamps = timestamps(1)+(0:nframes-1)/fps;
end
dt = diff(timestamps);
if ~overridefps,
  fps = 1/median(dt);
end

if strcmpi(arenatype,'None'),
  arenacenterx = 0;
  arenacentery = 0;
end

%% create new trx

trx = struct('x',{tmp.astrctTrackers.m_afX},...
  'y',{tmp.astrctTrackers.m_afY},...
  'theta',cellfun(@(x) -x, {tmp.astrctTrackers.m_afTheta},'UniformOutput',false),...
  'a',cellfun(@(x) x/2,{tmp.astrctTrackers.m_afA},'UniformOutput',false),...
  'b',cellfun(@(x) x/2,{tmp.astrctTrackers.m_afB},'UniformOutput',false),...
  'firstframe',num2cell(ones(1,nmice)),...
  'arena',cell(1,nmice),...
  'off',num2cell(zeros(1,nmice)),...
  'nframes',num2cell(nframes(ones(1,nmice))),...
  'endframe',num2cell(nframes(ones(1,nmice))),...
  'timestamps',repmat({timestamps},[1,nmice]),...
  'moviename',repmat({inmoviefile},[1,nmice]),...
  'annname',repmat({intrxfile},[1,nmice]),...
  'matname',repmat({trxfile},[1,nmice]),...
  'x_mm',cellfun(@(x) (x-arenacenterx)/pxpermm,{tmp.astrctTrackers.m_afX},'UniformOutput',false),...
  'y_mm',cellfun(@(x) (x-arenacentery)/pxpermm,{tmp.astrctTrackers.m_afY},'UniformOutput',false),...
  'a_mm',cellfun(@(x) x/pxpermm/2,{tmp.astrctTrackers.m_afA},'UniformOutput',false),...
  'b_mm',cellfun(@(x) x/pxpermm/2,{tmp.astrctTrackers.m_afB},'UniformOutput',false),...
  'theta_mm',cellfun(@(x) -x, {tmp.astrctTrackers.m_afTheta},'UniformOutput',false),...
  'dt',repmat({dt},[1,nmice]),...
  'fps',repmat({fps},[1,nmice]),...
  'pxpermm',repmat({pxpermm},[1,nmice]));
for i = 1:numel(sex),
  trx(i).sex = sex{i};
end

arenaradius_mm = arenaradius / pxpermm;
arenawidth_mm = arenawidth / pxpermm;
arenaheight_mm = arenaheight / pxpermm;
arenacenterx_mm = 0;
arenacentery_mm = 0;
trx = SetLandmarkParameters(trx,arenatype,arenacenterx_mm,arenacentery_mm,...
  arenaradius_mm,arenawidth_mm,arenaheight_mm); %#ok<NASGU>


%% create the experiment directory
if ~exist(expdir,'dir'),
  [success1,msg1] = mkdir(expdir);
  if ~success1,
    msg = msg1;
    return;
  end
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

%% copy/soft-link movie
if dosoftlink,
  if isunix,
    cmd = sprintf('ln -s %s %s',inmoviefile,moviefile);
    unix(cmd);
    % test to make sure it worked
    [status,result] = unix(sprintf('readlink %s',moviefile));
    result = strtrim(result);
    if status ~= 0 || ~strcmp(result,inmoviefile),
      warndlg(sprintf('Failed to make soft link, copying %s to %s instead',inmoviefile,moviefile));
      dosoftlink = false;
    end
  elseif ispc,
    cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',inmoviefile,moviefile);
    fprintf('Making a Windows shortcut file at "%s" with target "%s"\n',inmoviefile,moviefile);
    system(cmd);
    % test to make sure that worked
    [equalmoviefile,didfind] = GetPCShortcutFileActualPath(moviefile);
    if ~didfind || ~strcmp(equalmoviefile,inmoviefile),
      res = questdlg(sprintf('Failed to make shortcut, copy %s to %s instead?',inmoviefile,moviefile));
      if ismember(res,{'No','Cancel'}),
        msg = sprintf('Failed to make shortcut from file %s to %s',inmoviefile,moviefile);
        return;
      end
      dosoftlink = false;
    end
  else
    warndlg(sprintf('Unknown OS, not soft-linking movie file %s',inmoviefile));
    dosoftlink = false;
  end  
end
  
if ~dosoftlink,
  [success1,msg1] = copyfile(inmoviefile,moviefile);
  if ~success1,
    msg = msg1;
    success = false;
    return;
  end
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