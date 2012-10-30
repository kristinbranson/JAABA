function [success,msg] = ConvertQtrax2JAABA(varargin)

success = false;
msg = '';

[inmoviefile,featfile,roifile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  fps,overridefps,...
  dosoftlink,...
  dosavecadabrafeats] = myparse(varargin,...
  'inmoviefile','','featfile','','roifile','',...
  'expdir','','moviefilestr','movie.ufmf','trxfilestr','trx.mat','perframedirstr','perframe',...
  'arenatype','None','arenacenterx',0,'arenacentery',0,...
  'arenaradius',123,'arenawidth',123,'arenaheight',123,...
  'fps',30,...
  'overridefps',false,...
  'dosoftlink',false,...
  'dosavecadabrafeats',false);

% check that required inputs are given
if isempty(inmoviefile),
  msg = 'Input movie file is empty';
  return;
end
if isempty(featfile),
  msg = 'Input trx mat file is empty';
  return;
end
if isempty(roifile),
  msg = 'Input roi file is empty';
  return;
end

% output file locations
moviefile = fullfile(expdir,moviefilestr);
trxfile = fullfile(expdir,trxfilestr);
perframedir = fullfile(expdir,perframedirstr);

%% load in data

feat = load(featfile);
roi = load(roifile);
if doflipud || dofliplr,

  % read frame size, moviename is either the name of the movie or an mmreader object
  [readframe,~,fid] = get_readframe_fcn(moviefile);
  im = readframe(1);
  movieheight = size(im,1);
  moviewidth = size(im,2);
  if fid > 1,
    fclose(fid);
  end
end

obj = [feat.fly_feat.obj1,feat.fly_feat.obj2];

%% convert

% scale, pxpermm, hopefully these are the same
pxpermm = 1/mean([roi.scale.x,roi.scale.y]);

% get timestamps
timestamps = feat.fly_feat.time;
if overridefps,
  timestamps = timestamps(1) + (0:numel(feat.fly_feat.time)-1)/fps;
else
  fps = 1/median(diff(timestamps));
end

% allocate
c = cell(1,2);
trx = struct('x',c,'y',c,'theta',c,'a',c,'b',c,...
  'id',c,'moviename',c,'firstframe',c,'arena',c,...
  'nframes',c,'endframe',c,'pxpermm',c,'fps',c,'x_mm',c,...
  'y_mm',c,'a_mm',c,'b_mm',c,'theta_mm',c,'dt',c);

firstframeoff = feat.fly_feat.frame(1) - 1;
for fly = 1:2,
  
  % frames for which fly is tracked
  
  % (x,y) = (0,0) for untracked frames
  badframes = obj(fly).pos_x == 0 & obj(fly).pos_y == 0;
  lastframe = find(~badframes,1,'last');
  firstframe = find(~badframes,1);
  if isempty(lastframe), 
    firstframe = 1;
    lastframe = 0;
  end
  nframes = lastframe - firstframe + 1;
  
  % allocate
  trx(fly).x = nan(1,nframes);
  trx(fly).y = nan(1,nframes);
  trx(fly).theta = nan(1,nframes);
  trx(fly).a = nan(1,nframes);
  trx(fly).b = nan(1,nframes);
  trx(fly).xwingl = nan(1,nframes);
  trx(fly).ywingl = nan(1,nframes);
  trx(fly).xwingr = nan(1,nframes);
  trx(fly).ywingr = nan(1,nframes);
  trx(fly).x_mm = nan(1,nframes);
  trx(fly).y_mm = nan(1,nframes);
  trx(fly).a_mm = nan(1,nframes);
  trx(fly).b_mm = nan(1,nframes);
  trx(fly).theta_mm = nan(1,nframes);
  trx(fly).dt = nan(1,nframes-1);
  
  % store parameters
  trx(fly).moviename = moviefile;
  trx(fly).firstframe = firstframe + firstframeoff;
  trx(fly).endframe = lastframe + firstframeoff;
  trx(fly).nframes = nframes;
  trx(fly).pxpermm = pxpermm;
  trx(fly).fps = fps;
  trx(fly).id = fly;
  
  % store data
  idx = feat.fly_feat.frame(firstframe:lastframe)-trx(fly).firstframe + 1;

  trx(fly).x_mm(idx) = obj(fly).pos_x;
  trx(fly).y_mm(idx) = obj(fly).pos_y;

  % we use the quarter major, minor axis length
  trx(fly).a_mm(idx) = obj(fly).FLength/4;
  % maybe area is major/2 * minor/2 * pi, store quarter minor axis length
  trx(fly).b_mm(idx) = (obj(fly).FArea./(obj(fly).FLength/2)/pi)/2;
  % convert degrees to radians
  trx(fly).theta(idx) = obj(fly).headdir*pi/180;
  trx(fly).theta_mm = trx(fly).theta;
  
  % convert mm to px, incorporate offset
  trx(fly).x(idx) = obj(fly).pos_x*roi.scale.x + roi.ROI.cols(1) - 1;
  trx(fly).y(idx) = obj(fly).pos_y*roi.scale.y + roi.ROI.rows(1) - 1 - 1;
  trx(fly).a = trx(fly).a_mm*pxpermm;
  trx(fly).b = trx(fly).b_mm*pxpermm;
    
  phil = obj(fly).phil*pi/180;
  wingll = obj(fly).wingll*pxpermm;
  xwl = trx(fly).x(idx) + wingll.*cos(-phil-trx(fly).theta(idx)-pi);
  ywl = trx(fly).y(idx) + wingll.*sin(-phil-trx(fly).theta(idx));
  phir = obj(fly).phir*pi/180;
  winglr = obj(fly).winglr*pxpermm;
  xwr = trx(fly).x(idx) + winglr.*cos(phir-trx(fly).theta(idx)-pi);
  ywr = trx(fly).y(idx) + winglr.*sin(phir-trx(fly).theta(idx));
  trx(fly).xwingl(idx) = xwl;
  trx(fly).ywingl(idx) = ywl;
  trx(fly).xwingr(idx) = xwr;
  trx(fly).ywingr(idx) = ywr;  
  
  trx(fly).dt(idx(1:end-1)) = diff(timestamps(firstframe:lastframe));
  
  % flipud if necessary
  if doflipud,
    trx(fly).y = movieheight - trx(fly).y;
    trx(fly).ywingl = movieheight - trx(fly).ywingl;
    trx(fly).ywingr = movieheight - trx(fly).ywingr;
    trx(fly).theta = -trx(fly).theta;
  end

  % fliplr if necessary
  if dofliplr,
    trx(fly).x = moviewidth - trx(fly).x;
    trx(fly).xwingl = moviewidth - trx(fly).xwingl;
    trx(fly).xwingr = moviewidth - trx(fly).xwingr;
    trx(fly).theta = modrange(pi-trx(fly).theta,-pi,pi);
  end
  
end

% set landmark parameters
arenaradius_mm = arenaradius / pxpermm;
arenawidth_mm = arenawidth / pxpermm;
arenaheight_mm = arenaheight / pxpermm;
arenacenterx_mm = (arenacenterx - roi.ROI.cols(1) + 1)/roi.ROI.scale.x;
arenacentery_mm = (arenacentery - roi.ROI.rows(1) + 1 + 1)/roi.ROI.scale.y;
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
      warndlg(sprintf('Failed to make shortcut, copying %s to %s instead',inmoviefile,moviefile));
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

%% make per-frame directory

if ~exist(perframedir,'dir'),
  [success1,msg1] = mkdir(perframedir);
  if ~success1,
    msg = msg1;
    return;
  end
end

%% save the CADABRA per-frame features
% TODO: figure out what these variables are and save them in a more
% reasonable way

if dosavecadabrafeats,
  
  % allocate variables that will be saved
  data = cell(1,2);
  units = parseunits('unit'); %#ok<NASGU>
  
  % frames for which fly is tracked
  firstframes = nan(1,2);
  lastframes = nan(1,2);
  for fly = 1:2,
    % (x,y) = (0,0) for untracked frames
    badframes = obj(fly).pos_x == 0 & obj(fly).pos_y == 0;
    lastframe = find(~badframes,1,'last');
    firstframe = find(~badframes,1);
    if isempty(lastframe),
      firstframe = 1;
      lastframe = 0;
    end
    firstframes(fly) = firstframe;
    lastframes(fly) = lastframe;
  end
  
  % joint features
  fns_ignore = {'ind1_count'
    'ind2_count'
    'obj1'
    'obj2'
    'frame'
    'time'
    };
  fns = setdiff(fieldnames(feat.fly_feat),fns_ignore);
  nframesall = numel(feat.fly_feat.frame);
  for i = 1:numel(fns),
    fn = fns{i};
    n = numel(feat.fly_feat.(fn));
    offend = ceil((nframesall-n)/2);
    offstart = nframesall-n-offend;
    for fly = 1:2,
      data{fly} = feat.fly_feat.(fn)(firstframes(fly)+offstart:lastframes(fly)-offend);
    end
    save(fullfile(perframedir,[fn,'.mat']),'data','units');
  end
  
  % per-fly features
  fns = fieldnames(obj);
  for i = 1:numel(fns),
    fn = fns{i};
    for fly = 1:2,
      data{fly} = obj(fly).(fn);
    end
    save(fullfile(perframedir,[fn,'.mat']),'data','units');
  end
  
end

%% done

success = true;