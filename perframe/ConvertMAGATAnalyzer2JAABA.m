function [success,msg] = ConvertMAGATAnalyzer2JAABA(varargin)

SCALE = 10; % centimeters instead of millimeters
MAXTIMESTAMPERR = .01;

success = false;
msg = '';

%% set parameters

[inmoviefile,expfile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenatype,arenacenterx,arenacentery,...
  arenaradius,arenawidth,arenaheight,...
  dosoftlink,dotransposeimage] = myparse(varargin,...
  'inmoviefile','','expfile','',...
  'expdir','','moviefilestr','movie.ufmf','trxfilestr','trx.mat','perframedirstr','perframe',...
  'arenatype','None','arenacenterx',0,'arenacentery',0,...
  'arenaradius',123,'arenawidth',123,'arenaheight',123,...
  'dosoftlink',false,...
  'dotransposeimage',false);

%% check for files

if isempty(expfile),
  msg = 'Experiment file is empty';
  return;
end

if isempty(inmoviefile),
  msg = 'Video file is empty';
  return;
end

if ~exist(expfile,'file'),
  msg = sprintf('Experiment file %s does not exist',expfile);
  return;
end

if ~exist(inmoviefile,'file'),
  msg = sprintf('Video file %s does not exist',inmoviefile);
  return;
end

% output file locations
moviefile = fullfile(expdir,moviefilestr);
trxfile = fullfile(expdir,trxfilestr);
perframedir = fullfile(expdir,perframedirstr);

%% load in experiment data

try
  expdata = load(expfile);
catch ME,
  msg = getReport(ME);
  return;
end

%% compute trajectories

if dotransposeimage,
  XIND = 1;
  YIND = 2;
else
  XIND = 2;
  YIND = 1;
end

nflies = numel(expdata.experiment_1.track);
trx = [];
perframedata = struct;
timestamps = expdata.experiment_1.elapsedTime';
dt = diff(timestamps);

allframeidx = cell(1,nflies);

for i = 1:nflies,
  
  fprintf('Larva %d...\n',i);
  
  trk = struct;
  
  % sometimes frames are skipped
  [errt,frameidx] = min(dist2(timestamps',[expdata.experiment_1.track(i).pt.et]'),[],1);
  allframeidx{i} = frameidx;
  maxerrt = max(errt);
  if maxerrt > MAXTIMESTAMPERR,
    msg = sprintf('Could not match timestamps for larva %d',i);
    return;
  end
  trk.firstframe = min(frameidx);
  trk.endframe = max(frameidx);
  trk.nframes = trk.endframe - trk.firstframe+1;
  
  if trk.nframes > numel(expdata.experiment_1.track(i).pt),
    fprintf('%d frames missing for larva %d\n',trk.nframes - numel(expdata.experiment_1.track(i).pt),i);
  end

  % centroid
  loc = double(cat(2,expdata.experiment_1.track(i).pt.loc));
  trk.x_mm = nan(1,trk.nframes);
  trk.x_mm(frameidx) = loc(1,:);
  trk.y_mm = nan(1,trk.nframes);
  trk.y_mm(frameidx) = loc(2,:);

  % grab other data about the larva
  mid = double(cat(2,expdata.experiment_1.track(i).pt.mid));
  head = double(cat(2,expdata.experiment_1.track(i).pt.head));
  tail = double(cat(2,expdata.experiment_1.track(i).pt.tail));
  spine = {expdata.experiment_1.track(i).pt.spine};
  contour = {expdata.experiment_1.track(i).pt.contour};
  
  
  % convert to pixels, as everything else will be for plotting purposes
  loc_px = cat(1,expdata.experiment_1.camcalinfo.r2cX(loc(1,:),loc(2,:)),...
    expdata.experiment_1.camcalinfo.r2cY(loc(1,:),loc(2,:)));
  trk.x = nan(1,trk.nframes);
  trk.x(frameidx) = loc_px(XIND,:)+1;
  trk.y = nan(1,trk.nframes);
  trk.y(frameidx) = loc_px(YIND,:)+1;

  % midpoint
  xmid_cm = mid(XIND,:);
  ymid_cm = mid(YIND,:);
  perframedata.xmid_mm{i} = nan(1,trk.nframes);
  perframedata.xmid_mm{i}(frameidx) = xmid_cm*SCALE;
  perframedata.ymid_mm{i} = nan(1,trk.nframes);
  perframedata.ymid_mm{i}(frameidx) = ymid_cm*SCALE;
  mid_px = cat(1,expdata.experiment_1.camcalinfo.r2cX(mid(1,:),mid(2,:)),...
    expdata.experiment_1.camcalinfo.r2cY(mid(1,:),mid(2,:)));
  trk.xmid = nan(1,trk.nframes);
  trk.xmid(frameidx) = mid_px(XIND,:)+1;
  trk.ymid = nan(1,trk.nframes);
  trk.ymid(frameidx) = mid_px(YIND,:)+1;

  % head
  xhead_cm = head(XIND,:);
  yhead_cm = head(YIND,:);
  perframedata.xhead_mm{i} = nan(1,trk.nframes);
  perframedata.xhead_mm{i}(frameidx) = xhead_cm*SCALE;
  perframedata.yhead_mm{i} = nan(1,trk.nframes);
  perframedata.yhead_mm{i}(frameidx) = yhead_cm*SCALE;
  head_px = cat(1,expdata.experiment_1.camcalinfo.r2cX(head(1,:),head(2,:)),...
    expdata.experiment_1.camcalinfo.r2cY(head(1,:),head(2,:)));
  trk.xhead = nan(1,trk.nframes);
  trk.xhead(frameidx) = head_px(XIND,:)+1;
  trk.yhead = nan(1,trk.nframes);  
  trk.yhead(frameidx) = head_px(YIND,:)+1;
  
  % tail
  xtail_cm = tail(XIND,:);
  ytail_cm = tail(YIND,:);
  perframedata.xtail_mm{i} = nan(1,trk.nframes);
  perframedata.xtail_mm{i}(frameidx) = xtail_cm*SCALE;
  perframedata.ytail_mm{i} = nan(1,trk.nframes);
  perframedata.ytail_mm{i}(frameidx) = ytail_cm*SCALE;
  tail_px = cat(1,expdata.experiment_1.camcalinfo.r2cX(tail(1,:),tail(2,:)),...
    expdata.experiment_1.camcalinfo.r2cY(tail(1,:),tail(2,:)));
  trk.xtail = nan(1,trk.nframes);
  trk.xtail(frameidx) = tail_px(XIND,:)+1;
  trk.ytail = nan(1,trk.nframes);
  trk.ytail(frameidx) = tail_px(YIND,:)+1;

  % spine
  nspinepts = max(cellfun(@(x) size(x,2),spine));
  isspine = ~cellfun(@isempty,spine);
  xspine_cm = nan(nspinepts,trk.nframes);
  yspine_cm = nan(nspinepts,trk.nframes);
  xspine_cm(:,frameidx(isspine)) = cell2mat(cellfun(@(x) double(x(1,:)'),spine(isspine),'UniformOutput',false));
  yspine_cm(:,frameidx(isspine)) = cell2mat(cellfun(@(x) double(x(2,:)'),spine(isspine),'UniformOutput',false));
  if dotransposeimage,
    trk.xspine = expdata.experiment_1.camcalinfo.r2cY(xspine_cm,yspine_cm)+1;
    trk.yspine = expdata.experiment_1.camcalinfo.r2cX(xspine_cm,yspine_cm)+1;
    trk.xspine_mm = xspine_cm*SCALE;
    trk.yspine_mm = yspine_cm*SCALE;
  else
    trk.yspine = expdata.experiment_1.camcalinfo.r2cX(xspine_cm,yspine_cm)+1;
    trk.xspine = expdata.experiment_1.camcalinfo.r2cY(xspine_cm,yspine_cm)+1;
    trk.xspine_mm = yspine_cm*SCALE;
    trk.yspine_mm = xspine_cm*SCALE;
  end
  
  % contour
  xcontour_cm = cellfun(@(x) double(x(1,:)),contour,'UniformOutput',false);
  ycontour_cm = cellfun(@(x) double(x(2,:)),contour,'UniformOutput',false);
  trk.xcontour = cell(1,trk.nframes);
  trk.ycontour = cell(1,trk.nframes);
  perframedata.xcontour_mm{i} = cell(1,trk.nframes);
  perframedata.ycontour_mm{i} = cell(1,trk.nframes);
  if dotransposeimage,
    trk.xcontour(frameidx) = cellfun(@(x,y) expdata.experiment_1.camcalinfo.r2cY(x,y)+1, xcontour_cm,ycontour_cm,'UniformOutput',false);
    trk.ycontour(frameidx) = cellfun(@(x,y) expdata.experiment_1.camcalinfo.r2cX(x,y)+1, xcontour_cm,ycontour_cm,'UniformOutput',false);
    perframedata.xcontour_mm{i}(frameidx) = cellfun(@(x) x*SCALE,xcontour_cm,'UniformOutput',false);
    perframedata.ycontour_mm{i}(frameidx) = cellfun(@(x) x*SCALE,ycontour_cm,'UniformOutput',false);
  else
    trk.ycontour(frameidx) = cellfun(@(x,y) expdata.experiment_1.camcalinfo.r2cX(x,y)+1, xcontour_cm,ycontour_cm,'UniformOutput',false);
    trk.xcontour(frameidx) = cellfun(@(x,y) expdata.experiment_1.camcalinfo.r2cY(x,y)+1, xcontour_cm,ycontour_cm,'UniformOutput',false);
    perframedata.xcontour_mm{i}(frameidx) = cellfun(@(x) x*SCALE,ycontour_cm,'UniformOutput',false);
    perframedata.ycontour_mm{i}(frameidx) = cellfun(@(x) x*SCALE,xcontour_cm,'UniformOutput',false);
  end

  % covariance matrix is already in pixels; use it to get major, minor,
  % orientation
  S = cellfun(@(x) [x(1),x(2);x(2),x(3)], {expdata.experiment_1.track(i).pt.cov},'UniformOutput',false);
  S = cat(3,S{:});
  if ~dotransposeimage,
    S = [S(2,2,:),S(2,1,:);S(1,2,:),S(1,1,:)];
  end
  [a,b,theta] = cov2ell(S);
  % note that we don't worry about the zero-indexing correction because the
  % difference will cancel out
  thetahead = atan2(head_px(YIND,:)-tail_px(YIND,:),head_px(XIND,:)-tail_px(XIND,:));
  thetaflip = modrange(thetahead+modrange(theta-thetahead,-pi/2,pi/2),-pi,pi);
  trk.a = nan(1,trk.nframes);
  trk.a(frameidx) = a/2;
  trk.b = nan(1,trk.nframes);
  trk.b(frameidx) = b/2;
  trk.theta = nan(1,trk.nframes);
  trk.theta(frameidx) = thetaflip;
  
  % convert a, b, theta to mm
  x = loc_px(XIND,:);
  y = loc_px(XIND,:);
  if dotransposeimage,
    x1 = x + a.*cos(thetaflip);
    x2 = x - a.*cos(thetaflip);
    y1 = y + a.*sin(thetaflip);
    y2 = y - a.*sin(thetaflip);
  else
    y1 = x + a.*cos(thetaflip);
    y2 = x - a.*cos(thetaflip);
    x1 = x + a.*sin(thetaflip);
    x2 = x - a.*sin(thetaflip);
  end
  x1_cm = expdata.experiment_1.camcalinfo.c2rX(x1,y1);
  y1_cm = expdata.experiment_1.camcalinfo.c2rY(x1,y1);
  x2_cm = expdata.experiment_1.camcalinfo.c2rX(x2,y2);
  y2_cm = expdata.experiment_1.camcalinfo.c2rY(x2,y2);
  trk.theta_mm = nan(1,trk.nframes);
  trk.theta_mm(frameidx) = atan2(y2_cm-y1_cm,x2_cm-x1_cm);
  trk.a_mm = nan(1,trk.nframes);
  trk.a_mm(frameidx) = sqrt((y2_cm-y1_cm).^2+(x2_cm-x1_cm).^2)/4;
  if dotransposeimage,
    x1 = x + b.*cos(thetaflip+pi/2);
    x2 = x - b.*cos(thetaflip+pi/2);
    y1 = y + b.*sin(thetaflip+pi/2);
    y2 = y - b.*sin(thetaflip+pi/2);
  else
    y1 = x + b.*cos(thetaflip+pi/2);
    y2 = x - b.*cos(thetaflip+pi/2);
    x1 = y + b.*sin(thetaflip+pi/2);
    x2 = y - b.*sin(thetaflip+pi/2);
  end
  x1_cm = expdata.experiment_1.camcalinfo.c2rX(x1,y1);
  y1_cm = expdata.experiment_1.camcalinfo.c2rY(x1,y1);
  x2_cm = expdata.experiment_1.camcalinfo.c2rX(x2,y2);
  y2_cm = expdata.experiment_1.camcalinfo.c2rY(x2,y2);
  trk.b_mm = nan(1,trk.nframes);
  trk.b_mm(frameidx) = sqrt((y2_cm-y1_cm).^2+(x2_cm-x1_cm).^2)/4;
    
  % all this was actually in centimeters, so convert to millimeters
  trk.x_mm = trk.x_mm*SCALE;
  trk.y_mm = trk.y_mm*SCALE;
  trk.a_mm = trk.a_mm*SCALE;
  trk.b_mm = trk.b_mm*SCALE;

  % approximate pxpermm
  trk.pxpermm = 1/expdata.experiment_1.camcalinfo.realUnitsPerPixel/SCALE;

  % frame rate
  trk.dt = dt(trk.firstframe:trk.endframe-1);
  %trk.dt = diff([expdata.experiment_1.track(i).pt.et]);
  
  trk.off = 1-trk.firstframe;
  trk.id = expdata.experiment_1.track(i).trackNum;
  trx = structappend(trx,trk);
  
end

%% arena parameters

switch lower(arenatype),

  case 'circle',
    arenacenterx_mm = expdata.experiment_1.camcalinfo.c2rX(arenacentery,arenacenterx)*SCALE;
    arenacentery_mm = expdata.experiment_1.camcalinfo.c2rY(arenacentery,arenacenterx)*SCALE;
    arenaradius_mm = arenaradius/trx(1).pxpermm;

    for i = 1:numel(trx),
      trx(i).arena = struct;
      trx(i).arena.arena_radius_mm = arenaradius_mm;
      trx(i).arena.arena_center_mm_x = arenacenterx_mm;
      trx(i).arena.arena_center_mm_y = arenacentery_mm;
    end
    
  case 'rectangle',
    
    tl = [arenacenterx - arenawidth/2,arenacentery - arenaheight/2];
    tr = [arenacenterx + arenawidth/2,arenacentery - arenaheight/2];
    bl = [arenacenterx - arenawidth/2,arenacentery + arenaheight/2];
    br = [arenacenterx + arenawidth/2,arenacentery + arenaheight/2];
    
    tl_mm = [expdata.experiment_1.camcalinfo.c2rX(tl(2),tl(1)),...
      expdata.experiment_1.camcalinfo.c2rY(tl(2),tl(1))]*SCALE;
    tr_mm = [expdata.experiment_1.camcalinfo.c2rX(tr(2),tr(1)),...
      expdata.experiment_1.camcalinfo.c2rY(tr(2),tr(1))]*SCALE;
    bl_mm = [expdata.experiment_1.camcalinfo.c2rX(bl(2),bl(1)),...
      expdata.experiment_1.camcalinfo.c2rY(bl(2),bl(1))]*SCALE;
    br_mm = [expdata.experiment_1.camcalinfo.c2rX(br(2),br(1)),...
      expdata.experiment_1.camcalinfo.c2rY(br(2),br(1))]*SCALE;

    for i = 1:numel(trx),
      trx(i).arena = struct;
      trx(i).arena.tl = tl_mm;
      trx(i).arena.tr = tr_mm;
      trx(i).arena.bl = bl_mm;
      trx(i).arena.br = br_mm;
    end

end


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

if strcmp(fullfile(inmoviefile),fullfile(moviefile)),
  fprintf('Input and out movie files are the same, not copying/linking.\n');
else
  
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
  
end

%% make per-frame directory
if ~exist(perframedir,'dir'),
  [success1,msg1] = mkdir(perframedir);
  if ~success1,
    msg = msg1;
    return;
  end
end

%% per-frame data

% save the perframedata we computed above
fns = fieldnames(perframedata);
% everything is currently in mm
units = parseunits('mm'); %#ok<NASGU>
for i = 1:numel(fns),
  fn = fns{i};
  data = perframedata.(fn); %#ok<NASGU>
  matfilename = fullfile(perframedir,[fn,'.mat']);
  fprintf('Saving data to %s...\n',matfilename);
  save(matfilename,'data','units');
end

%
% % for some reason, the derived measurements have different indices
% idx = cell(1,nflies);
% for i = 1:nflies,
%   idx{i} = expdata.experiment_1.track(i).getDerivedQuantity('mapptstointerped',false);
% end
% 
% fielddict = struct;
% fielddict.vel = 'ma_vel';
% fielddict.speed = 'ma_speed';
% fielddict.vnorm = 'ma_vnorm';
% fielddict.theta = 'ma_theta';
% fielddict.adjspeed = 'ma_adjspeed';
% fielddict.lrdtheta = 'ma_lrdtheta';
% fielddict.pathLength = 'ma_pathLength';
% fielddict.covRatio = 'ma_covRatio';
% fielddict.covTheta = 'ma_covTheta';
% fielddict.covMinor = 'ma_covMinor';
% fielddict.covMajor = 'ma_covMajor';
% fielddict.dcovRatio = 'ma_dcovRatio';
% fielddict.sarea = 'ma_aarea';
% 
% fnsin = fieldnames(fielddict);
% 
% for j = 1:numel(fnsin),
%   fnin = fnsin{j};
%   fnout = fielddict.(fnin);
%   matfilename = fullfile(perframedir,[fnout,'.mat']);
%   fprintf('%s -> %s...\n',fnin,matfilename);
% 
%   data = cell(1,nflies);  
%   for i = 1:nflies,
%     tmp = expdata.experiment_1.track(i).getDerivedQuantity(fnin);
%     data{i} = nan(size(tmp,1),trx(i).nframes);
%     data{i}(:,allframeidx{i}) = tmp(:,idx{i});
%   end
%   save(matfilename,'data','units');  
% end

% also save unsmoothed area
data = cell(1,nflies);
for i = 1:nflies,
  data{i} = nan(1,trx(i).nframes);
  data{i}(allframeidx{i}) = [expdata.experiment_1.track(i).pt.area];
end
matfilename = fullfile(perframedir,['area','.mat']);
save(matfilename,'data','units');

success = true;
