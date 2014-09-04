function [success,msg] = ConvertTreadmill2JAABA(varargin)

success = false;
msg = {};

QTRMAJ = 0.6556;
QTRMIN = 0.2269;


[inmoviefile,intrxfile,...
  expdir,moviefilestr,trxfilestr,perframedirstr,...
  arenafile,dosoftlink,uniqueframesonly,resamplerate,roiradius_mm] = myparse(varargin,...
  'inmoviefile','','intrxfile','',...
  'expdir','','moviefilestr','movie.avi','trxfilestr','trx.mat','perframedirstr','perframe',...
  'arenafile','','dosoftlink',false,...
  'uniqueframesonly',false,'resamplerate',10,...
  'roiradius_mm',4);

%% check that required inputs are given
% if isempty(inmoviefile),
%   msg = 'Input movie file is empty';
%   return;
% end
if isempty(intrxfile),
  msg = 'Input trx mat file is empty';
  return;
end
if isempty(arenafile),
  msg = 'Input arena file is empty';
  return;
end
if isempty(expdir),
  msg = 'Input experiment directory is empty';
  return;
end

% output file locations
moviefile = fullfile(expdir,moviefilestr);
trxfile = fullfile(expdir,trxfilestr);
perframedir = fullfile(expdir,perframedirstr);

%% load in the trx

fid = fopen(intrxfile,'r');
isfirst = true;
fns = {'time','x','y','orientation'};
trx = struct('x',[],'y',[],'theta',[],'a',[],'b',[],'id',1,...
  'moviename',moviefile,'firstframe',1,'off',0,...
  'nframes',nan,'endframe',nan,'timestamps',[],...
  'x_mm',[],'y_mm',[],'a_mm',[],'b_mm',[],'theta_mm',[],...
  'dt',[],'fps',nan,'pxpermm',1);
while true,
  
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  if isempty(s),
    continue;
  end
  ss = regexp(s,',','split');
  if isfirst,
    headers = ss;
    for i = 1:numel(fns),
      idx.(fns{i}) = find(strcmpi(headers,fns{i}),1);
      if isempty(idx.(fns{i})),
        msg = sprintf('Could not find %s in header list',fns{i});
        return;
      end
    end
    isfirst = false;
    continue;
  end

  ss = str2double(ss);
  
  trx.x(end+1) = ss(idx.x);
  trx.y(end+1) = ss(idx.y);
  trx.theta(end+1) = modrange(ss(idx.orientation)*pi/180,-pi,pi);
  trx.a(end+1) = QTRMAJ;
  trx.b(end+1) = QTRMIN;
  trx.timestamps(end+1) = ss(idx.time);
  
end

if uniqueframesonly || ~isempty(resamplerate),
  isunique = [true,diff(trx.x) ~= 0 & diff(trx.y) ~= 0];
  trx.x = trx.x(isunique);
  trx.y = trx.y(isunique);
  trx.theta = trx.theta(isunique);
  trx.a = trx.a(isunique);
  trx.b = trx.b(isunique);
  trx.timestamps = trx.timestamps(isunique);
  trx.tidx = find(isunique);
end
if ~isempty(resamplerate),
  
  T0 = trx.timestamps(1);
  T1 = trx.timestamps(end);
  ts = T0:1/resamplerate:T1;

  trx.x = interp1(trx.timestamps,trx.x,ts);
  trx.y = interp1(trx.timestamps,trx.y,ts);
  trx.theta = modrange(interp1(trx.timestamps,unwrap(trx.theta),ts),-pi,pi);
  trx.nframes = numel(ts);
  trx.a = repmat(QTRMAJ,[1,trx.nframes]);
  trx.b = repmat(QTRMIN,[1,trx.nframes]);
  off = 1;
  % get timestamp index
  trx.tidx = nan(1,trx.nframes);
  for i = 1:trx.nframes,
    for j = off:numel(trx.timestamps),
      if trx.timestamps(j)>=ts(i),
        d2 = trx.timestamps(j)-ts(i);
        if j > 1,
          d1 = ts(i)-trx.timestamps(j-1);
        else
          d1 = inf;
        end
        if d2 < d1,
          off = j;
        else
          off = j-1;
        end
        trx.tidx(i) = off;
        break;
      end
    end
    if isnan(trx.tidx(i)),
      error('Sanity check: tidx(%d) did not get assigned',i);
    end
  end
  trx.timestamps = ts;
  
end

trx.x_mm = trx.x;
trx.y_mm = trx.y;
trx.theta_mm = trx.theta;
trx.a_mm = trx.a;
trx.b_mm = trx.b;
trx.dt = diff(trx.timestamps);
trx.nframes = numel(trx.x);
trx.endframe = trx.nframes;
trx.fps = 1/mean(trx.timestamps);

fclose(fid);

%% read in arena parameters

fid = fopen(arenafile,'r');

rois = nan(0,4);
roinames = {};
arenatype = 'circle';
arenaradius = nan;
arenacenterx = nan;
arenacentery = nan;

while true,
  s = fgetl(fid);
  if ~ischar(s),
    break;
  end
  s = strtrim(s);
  if isempty(s),
    continue;
  end
  i = find(s == '-',1);
  if isempty(i),
    error('Could not find - in line %s',s);
  end
  name = strtrim(s(1:i-1));
  vals = sscanf(s(i+1:end),'%f');

  switch name,
    case 'Origin',
      origin = vals;
    case 'Arena radius'
      if numel(vals) ~= 1,
        error('Did not find 1 value for %s: %s',name,s);
      end
      arenaradius = vals;
    case 'Arena center'
      if numel(vals) < 2,
        error('Did not find at least 2 values for %s: %s',name,s);
      end
      arenacenterx = vals(1);
      arenacentery = vals(2);
    otherwise
      if numel(vals) ~= 3,
        error('Did not find 3 values for %s: %s',name,s);
      end
      if ~isempty(roiradius_mm),
        vals(3) = roiradius_mm;
      end
      rois(end+1,:) = [vals',nan];
      roinames{end+1} = name;
  end
    
end
fclose(fid);

% update positions based on origin
trx.x_mm = trx.x_mm+origin(1);
trx.y_mm = trx.y_mm+origin(2);
trx.x = trx.x+origin(1);
trx.y = trx.y+origin(2);

trx = SetLandmarkParameters(trx,arenatype,arenacenterx,arenacentery,...
  arenaradius,[],[],rois);
trx.roinames = roinames;

timestamps = trx.timestamps;

%% create the experiment directory
if ~exist(expdir,'dir'),
  [success1,msg1] = mkdir(expdir);
  if ~success1,
    msg = msg1;
    return;
  end
  msg{end+1} = sprintf('Created experiment directory %s.',expdir);
end

% save the trx file
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

% % copy/soft-link movie
% if strcmp(fullfile(inmoviefile),fullfile(moviefile)),
%   msg{end+1} = 'Input and out movie files are the same, not copying/linking.';
% else
%   
%   if dosoftlink,
%     if exist(moviefile,'file'),
%       delete(moviefile);
%     end
%     if isunix,
%       cmd = sprintf('ln -s %s %s',inmoviefile,moviefile);
%       unix(cmd);
%       % test to make sure it worked
%       [status,result] = unix(sprintf('readlink %s',moviefile));
%       result = strtrim(result);
%       if status ~= 0 || ~strcmp(result,inmoviefile),
%         res = questdlg(sprintf('Failed to make soft link. Copy %s to %s instead?',inmoviefile,moviefile));
%         if ~strcmpi(res,'Yes'),
%           msg = sprintf('Failed to make soft link from %s to %s.',inmoviefile,moviefile);
%           return;
%         end
%         dosoftlink = false;
%       end
%     elseif ispc,
%       if exist([moviefile,'.lnk'],'file'),
%         delete([moviefile,'.lnk']);
%       end
%       cmd = sprintf('mkshortcut.vbs /target:"%s" /shortcut:"%s"',inmoviefile,moviefile);
%       fprintf('Making a Windows shortcut file at "%s" with target "%s"\n',inmoviefile,moviefile);
%       system(cmd);
%       % test to make sure that worked
%       [equalmoviefile,didfind] = GetPCShortcutFileActualPath(moviefile);
%       if ~didfind || ~strcmp(equalmoviefile,inmoviefile),
%         res = questdlg(sprintf('Failed to make shortcut. Copy %s to %s instead?',inmoviefile,moviefile));
%         if ~strcmpi(res,'Yes'),
%           msg = sprintf('Failed to make shortcut from %s to %s.',inmoviefile,moviefile);
%           return;
%         end
%         dosoftlink = false;
%       end
%     else
%       res = questdlg(sprintf('Unknown OS, cannot soft-link movie file %s. Copy instead?',inmoviefile));
%       if ~strcmpi(res,'Yes'),
%         msg = sprintf('Failed to make softlink from %s to %s.',inmoviefile,moviefile);
%         return;
%       end
%       dosoftlink = false;
%     end
%     if dosoftlink,
%       msg{end+1} = sprintf('Made a link to movie file %s at %s',inmoviefile,moviefile);
%     end
%   end
%   
%   if ~dosoftlink,
%     if ispc,
%       if exist([moviefile,'.lnk'],'file'),
%         delete([moviefile,'.lnk']);
%       end
%     end
%     if exist(moviefile,'file'),
%       delete(moviefile);
%     end
%     [success1,msg1] = copyfile(inmoviefile,moviefile);
%     if ~success1,
%       msg = msg1;
%       success = false;
%       return;
%     end
%     msg{end+1} = sprintf('Copied movie file %s to %s',inmoviefile,moviefile);
%   end
% 
% end
  
% make per-frame directory
if ~exist(perframedir,'dir'),
  [success1,msg1] = mkdir(perframedir);
  if ~success1,
    msg = msg1;
    return;
  end
end

success = true;
