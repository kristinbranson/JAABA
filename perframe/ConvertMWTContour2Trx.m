function [outexpdir] = ConvertMWTContour2Trx(experiment_name,rootoutputdir,varargin)

%% set parameters

% default parameters
rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/data/larvae_mwt/rawdata';
contourext = 'outline';
trxfilestr = 'trx.mat';
perframedirstr = 'perframe';
MINDT = .001;
pxpermm = 11.111;
maxlength_mm = 5;
DEBUG = false;
makesoftlink = true;


[rootdatadir,contourext,trxfilestr,perframedirstr,DEBUG,maxlength_mm,pxpermm,makesoftlink] = ...
  myparse(varargin,...
  'rootdatadir',rootdatadir,...
  'contourext',contourext,...
  'trxfilestr',trxfilestr,...
  'perframedirstr',perframedirstr,...
  'debug',DEBUG,...
  'maxlength',maxlength_mm,...
  'pxpermm',pxpermm,...
  'makesoftlink',makesoftlink);

maxlength_px = ceil(maxlength_mm * pxpermm);

%% names of files

% find contour files for each trajectory
files = dir(fullfile(rootdatadir,[experiment_name,'.*.',contourext]));
contournames = {};
ids = [];
for i = 1:numel(files),
  match = regexp(files(i).name,['\.(\d+)\.',contourext,'$'],'tokens','once');
  if isempty(match),
    continue;
  end
  contournames{end+1} = fullfile(rootdatadir,files(i).name);
  ids(end+1) = str2double(match{1});
end

if isempty(contournames),
  error('Could not find any contours matching the pattern %s.*.%s',fullfile(rootdatadir,experiment_name),contourext);
end

ncontours = numel(contournames);

%% put all the contours into one trx file

% read everthing in from the contour files
alltimestamps = {};
trx = [];
arena = struct('x',nan,'y',nan,'r',nan);
for i = 1:ncontours,

  trk = struct('id',ids(i),'xcontour',{{}},'ycontour',{{}});
  timestamps = [];
  fid = fopen(contournames{i},'r');
  j = 1;
  while true,
    s = fgetl(fid);
    if ~ischar(s),
      break;
    end
    d = sscanf(s,'%f');
    if numel(d) < 3,
      warning('Contour line %s length < 3, skipping this line',s);
      continue;
    end
    if mod(numel(d),2) ~= 1,
      warning('Contour line %s length is not odd, skipping this line',s);
      continue;
    end
    d = d';
    timestamp = d(1);
    x = d(2:2:end-1);
    y = d(3:2:end);
    trk.xcontour{j} = x;
    trk.ycontour{j} = y;
    timestamps(j) = timestamp;
    j = j + 1;
  end
  fclose(fid);
  trk.arena = arena;
  alltimestamps{i} = timestamps;
  trx = structappend(trx,trk);
end

% make all the timestamps agree with each other
timestamps = union(alltimestamps{:});
dts = diff(timestamps);
% sanity check to look for rounding errors
if any(dts < MINDT),
  error('Rounding error check failed: there are timestamps < %f apart',MINDT);
end
t0s = cellfun(@(x) x(1), alltimestamps);
t1s = cellfun(@(x) x(end), alltimestamps);
fps = nanmean(dts);
for i = 1:ncontours,
  trx(i).firstframe = find(t0s(i) >= timestamps,1);
  trx(i).endframe = find(t1s(i) >= timestamps,1,'last');
  trx(i).off = 1 - trx(i).firstframe;
  % sanity check to make sure this matches the number of frames
  n = trx(i).endframe - trx(i).firstframe + 1;
  if numel(trx(i).xcontour) ~= n,
    error('Number of frames check failed: Number of frames of contour data does not match start and end timestamps');
  end
  trx(i).dt = dts(trx(i).firstframe:trx(i).endframe-1);
  trx(i).nframes = n;
  trx(i).fps = fps;
end

%% fit ellipses

[XGRID,YGRID] = meshgrid(0:maxlength_px-1,0:maxlength_px-1);

for i = 1:ncontours,
  
  trx(i).x = nan(1,trx(i).nframes);
  trx(i).y = nan(1,trx(i).nframes);
  trx(i).x_mm = nan(1,trx(i).nframes);
  trx(i).y_mm = nan(1,trx(i).nframes);
  trx(i).a = nan(1,trx(i).nframes);
  trx(i).b = nan(1,trx(i).nframes);
  trx(i).a_mm = nan(1,trx(i).nframes);
  trx(i).b_mm = nan(1,trx(i).nframes);
  trx(i).theta = nan(1,trx(i).nframes);
  trx(i).theta_mm = nan(1,trx(i).nframes);
  
  for t = 1:trx(i).nframes,
    
    % create an image 
    xc = round(trx(i).xcontour{t}*pxpermm);
    yc = round(trx(i).ycontour{t}*pxpermm);
    minx = min(xc);
    maxx = max(xc);
    miny = min(yc);
    maxy = max(yc);
    nx = maxx-minx+1;
    ny = maxy-miny+1;
    bw = [inpolygon(XGRID(1:ny,1:nx),YGRID(1:ny,1:nx),[xc,xc(1)]-minx,[yc,yc(1)]-miny);false(maxlength_px-ny,nx)];
    
    % fit an ellipse
    x = XGRID(bw)+minx;
    y = YGRID(bw)+miny;
    mux = nanmean(x);
    muy = nanmean(y);
    S = cov([x,y],1);
    [a,b,theta] = cov2ell(S);
    % try and make theta be as continuous as possible
    if t > 1,
      dsame = abs(modrange(trx(i).theta(t-1)-theta,-pi,pi));
      dflip = abs(modrange(trx(i).theta(t-1)-theta+pi,-pi,pi));
      if dflip < dsame,
        theta = modrange(theta+pi,-pi,pi);
      end
    end
    
    % store quarter-major, quarter-minor axes
    a = a/2;
    b = b/2;
    mux_mm = mux / pxpermm;
    muy_mm = muy / pxpermm;
    a_mm = a / pxpermm;
    b_mm = b / pxpermm;
    trx(i).x(t) = mux_mm;
    trx(i).y(t) = muy_mm;
    trx(i).x_mm(t) = mux_mm;
    trx(i).y_mm(t) = muy_mm;
    trx(i).a(t) = a_mm;
    trx(i).b(t) = b_mm;
    trx(i).a_mm(t) = a_mm;
    trx(i).b_mm(t) = b_mm;
    trx(i).theta(t) = theta;
    trx(i).theta_mm(t) = theta;
    
    if DEBUG,
      hold off;
      plot(trx(i).xcontour{t},trx(i).ycontour{t},'k.-');
      hold on;
      drawellipse(trx(i).x(t),trx(i).y(t),trx(i).theta(t),trx(i).a(t)*2,trx(i).b(t)*2,'r');
      axis equal;
      title(num2str(t));
      drawnow;
    end
      
  end
  
end

%% plot the trajectories

if DEBUG,

nframes = numel(timestamps);
colors = jet(ncontours)*.7;
for t = 1:nframes,
  hold off;
  for i = 1:ncontours,
    if t > trx(i).endframe || t < trx(i).firstframe, 
      continue;
    end
    j = t + trx(i).off;
    plot(trx(i).xcontour{j},trx(i).ycontour{j},'.-','color',colors(i,:));
    hold on;
    drawellipse(trx(i).x(j),trx(i).y(j),trx(i).theta(j),trx(i).a(j)*2,trx(i).b(j)*2,'color','r');
  end
  axis equal;
  drawnow;
end

end

%% create output directory

% make sure the root output directory exists
if ~exist(rootoutputdir,'dir'),
  mkdir(rootoutputdir);
end

% create the experiment directory
outexpdir = fullfile(rootoutputdir,experiment_name);
if ~exist(outexpdir,'dir'),
  mkdir(outexpdir);
end

%% create the trx file

outmatname = fullfile(outexpdir,trxfilestr);
save(outmatname,'trx');

%% copy over the contour files
for i = 1:ncontours,
  [~,basename] = myfileparts(contournames{i});
  outcontourname = fullfile(outexpdir,basename);
  if isunix && makesoftlink,
    cmd = sprintf('ln -s %s %s',contournames{i},outcontourname);
    unix(cmd);
  else
    [success,msg] = copyfile(contournames{i},outcontourname);
    if ~success,
      error('Error copying file %s to %s: %s',moviename,outmoviename,msg);
    end
  end
end

%% save the per-frame features
perframedir = fullfile(outexpdir,perframedirstr);
if ~exist(perframedir,'dir'),
  mkdir(perframedir);
end

% to do: save choreography's per-frame features
