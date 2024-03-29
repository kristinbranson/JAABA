% obj.AddExpDir(expdir)
% Adds data associated with experiment expdir to the data represented by obj.
function AddExpDir(obj,expdir,varargin)

[dooverwrite,openmovie,aptInfo,fromAPT] = myparse(varargin,'dooverwrite',true,...
  'openmovie',true,'aptInfo',[],'fromAPT',false);

% remove trailing /
if expdir(end) == '/' || (ispc && expdir(end) == '\'),
  expdir = expdir(1:end-1);
end

if ismember(expdir,obj.expdirs),
  obj.RemoveExpDir(expdir);
end

obj.nexpdirs = obj.nexpdirs + 1;
n = obj.nexpdirs;

obj.expdirs{n} = expdir;
[~,expname] = myfileparts(expdir);
if ischar(obj.rootwritedir),
  obj.outexpdirs{n} = fullfile(obj.rootwritedir,expname);
  if ~exist(obj.outexpdirs{n},'dir'),
    [success,msg,~] = mkdir(obj.rootwritedir,expname);
    if ~success,
      error('Could not create output write directory %s: %s',obj.outexpdirs{n},msg);
    end
  end
else
  obj.outexpdirs{n} = expdir;
end

if openmovie && ~isempty(obj.moviefilestr),
  moviename = fullfile(obj.expdirs{n},obj.moviefilestr);
  outmoviename = fullfile(obj.outexpdirs{n},obj.moviefilestr);
  if ~exist(moviename,'file') && exist(outmoviename,'file'),
    moviename = outmoviename;
  end
  obj.movienames{n} = moviename;
  if ~exist(moviename,'file'),
    error('Movie %s does not exist',moviename);
  end
  [~,nframes,fid,vidinfo] = get_readframe_fcn(moviename);
  if ~isempty(fid) && ~isnan(fid) && fid > 1,
    fclose(fid);
  end

  % store video info
  obj.nrs(n) = vidinfo.nr;
  obj.ncs(n) = vidinfo.nc;
  if isfield(vidinfo,'ncolors'),
    obj.ncolors(n) = vidinfo.ncolors;
  else
    obj.ncolors(n) = 1;
  end
  obj.movie_nframes(n) = nframes;
end

% read trajectories
trxfilename = fullfile(obj.expdirs{n},obj.trxfilestr);
outtrxfilename = fullfile(obj.outexpdirs{n},obj.trxfilestr);
% deal with windows shortcuts
if ispc && ~exist(trxfilename,'file'),
  [actualfilename,didfind] = GetPCShortcutFileActualPath(trxfilename);
  if didfind,
    trxfilename = actualfilename;
  end
end
if ~exist(trxfilename,'file') && exist(outtrxfilename,'file'),
  trxfilename = outtrxfilename;
end
if ispc && ~exist(trxfilename,'file') && ~exist(outtrxfilename,'file'),
  [actualfilename,didfind] = GetPCShortcutFileActualPath(outtrxfilename);
  if didfind,
    trxfilename = actualfilename;
  end
end

obj.trxfiles{n} = trxfilename;
if ~exist(obj.trxfiles{n},'file'),
  error('Trajectory file %s does not exist',obj.trxfiles{n});
end
traj = load_tracks(obj.trxfiles{n});

% added to deal with .dt field being stored incorrectly by FlyTracker
lastframe = max([traj.endframe]);
for i = 1:numel(traj),
  if numel(traj(i).dt) > traj(i).nframes-1,
    % dt is for every frame: know how to fix this
    if numel(traj(i).dt) == lastframe-1,
      traj(i).dt = traj(i).dt(traj(i).firstframe:traj(i).endframe-1);
    else
      warning('dt for experiment %s target %d is of an unexpected size, trying anyways!',expdir,i);
      meandt = nanmedian(traj(i).dt);
      traj(i).dt = meandt + zeros(1,traj(i).nframes-1);
    end
  end
end

% set movie properties when there is no movie
if ~openmovie || isempty(obj.moviefilestr),
  obj.nrs(n) = max([traj.y]);
  obj.ncs(n) = max([traj.x]);
  obj.ncolors(n) = 0;
  obj.movie_nframes(n) = max([traj.endframe]);
end

% number of flies
obj.nfliespermovie(n) = length(traj);
nfliesold = obj.nflies;
obj.nflies = obj.nflies + obj.nfliespermovie(n);

% indexing flies by movie
obj.exp2flies{n} = nfliesold+1:obj.nflies;
obj.fly2exp(nfliesold+1:obj.nflies) = n;

% initialize data cache
obj.ndatacachedperexp(n) = 0;
obj.datacached{n} = {};
obj.fnscached{n} = {};
obj.nfnscached(n) = 0;

% Add apt stuff if it is missing
if fromAPT && ~any(strcmp(fieldnames(traj),'kpts'))
  trkfilename = fullfile(obj.expdirs{n},obj.trkfilestr);
  for view = 1:numel(trkfilename)
    traj = addAPTTrk2Trx(traj,trkfilename{view},'aptInfo',aptInfo,'view',view);
  end
end

% store trajectories
obj.StoreTrajectories(n,traj,dooverwrite);

% movie size in mm
obj.width_mms(n) = obj.ncs(n) / obj.pxpermm(n);
obj.height_mms(n) = obj.nrs(n) / obj.pxpermm(n);

