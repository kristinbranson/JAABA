% obj.AddExpDir(expdir)
% Adds data associated with experiment expdir to the data represented by obj.
function AddExpDir(obj,expdir,varargin)

[dooverwrite] = myparse(varargin,'dooverwrite',true);

if ismember(expdir,obj.expdirs),
  obj.RemoveExpDir(expdir);
end

obj.nexpdirs = obj.nexpdirs + 1;
n = obj.nexpdirs;

obj.expdirs{n} = expdir;

moviename = fullfile(obj.expdirs{n},obj.moviefilestr);
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

% read trajectories
obj.trxfiles{n} = fullfile(obj.expdirs{n},obj.trxfilestr);
if ~exist(obj.trxfiles{n},'file'),
  error('Trajectory file %s does not exist',obj.trxfiles{n});
end
traj = load_tracks(obj.trxfiles{n});

% number of flies
obj.nfliespermovie(n) = length(traj);
nfliesold = obj.nflies;
obj.nflies = obj.nflies + obj.nfliespermovie(n);

% indexing flies by movie
obj.exp2flies{n} = nfliesold+1:obj.nflies;
obj.fly2exp(nfliesold+1:obj.nflies) = n;

% initialize data cache
obj.ndatacachedperexp(n) = 0;
obj.datacached{n} = struct;

% store trajectories
obj.StoreTrajectories(n,traj,dooverwrite);

% movie size in mm
obj.width_mms(n) = obj.ncs(n) / obj.pxpermm(n);
obj.height_mms(n) = obj.nrs(n) / obj.pxpermm(n);