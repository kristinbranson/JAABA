function nthreads = SetUpParPool()

% limitations:
% computation_threads >= 1
% framecache_threads >= 1
% computation_threads + framecache_threads <= NumWorkers
% computation_threads <= numCores
% framecache_threads <= numCores

c=parcluster;
numCores = feature('numCores');
NumWorkers = c.NumWorkers;
pool_exists = ~isempty(gcp('nocreate'));

if pool_exists,
  curpool = gcp('nocreate');
  curpool_size = curpool.NumWorkers;
else
  curpool_size = 0;
end

% load in rc file
isrcfile = exist('.JLabelrc.mat','file');
if isrcfile,
  rc=load('.JLabelrc.mat');  % would've been nice to be able to just call LoadRC
else
  rc = struct;
end

% if matlabpool has already been started, use this as the number of
% computation threads
if curpool_size > 0,
  rc.computation_threads = min(curpool_size,NumWorkers-1);
end

% if both computation_threads and framecache_threads are set in rc file,
% just adjust them so that they are legal values
if isfield(rc,'computation_threads') && isfield(rc,'framecache_threads'),
  rc.computation_threads = max(1,min(rc.computation_threads,numCores));
  rc.framecache_threads = max(1,min(rc.framecache_threads,numCores));
  if rc.computation_threads + rc.framecache_threads > NumWorkers,
    if NumWorkers >= 4,
      rc.computation_threads = NumWorkers-2;
      rc.framecache_threads = 2;
    else
      rc.computation_threads = max(1,NumWorkers-1);
      rc.framecache_threads = 1;
    end
  end
% if only computation threads is set, make sure computation threads is
% legal, then set framecache_threads to min numCores-1 and rest of possible
% workers
elseif isfield(rc,'computation_threads'),
  rc.computation_threads = max(1,min(rc.computation_threads,min(numCores,NumWorkers-1)));
  rc.framecache_threads = max(1,min(numCores-1,NumWorkers-rc.computation_threads));
% if only framecache threads is set, make sure framecache threads is
% legal, then set computation_threads to min numCores and rest of possible
% workers
elseif isfield(rc,'framecache_threads'),
  rc.framecache_threads = max(1,min(rc.framecache_threads,min(numCores,NumWorkers-1)));
  rc.computation_threads = max(1,min(numCores,NumWorkers-rc.framecache_threads));
else
  if NumWorkers >= 4,
    rc.framecache_threads = min(numCores,2);
    rc.computation_threads = min(numCores,NumWorkers - rc.framecache_threads);
  else
    rc.framecache_threads = 1;
    rc.computation_threads = min(numCores,NumWorkers - 1);
  end
end

nthreads.framecache_threads = rc.framecache_threads;
nthreads.computation_threads = rc.computation_threads;
fprintf('Number of threads allocated for computation: %d\n',nthreads.computation_threads);
fprintf('Number of threads allocated for display: %d\n',nthreads.framecache_threads);

if curpool_size>0 && curpool_size ~= rc.computation_threads,
  delete(gcp);
  curpool_size = numel(gcp('nocreate'));
end

% start up matlabpool if necessary
if (numCores>1) && (curpool_size<1),
  parpool(rc.computation_threads);
end

if isrcfile,
  save('.JLabelrc.mat','-struct','rc')
end
