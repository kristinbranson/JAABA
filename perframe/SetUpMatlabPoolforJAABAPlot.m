function nthreads = SetUpMatlabPoolforJAABAPlot()

% limitations:
% computation_threads >= 1
% framecache_threads >= 1
% computation_threads + framecache_threads <= NumWorkers
% computation_threads <= numCores
% framecache_threads <= numCores

if verLessThan('matlab','8.3.0.532'),
  
  c=parcluster;
  numCores = feature('numCores');
  NumWorkers = c.NumWorkers;
  matlabpool_size = matlabpool('size');
  
  % load in rc file
  isrcfile = exist('most_recent_config.mat','file');
  if isrcfile,
    rc=load('most_recent_config.mat');  % would've been nice to be able to just call LoadRC
  else
    rc = struct;
  end
  
  % if matlabpool has already been started, use this as the number of
  % computation threads
  if matlabpool_size > 0,
    rc.handles.computation_threads = min(matlabpool_size,NumWorkers);
  end
  
  % if both computation_threads and framecache_threads are set in rc file,
  % just adjust them so that they are legal values
  if isfield(rc,'handles') && isfield(rc.handles,'computation_threads'),
    rc.handles.computation_threads = max(1,min([rc.handles.computation_threads,numCores,NumWorkers]));
  else
    rc.handles.computation_threads = min(numCores,NumWorkers);
  end
  
  
  nthreads = rc.handles.computation_threads;
  fprintf('Number of threads allocated for computation: %d\n',nthreads);
  
  
  if matlabpool_size>0 && matlabpool_size ~= rc.handles.computation_threads,
    matlabpool close;
    matlabpool_size = matlabpool('size');
  end
  
  % start up matlabpool if necessary
  if (numCores>1) && (matlabpool_size<1),
    matlabpool('open',rc.handles.computation_threads);
  end
  
  if isrcfile,
    save('most_recent_config.mat','-struct','rc')
  end
else
  
  c=parcluster;
  numCores = feature('numCores');
  NumWorkers = c.NumWorkers;
  curpool = gcp('nocreate');
  if isempty(curpool)
    curpool_size = 0;
  else
    curpool_size = curpool.NumWorkers;
  end
  % load in rc file
  isrcfile = exist('most_recent_config.mat','file');
  if isrcfile,
    rc=load('most_recent_config.mat');  % would've been nice to be able to just call LoadRC
  else
    rc = struct;
  end
  
  % if matlabpool has already been started, use this as the number of
  % computation threads
  if curpool_size > 0,
    rc.handles.computation_threads = min(curpool_size,NumWorkers);
  end
  
  % if both computation_threads and framecache_threads are set in rc file,
  % just adjust them so that they are legal values
  if isfield(rc,'handles') && isfield(rc.handles,'computation_threads'),
    rc.handles.computation_threads = max(1,min([rc.handles.computation_threads,numCores,NumWorkers]));
  else
    rc.handles.computation_threads = min(numCores,NumWorkers);
  end
  
  
  nthreads = rc.handles.computation_threads;
  fprintf('Number of threads allocated for computation: %d\n',nthreads);
  
  
  if curpool_size>0 && curpool_size ~= rc.handles.computation_threads,
    delete(gcp);
    curpool_size = 0;
  end
  
  % start up matlabpool if necessary
  if (numCores>1) && (curpool_size<1),
    parpool(rc.handles.computation_threads);
  end
  
  if isrcfile,
    save('most_recent_config.mat','-struct','rc')
  end
  
end
  
  