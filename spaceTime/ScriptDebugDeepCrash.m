moviename = '/home/mayank/Work/FlySpaceTime/walkMovies/SS03500_test/movie.ufmf';
trackfilename = '/home/mayank/Work/FlySpaceTime/walkMovies/SS03500_test/trx.mat';
stationary = true;
method = 'deep-sup';
flowname = 'DS';

ff = fopen('DebugDeepCrash.txt','w');

tracks = load(trackfilename);
tracks = tracks.trx;

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
fclose(fid);

blocksize = 50;

minfirst = min([tracks.firstframe]);
maxlast = max([tracks.endframe]);
nframes = maxlast-minfirst+1;
nblocks = ceil((nframes-1)/blocksize);

allftrs ={};

% compute features in parallel for different intervals of frames.
% parfor ndx = 1:nblocks
for ndx = 41:80
  [readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
  fstart = minfirst + (ndx-1)*blocksize;
  fend = min(maxlast,ndx*blocksize);
  tic;
  allftrs{ndx} = genFeatures(readfcn,headerinfo,fstart,fend,tracks,stationary,method);
  telapsed = toc;
  fclose(fid);
  fprintf(ff,'%d\n',ndx);
  
end

