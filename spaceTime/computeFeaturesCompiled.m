function computeFeaturesCompiled(moviename,trackfilename,stationary,...
  method,blocknum,blocksize,savename)

stationary = str2double(stationary);
blocknum = str2double(blocknum);
blocksize = str2double(blocksize);

tracks = load(trackfilename);
tracks = tracks.trx;

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
fclose(fid);

ff = [tracks.firstframe];
ee = [tracks.endframe];

minfirst = min([tracks.firstframe]);
maxlast = max([tracks.endframe]);
nframes = maxlast-minfirst+1;

% compute features in parallel for different intervals of frames.
ndx = blocknum;
[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
fstart = minfirst + (ndx-1)*blocksize;
fend = min(maxlast,ndx*blocksize+minfirst-1);
curftrs = genFeatures(readfcn,headerinfo,fstart,fend,tracks,stationary,method);
fclose(fid);

save(sprintf('%s_%d.mat',savename,blocknum),'curftrs');
