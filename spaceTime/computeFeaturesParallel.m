function ftrs = computeFeaturesParallel(moviename,trackfilename,stationary,method)

if nargin<3,
  stationary = false;
end

tracks = load(trackfilename);
tracks = tracks.trx;

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
if fid>0
  fclose(fid);
end

params = getParams;
blocksize = params.blocksize;

minfirst = min([tracks.firstframe]);
maxlast = max([tracks.endframe]);
nframes = maxlast-minfirst+1;
nblocks = ceil((nframes-1)/blocksize);

allftrs ={};

% compute features in parallel for different intervals of frames.
parfor ndx = 1:nblocks
% for ndx = 1:nblocks
  [readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
  fstart = minfirst + (ndx-1)*blocksize;
  fend = min(maxlast,ndx*blocksize);
  tic;
  allftrs{ndx} = genFeatures(readfcn,headerinfo,fstart,fend,tracks,stationary,method,params);
  telapsed = toc;
  if fid>0
    fclose(fid);
  end
  fprintf('.');
  if mod(ndx,20)==0, fprintf('\n%.2f\n',telapsed); end
end

ff = fields(allftrs{1});

% Initialize the struct for features of all the frames
ftrs = struct;
for fly = 1:numel(tracks)
  frames_fly = tracks(fly).nframes;
  for fnum = 1:numel(ff)
    ftrsz = size(allftrs{1}.(ff{fnum}){1});
    ftrsz(end) = [];
    ftrs.(ff{fnum}){fly} = zeros([ftrsz frames_fly],'single');
  end
end

% assign the features from each block.
for fnum = 1:numel(ff)
  for flynum = 1:numel(tracks)
    count = 1;
    for bnum = 1:numel(allftrs)
      numframes = size(allftrs{bnum}.(ff{fnum}){flynum});
      numframes = numframes(end);
      if numframes < 1, continue; end
%       ftrs.(ff{fnum}){flynum}(:,count:count+numframes-1) = ...
%         single(reshape(allftrs{bnum}.(ff{fnum}){flynum},[],numframes));
      ftrs.(ff{fnum}){flynum}(:,:,:,count:count+numframes-1) = ...
        single(allftrs{bnum}.(ff{fnum}){flynum});
      count = count + numframes;
    end
  end
end
