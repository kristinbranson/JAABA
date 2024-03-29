function ftrs = computeSTFeaturesParallel(moviename,trackfilename,stationary,method,params)

if nargin<3,
  stationary = false;
end

tracks = load(trackfilename);
tracks = tracks.trx;

[readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
if fid>0
  fclose(fid);
end

blocksize = params.blocksize;

minfirst = min([tracks.firstframe]);
maxlast = max([tracks.endframe]);
nframes = maxlast-minfirst+1;
nblocks = ceil((nframes-1)/blocksize);

allftrs ={};
parfor_progress(nblocks);
% compute features in parallel for different intervals of frames.
parfor ndx = 1:nblocks
% for ndx = 1:nblocks
  [readfcn,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
  fstart = minfirst + (ndx-1)*blocksize;
  fend = min(maxlast,ndx*blocksize);
  tic;
  allftrs{ndx} = genSTFeatures(readfcn,headerinfo,fstart,fend,tracks,stationary,method,params);
  telapsed = toc;
  if fid>0
    fclose(fid);
  end
  parfor_progress;
%   fprintf('.');
%   if mod(ndx,20)==0, fprintf('\n%.2f\n',telapsed); end
end
parfor_progress(0);

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
      if numel(numframes)<4
        numframes = 1;
      else
        numframes = numframes(end);
      end
      if numframes < 1, continue; end
%       ftrs.(ff{fnum}){flynum}(:,count:count+numframes-1) = ...
%         single(reshape(allftrs{bnum}.(ff{fnum}){flynum},[],numframes));
      ftrs.(ff{fnum}){flynum}(:,:,:,count:count+numframes-1) = ...
        single(allftrs{bnum}.(ff{fnum}){flynum});
      count = count + numframes;
    end
  end
end
