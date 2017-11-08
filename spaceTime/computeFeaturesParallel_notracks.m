function ftrs = computeFeaturesParallel_notracks(moviename,method)


[~,nframes,fid,~] = get_readframe_fcn(moviename);
if fid>0
  fclose(fid);
end

params = getParamsKatie;
blocksize = params.blocksize;

minfirst = 1;
maxlast = nframes;
nframes = maxlast-minfirst+1;
nblocks = ceil((nframes-1)/blocksize);

allftrs ={};

% compute features in parallel for different intervals of frames.
parfor ndx = 1:nblocks
% for ndx = 1:nblocks
  [readfcn,~,fid,headerinfo] = get_readframe_fcn(moviename);
  fstart = minfirst + (ndx-1)*blocksize;
  fend = min(maxlast,ndx*blocksize);
  tic;
  allftrs{ndx} = genFeatures_notracks(readfcn,headerinfo,fstart,fend,method,params);
  telapsed = toc;
  if fid>0,
    fclose(fid);
  end
  fprintf('.');
  if mod(ndx,20)==0, fprintf('\n%.2f\n',telapsed); end
end

ff = fields(allftrs{1});

% Initialize the struct for features of all the frames
ftrs = struct;
frames_fly = nframes;
for fnum = 1:numel(ff)
  ftrsz = size(allftrs{1}.(ff{fnum}){1});
  ftrsz(end) = [];
  ftrs.(ff{fnum}){1} = zeros([ftrsz frames_fly],'single');
end

% assign the features from each block.
for fnum = 1:numel(ff)
  count = 1;
  for bnum = 1:numel(allftrs)
    numframes = size(allftrs{bnum}.(ff{fnum}){1});
    numframes = numframes(end);
    if numframes < 1, continue; end
%       ftrs.(ff{fnum}){flynum}(:,count:count+numframes-1) = ...
%         single(reshape(allftrs{bnum}.(ff{fnum}){flynum},[],numframes));
    ftrs.(ff{fnum}){1}(:,:,:,count:count+numframes-1) = ...
      single(allftrs{bnum}.(ff{fnum}){1});
    count = count + numframes;
  end
end
