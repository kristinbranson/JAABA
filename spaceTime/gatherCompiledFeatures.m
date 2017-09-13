function ftrs = gatherCompiledFeatures(outdir,expdir,stationary,varargin)
% function ftrs = gatherCompiledFeatures(outdir,expdir,stationary,varargin)

[moviename,trxfilename,method] = myparse(varargin,...
  'moviename','movie.ufmf','trxfilename','trx.mat','method','deep-sup');

% [~,nframes] = get_readframe_fcn(fullfile(expdir,moviename));
[~,expname] = fileparts(expdir);
savename = fullfile(outdir,expname);
params = getParams;
allftrs = {};

tt = load(fullfile(expdir,trxfilename));
tracks = tt.trx;

ff = [tracks.firstframe];
ee = [tracks.endframe];

minfirst = min([tracks.firstframe]);
maxlast = max([tracks.endframe]);
nframes = maxlast-minfirst+1;

numblocks = ceil(nframes/params.blocksize);
for ndx = 1:numblocks
  Q = load(sprintf('%s_%d.mat',savename,ndx));
  allftrs{ndx} = Q.curftrs;
end

ff = fields(allftrs{1});

params = getParams;
mndx = find(strcmp(params.methods,method));
flowname = params.flownames{mndx};


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

extractPerframeFtrs(fullfile(expdir,'perframe'),ftrs,stationary,flowname);
