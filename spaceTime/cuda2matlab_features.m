function outfile = cuda2matlab_features(varargin)

[outdir,outfile,outfilestr,...
  h5file,indir,infrontfilestr,insidefilestr] = ...
  myparse(varargin,'outdir',[],'outfile',[],...
  'outfilestr','features.mat','h5file',[],...
  'indir',[],'infrontfilestr','hoghof_avg_front_biasjaaba.csv',...
  'insidefilestr','hoghof_avg_side_biasjaaba.csv');

if ischar(h5file),
  indir = fileparts(h5file);
  datasets = h5info(h5file).Datasets;
  cudafeat = struct;
  for i = 1:numel(datasets),
    n = datasets(i).Name;
    cudafeat.(n) = h5read(h5file,['/',n]);
  end
  curFeatures = combine_hog_hof_features(cudafeat,true);
else
  assert(ischar(indir));
  cudafeat.front = csvread(fullfile(indir,infrontfilestr));
  cudafeat.side = csvread(fullfile(indir,insidefilestr));

  nfront = size(cudafeat.front,2);
  nside = size(cudafeat.side,2);
  curFeatures = [cudafeat.side(:,nside/2+1:end),...
    cudafeat.front(:,nside/2+1:end),...
    cudafeat.side(:,1:nside/2),...
    cudafeat.front(:,1:nside/2)];
end

if ~ischar(outfile),
  if ~ischar(outdir),
    outdir = indir;
  end
  outfile = fullfile(outdir,outfilestr);
else
  outdir = fileparts(outfile);
end

if ~exist(outdir,'dir'),
  mkdir(outdir);
end


save(outfile,'curFeatures');