function outfile = cuda2matlab_features(h5file,varargin)

[outdir,outfile,outfilestr] = myparse(varargin,'outdir',[],'outfile',[],'outfilestr','features.mat');

if ~ischar(outfile),
  if ~ischar(outdir),
    outdir = fileparts(h5file);
  end
  outfile = fullfile(outdir,outfilestr);
else
  outdir = fileparts(outfile);
end

if ~exist(outdir,'dir'),
  mkdir(outdir);
end

datasets = h5info(h5file).Datasets;
cudafeat = struct;
for i = 1:numel(datasets),
  n = datasets(i).Name;
  cudafeat.(n) = h5read(h5file,['/',n]);
end

curFeatures = combine_hog_hof_features(cudafeat,true);
save(outfile,'curFeatures');