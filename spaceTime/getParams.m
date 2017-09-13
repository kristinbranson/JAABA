function params = getParams

params = struct;
params.scale = 6;
params.npatches = 8;
params.psize = 12; % patch size for hog/hof
params.nbins = 8; % number of bins in hog/hof
params.patchsz = params.psize*params.npatches;
params.wd = 0.5; % width of hog features for visualization

params.optflowwinsig = 3;
params.optflowsig = 2;
params.optreliability = 1e-4;

params.blocksize = 500; % block size for parallel features computation.

params.flow_thres = 1;
params.deepscale = 4;
params.deep_thres = 3.5;
params.deepmatching_step = 8;
params.warping_iters = 2;

params.methods = {'deep-sup','hs-sup'};
params.flownames = {'DS','hs_sup'};