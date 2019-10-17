function params = getSTParams()
% Generate default params for space time features

params = struct;
params.scale = 6;
params.npatches_x = 8;
params.npatches_y = 8;
params.psize = 12; % patch size for hog/hof
params.nbins = 8; % number of bins in hog/hof
params.wd = 0.5; % width of hog features for visualization
params.is_stationary = false;

params.optflowwinsig = 3;
params.optflowsig = 2;
params.optreliability = 1e-4;

params.blocksize = 50; % block size for parallel features computation.

params.flow_thres = 1;
params.deepscale = 4;
params.deep_thres = 3.5;
params.deepmatching_step = 8;
params.warping_iters = 2;

params.methods = {'deep-sup','hs-sup'};
params.flownames = {'DS','ff'};