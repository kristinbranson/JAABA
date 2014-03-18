%% plot the relative locations of other flies when a fly performs a given behavior

%% set up path

if ispc,
  addpath ../../FlyBowlAnalysis;
else
  addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis;
end
addpath ../perframe;
addpath ../perframe/compute_perframe_features;
addpath ../misc;
addpath ../filehandling;
outfigdir = '../figures/BehaviorStimulus';
if ~exist(outfigdir,'dir'),
  mkdir(outfigdir);
end
rootdatadir0 = '../experiments/age';

if ispc,
  addpath C:\Code\FlyBowlAnalysis;
else
  addpath /groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis;
end

%%