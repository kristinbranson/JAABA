%% script for making plots used to demonstrate accuracy of detectors

%% set up path

addpath ../perframe;
addpath ../perframe/compute_perframe_features;
addpath ../misc;
addpath ../filehandling;
outfigdir = 'C:/Code/Jdetect/figures/AccuracyOut';
if ~exist(outfigdir,'dir'),
  mkdir(outfigdir);
end

%% touch

behavior = 'Touch';
scoresfilestr = 'scoresTouch.mat';
rootdatadir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';
experiment_name = 'pBDPGAL4U_TrpA_Rig2Plate17BowlC_20110811T133805';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Touch.xml';
mainfly = 3;
otherflies = 18;
ts = [13879,13882,13888,13892,13894,13899];

colorpos = [.7,0,0];
colorneg = [0,0,.7];

%% load the data

trx = Trx('trxfilestr','registered_trx.mat','moviefilestr','movie.ufmf','perframedir','perframe');
expdir = fullfile(rootdatadir,experiment_name);
trx.AddExpDir(expdir);
moviename = trx.movienames{1};
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
scoresdata = load(fullfile(expdir,scoresfilestr));
predictions = cellfun(@(x) x > 0,scoresdata.allScores.scores,'UniformOutput',false);

%% draw example frames

hfig = 1;

hfig = PlotSampleFrames(trx,readframe,predictions,mainfly,otherflies,ts,...
  'colorpos',colorpos,...
  'colorneg',colorneg,...
  'hfig',hfig,...
  'figpos',[]);

SaveFigLotsOfWays(hfig,fullfile(outfigdir,sprintf('Example%s_%s_fly%02d_frame%02d',behavior,experiment_name,mainfly,ts(1))));

%% parameters

roottraindir = 'C:/Data/JAABA/FlyBowl';
rootgroundtruthdir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';

train_experiment_names = {
  'GMR_71G01_AE_01_TrpA_Rig2Plate17BowlD_20110921T084818'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928'
  'pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545'
  'pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110707T085804'};

groundtruth_experiment_names = {
  'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110609T091905'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlD_20120210T152130'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20120202T145723'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlC_20110811T133805'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlD_20111130T094443'};


train_expdirs = cellfun(@(x) fullfile(roottraindir,x),train_experiment_names,'UniformOutput',false);
groundtruth_expdirs = cellfun(@(x) fullfile(rootgroundtruthdir,x),groundtruth_experiment_names,'UniformOutput',false);

%% plot ground truth error as a function of training set size

ifpcolor = [0,0,.7];
fpcolor = [.3,.3,1];
ifncolor = [.7,0,0];
fncolor = [1,.3,.3];
ierrorcolor = [0,0,0];
errorcolor = [.3,.3,.3];
chancecolor = [.7,.7,.7];

[hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[2,3],...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);

%% crabwalk

behavior = 'Crabwalk';
scoresfilestr = 'scoresCrabwalk.mat';
rootdatadir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';
experiment_name = 'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110609T091905';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Crabwalk.xml';
mainfly = 19;
otherflies = [];
ts_sample = 15383:4:15403;
ts_overlay = 15383:4:15423;
%ts = round(linspace(15370,15407,6));
border = 16;
colorpos = [.7,0,0];
colorneg = [0,0,.7];

%% load the data

trx = Trx('trxfilestr','registered_trx.mat','moviefilestr','movie.ufmf','perframedir','perframe');
expdir = fullfile(rootdatadir,experiment_name);
trx.AddExpDir(expdir);
moviename = trx.movienames{1};
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviename);
scoresdata = load(fullfile(expdir,scoresfilestr));
predictions = cellfun(@(x) x > 0,scoresdata.allScores.scores,'UniformOutput',false);

%% draw example frames

hfig = 4;

hfig = PlotSampleFrames(trx,readframe,predictions,mainfly,otherflies,ts_sample,...
  'colorpos',colorpos,...
  'colorneg',colorneg,...
  'hfig',hfig,...
  'figpos',[],...
  'border',border);

SaveFigLotsOfWays(hfig,fullfile(outfigdir,sprintf('Example%s_%s_fly%02d_frame%02d',behavior,experiment_name,mainfly,ts_sample(1))));

%% overlay example frames

hfig = 5;
figpos = [97 886 523 219];
maxv_foreground = .75;
max_weight_color = .25;

hfig = PlotMeanFrames(trx,readframe,predictions,mainfly,otherflies,ts_overlay,...
  'colorpos',colorpos,...
  'colorneg',colorneg,...
  'hfig',hfig,...
  'figpos',figpos,...
  'maxv_foreground',maxv_foreground,...
  'max_weight_color',max_weight_color,...
  'border',border);

SaveFigLotsOfWays(hfig,fullfile(outfigdir,sprintf('Example%s_Overlay_%s_fly%02d_frame%02d',behavior,experiment_name,mainfly,ts_sample(1))));

%% parameters

roottraindir = 'C:/Data/JAABA/FlyBowl';
rootgroundtruthdir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';

train_experiment_names = {
  'GMR_71G01_AE_01_TrpA_Rig2Plate17BowlD_20110921T084818'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928'
  'pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545'
  'pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110707T085804'};

groundtruth_experiment_names = {
  'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110609T091905'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlD_20120210T152130'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20120202T145723'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlC_20110811T133805'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlD_20111130T094443'};


train_expdirs = cellfun(@(x) fullfile(roottraindir,x),train_experiment_names,'UniformOutput',false);
groundtruth_expdirs = cellfun(@(x) fullfile(rootgroundtruthdir,x),groundtruth_experiment_names,'UniformOutput',false);

%% plot ground truth error as a function of training set size

ifpcolor = [0,0,.7];
fpcolor = [.3,.3,1];
ifncolor = [.7,0,0];
fncolor = [1,.3,.3];
ierrorcolor = [0,0,0];
errorcolor = [.3,.3,.3];
chancecolor = [.7,.7,.7];

[hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[2,3],...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);


