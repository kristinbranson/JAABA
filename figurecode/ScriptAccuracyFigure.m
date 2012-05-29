%% script for making plots used to demonstrate accuracy of detectors

%% set up path

addpath ../perframe;
addpath ../perframe/compute_perframe_features;
addpath ../misc;
addpath ../filehandling;
addpath C:\Code\Ctrax\trunk\matlab\netlab;
outfigdir = 'C:/Code/Jdetect/figures/AccuracyOut';
if ~exist(outfigdir,'dir'),
  mkdir(outfigdir);
end


%% global parameters

% sample frames
annfilestr = 'movie.ufmf.ann';
colorpos = [.7,0,0];
colorneg = [0,0,.7];
border = 16;
bg_thresh = 10/255;
wmah = .5;
frac_a_back = 1;
dist_epsilon = .1;
figpos_example = [97 886 1800 650];

% groundtruthing
ifpcolor = [0,0,.7];
fpcolor = [.3,.3,1];
ifncolor = [.7,0,0];
fncolor = [1,.3,.3];
ierrorcolor = [0,0,0];
errorcolor = [.3,.3,.3];
chancecolor = [.7,.7,.7];

TOUCH = 1;
CRABWALK = 2;
BACKUP = 3;
CHASE = 4;


%% TOUCH parameters

behavior = 'Touch';
scoresfilestr = 'scoresTouch.mat';
rootdatadir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';
experiment_name = 'pBDPGAL4U_TrpA_Rig2Plate17BowlC_20110811T133805';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Touch.xml';
mainfly = 3;
otherflies = 18;
ts = [13879,13882,13888,13892,13894,13899];
ts_overlay = [13875,13882,13888,13899];

%% TOUCH draw example frames

PlotSampleFramesWrapper(behavior,rootdatadir,experiment_name,...
  scoresfilestr,mainfly,otherflies,ts,ts_overlay,...
  'outfigdir',outfigdir,...
  'annfilestr',annfilestr,...
  'colorpos',colorpos,...
  'colorneg',colorneg,...
  'border',border,...
  'bg_thresh',bg_thresh,...
  'wmah',wmah,...
  'frac_a_back',frac_a_back,...
  'dist_epsilon',dist_epsilon,...
  'figpos',figpos_example,...
  'hfig_base',TOUCH*100);

%% TOUCH plot ground truth error as a function of training set size

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

[hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[3,4]+100*TOUCH,...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);

%% CRABWALK

behavior = 'Crabwalk';
scoresfilestr = 'scoresCrabwalk.mat';
rootdatadir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';
experiment_name = 'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110609T091905';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Crabwalk.xml';
mainfly = 19;
otherflies = [];
ts = 15383:4:15403;
ts_overlay = 15383:4:15423;
%ts = round(linspace(15370,15407,6));

%% CRABWALK draw example frames

PlotSampleFramesWrapper(behavior,rootdatadir,experiment_name,...
  scoresfilestr,mainfly,otherflies,ts,ts_overlay,...
  'outfigdir',outfigdir,...
  'annfilestr',annfilestr,...
  'colorpos',colorpos,...
  'colorneg',colorneg,...
  'border',border,...
  'bg_thresh',bg_thresh,...
  'wmah',wmah,...
  'frac_a_back',frac_a_back,...
  'dist_epsilon',dist_epsilon,...
  'figpos',figpos_example,...
  'hfig_base',CRABWALK*100);

%% CRABWALK plot ground truth error as a function of training set size

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

[hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[3,4]+100*CRABWALK,...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);


%% backup

behavior = 'Backup';
scoresfilestr = 'scoresBackup.mat';
rootdatadir = 'C:/Data/JAABA/FlyBowl';
experiment_name = 'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Backup.xml';
mainfly = 1;
otherflies = [];
ts = 12390:4:12406;
ts_overlay = [12390,12398,12406];

%% BACKUP draw example frames

PlotSampleFramesWrapper(behavior,rootdatadir,experiment_name,...
  scoresfilestr,mainfly,otherflies,ts,ts_overlay,...
  'outfigdir',outfigdir,...
  'annfilestr',annfilestr,...
  'colorpos',colorpos,...
  'colorneg',colorneg,...
  'border',border,...
  'bg_thresh',bg_thresh,...
  'wmah',wmah,...
  'frac_a_back',frac_a_back,...
  'dist_epsilon',dist_epsilon,...
  'figpos',figpos_example,...
  'hfig_base',BACKUP*100);


%% BACKUP plot ground truth error as a function of training set size

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

[hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',3:4+BACKUP*100,...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);

%% chase

behavior = 'Chase';
scoresfilestr = 'scoresChasev7.mat';
rootdatadir = 'C:/Data/JAABA/FlyBowl';
experiment_name = 'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Chase.xml';
mainfly = 20;
otherflies = [12];
ts_overlay = [6440,6461,6470,6480,6490];
ts = ts_overlay;

%% CHASE plot example frames

PlotSampleFramesWrapper(behavior,rootdatadir,experiment_name,...
  scoresfilestr,mainfly,otherflies,ts,ts_overlay,...
  'outfigdir',outfigdir,...
  'annfilestr',annfilestr,...
  'colorpos',colorpos,...
  'colorneg',colorneg,...
  'border',border,...
  'bg_thresh',bg_thresh,...
  'wmah',wmah,...
  'frac_a_back',frac_a_back,...
  'dist_epsilon',dist_epsilon,...
  'figpos',figpos_example,...
  'hfig_base',CHASE*100);


%% CHASE plot ground truth error as a function of training set size


roottraindir = 'C:/Data/JAABA/FlyBowl';
rootgroundtruthdir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';

train_experiment_names = {
  'GMR_14C07_AE_01_TrpA_Rig1Plate15BowlA_20120404T141155'
  'GMR_53B02_AE_01_TrpA_Rig1Plate15BowlC_20110930T140347'
  'GMR_71G01_AE_01_TrpA_Rig2Plate14BowlC_20110707T154934'
  'GMR_94B10_AE_01_TrpA_Rig1Plate15BowlC_20111007T155325'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928'
  'pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545'};

groundtruth_experiment_names = {
  'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110609T091905'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlD_20120210T152130'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20120202T145723'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlC_20110811T133805'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlD_20111130T094443'};


train_expdirs = cellfun(@(x) fullfile(roottraindir,x),train_experiment_names,'UniformOutput',false);
groundtruth_expdirs = cellfun(@(x) fullfile(rootgroundtruthdir,x),groundtruth_experiment_names,'UniformOutput',false);

[hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[3,4]+CHASE*100,...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);

%% wing grooming

behavior = 'WingGrooming';
scoresfilestr = 'scoresWingGrooming.mat';
rootdatadir = 'C:/Data/JAABA/FlyBowl';
experiment_name = 'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928';
configfilename = '../perframe/params/FlyBowl_WingGrooming_params.xml';
mainfly = 1;
otherflies = [];
ts_sample = 12390:4:12406;
ts_overlay = ts_sample;
%ts = round(linspace(15370,15407,6));
border = 16;
colorpos = [.7,0,0];
colorneg = [0,0,.7];
outfigdir = 'C:/Code/Jdetect/figures/AccuracyOut';


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
  'hfigs',[10,11],...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);

%% chase

behavior = 'Jump';
scoresfilestr = 'scores_Jump.mat';
rootdatadir = 'C:/Data/JAABA/FlyBowl';
experiment_name = 'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Jump.xml';
mainfly = 1;
otherflies = [];
ts_sample = 12390:4:12406;
ts_overlay = ts_sample;
%ts = round(linspace(15370,15407,6));
border = 16;
colorpos = [.7,0,0];
colorneg = [0,0,.7];
outfigdir = 'C:/Code/Jdetect/figures/AccuracyOut';


%% parameters

roottraindir = 'C:/Data/JAABA/FlyBowl';
rootgroundtruthdir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';

train_experiment_names = {
'GMR_14C08_AE_01_TrpA_Rig1Plate15BowlB_20110914T113113'
'GMR_20B07_AE_01_TrpA_Rig2Plate17BowlC_20120216T143940'
'GMR_71G01_AE_01_TrpA_Rig2Plate14BowlC_20110707T154934'
'GMR_82D11_AE_01_TrpA_Rig1Plate15BowlB_20111028T093039'
'pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110610T153218'
'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20110929T143440'
};

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
  'hfigs',[10,11],...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);

