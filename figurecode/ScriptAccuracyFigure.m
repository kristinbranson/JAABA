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
WINGGROOMING = 5;
JUMP = 6;
STOP = 7;
RIGHTING = 8;
TAILPIVOTTURN = 9;
WALK = 10;

MOUSE_WALK = 11;
MOUSE_MALEFOLLOW = 12;

LARVA_HEADCAST = 13;
LARVA_CRAWL = 14;


%% TOUCH 

behavior = 'Touch';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Touch.xml';
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
configfilename = '../perframe/params/JLabelParams_FlyBowl_Crabwalk.xml';

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

try

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

catch ME  
  warning('%s: %s',behavior,getReport(ME));
end


%% backup

behavior = 'Backup';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Backup.xml';

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

try
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
catch ME
  warning('%s: %s',behavior,getReport(ME));
end

%% chase

behavior = 'Chase';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Chase.xml';

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

try
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
catch ME
  warning('%s: %s',behavior,getReport(ME));
end

%% wing grooming

behavior = 'WingGrooming';
configfilename = '../perframe/params/JLabelParams_FlyBowl_WingGrooming.xml';

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

try
  [hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[3,4]+WINGGROOMING*100,...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);
catch ME
  warning('%s: %s',behavior,getReport(ME));
end

%% jump

behavior = 'Jump';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Jump.xml';

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

try
  [hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[3,4]+JUMP*100,...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);
catch ME
  warning('%s: %s',behavior,getReport(ME));
end

%% stop

behavior = 'Stop';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Stop.xml';

roottraindir = 'C:/Data/JAABA/FlyBowl';
rootgroundtruthdir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';

% find train directories by loading in the classifier
train_experiment_names = {
  'pBDPGAL4U_TrpA_Rig2Plate17BowlB_20110805T101647'
};

% always doing groundtruthing on these
groundtruth_experiment_names = {
  'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110609T091905'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlD_20120210T152130'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20120202T145723'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlC_20110811T133805'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlD_20111130T094443'};


train_expdirs = cellfun(@(x) fullfile(roottraindir,x),train_experiment_names,'UniformOutput',false);
groundtruth_expdirs = cellfun(@(x) fullfile(rootgroundtruthdir,x),groundtruth_experiment_names,'UniformOutput',false);

try
  [hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[3,4]+STOP*100,...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);
catch ME
  warning('%s: %s',behavior,getReport(ME));
end

%% righting

behavior = 'Righting';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Righting.xml';

roottraindir = 'C:/Data/JAABA/FlyBowl';
rootgroundtruthdir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';

% find train directories by loading in the classifier
train_experiment_names = {
  'pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110610T153218'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20110929T143440'
  'GMR_14C08_AE_01_TrpA_Rig1Plate15BowlB_20110914T113113'
  'GMR_71G01_AE_01_TrpA_Rig2Plate14BowlC_20110707T154934'
  'GMR_82D11_AE_01_TrpA_Rig1Plate15BowlB_20111028T093039'
};

% always doing groundtruthing on these
groundtruth_experiment_names = {
  'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110609T091905'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlD_20120210T152130'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20120202T145723'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlC_20110811T133805'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlD_20111130T094443'};


train_expdirs = cellfun(@(x) fullfile(roottraindir,x),train_experiment_names,'UniformOutput',false);
groundtruth_expdirs = cellfun(@(x) fullfile(rootgroundtruthdir,x),groundtruth_experiment_names,'UniformOutput',false);

try
  [hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[3,4]+RIGHTING*100,...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);
catch ME
  warning('%s: %s',behavior,getReport(ME));
end

%% tail pivot turn

behavior = 'TailPivotTurn';
configfilename = '../perframe/params/JLabelParams_FlyBowl_TailPivotTurn.xml';

roottraindir = 'C:/Data/JAABA/FlyBowl';
rootgroundtruthdir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';

% find train directories by loading in the classifier
train_experiment_names = {
  'GMR_82D11_AE_01_TrpA_Rig1Plate15BowlB_20111028T093039'
  'GMR_80A01_AE_01_TrpA_Rig2Plate14BowlD_20110408T140618'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928'
  'pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545'
  'FCF_attP2_1500062_None_Rig2Plate17BowlB_20120519T165213'
};

% always doing groundtruthing on these
groundtruth_experiment_names = {
  'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110609T091905'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlD_20120210T152130'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20120202T145723'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlC_20110811T133805'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlD_20111130T094443'};


train_expdirs = cellfun(@(x) fullfile(roottraindir,x),train_experiment_names,'UniformOutput',false);
groundtruth_expdirs = cellfun(@(x) fullfile(rootgroundtruthdir,x),groundtruth_experiment_names,'UniformOutput',false);

try
  [hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[3,4]+TAILPIVOTTURN*100,...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);
catch ME
  warning('%s: %s',behavior,getReport(ME));
end

%% walk

behavior = 'Walk';
configfilename = '../perframe/params/JLabelParams_FlyBowl_Walk.xml';

roottraindir = 'C:/Data/JAABA/FlyBowl';
rootgroundtruthdir = 'C:/Data/JAABA/groundtruth_pBDPGAL4U_data';

% find train directories by loading in the classifier
train_experiment_names = {
 'pBDPGAL4U_TrpA_Rig2Plate17BowlB_20110805T101647'
 'pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110707T085804'
 'pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545'
 'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928'
 'GMR_80A01_AE_01_TrpA_Rig2Plate14BowlD_20110408T140618'
 'GMR_82D11_AE_01_TrpA_Rig1Plate15BowlB_20111028T093039'
};

% always doing groundtruthing on these
groundtruth_experiment_names = {
  'pBDPGAL4U_TrpA_Rig1Plate10BowlB_20110609T091905'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlD_20120210T152130'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlA_20120202T145723'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlC_20110811T133805'
  'pBDPGAL4U_TrpA_Rig2Plate17BowlD_20111130T094443'};


train_expdirs = cellfun(@(x) fullfile(roottraindir,x),train_experiment_names,'UniformOutput',false);
groundtruth_expdirs = cellfun(@(x) fullfile(rootgroundtruthdir,x),groundtruth_experiment_names,'UniformOutput',false);

try
  [hfigs,errordata,figfilenames,outmatfilename] = ...
  Compute_GroundTruthError_vs_TrainingSetSize(behavior,configfilename,train_expdirs,groundtruth_expdirs,...
  'hfigs',[3,4]+WALK*100,...
  'figpos',[95 550 1246 443],...
  'ifpcolor',ifpcolor,...
  'fpcolor',fpcolor,...
  'ifncolor',ifncolor,...
  'fncolor',fncolor,...
  'ierrorcolor',ierrorcolor,...
  'errorcolor',errorcolor,...
  'chancecolor',chancecolor,...
  'outfigdir',outfigdir);
catch ME
  warning('%s: %s',behavior,getReport(ME));
end