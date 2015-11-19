%% set up path

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',

    JCtrax_path = 'E:\Code\JCtrax';
    %FlyBowlAnalysis_path = 'E:\Code\FlyBowlAnalysis';
    %rootdatadir = 'E:\Code\Jdetect\larva\fly_data\TrainingData';
    %settingsdir = 'E:\Code\FlyBowlAnalysis\settings';
    configfilename = 'params\JLabelParams.xml';
    npool = 4;

  case 'bransonk-lw2',

    JCtrax_path = 'C:\Code\JCtrax';
    %FlyBowlAnalysis_path = 'C:\Code\FlyBowlAnalysis';
    %rootdatadir = 'C:\Code\Jdetect\larva\fly_data\TrainingData';
    %settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    configfilename = 'params\JLabelParamsKristinChase.xml';
    npool = 4;
    
  case 'bransonk-desktop',

    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    %FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    %rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/larva/fly_data/TrainingData/';
    %settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    configfilename = 'params/JLabelParams_bransonk-desktop.xml';
    npool = 8;
    
  case 'kabram-ws.janelia.priv',

    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    %FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    %rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/larva/fly_data/TrainingData/';
    %settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    configfilename = 'params/JLabelParamsMayank.xml';
    npool = 8;
    
  case 'robiea-ws'
    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    %FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    %rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    %settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    configfilename = '/groups/branson/home/robiea/Projects_data/JLabel/JLabelParams.xml';
    npool = 8;
    
  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    %FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    %rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/larva/fly_data/TrainingData/';
    %settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    configfilename = 'params/JLabelParams_bransonk-desktop.xml';
    npool = 1;
    
end

addpath(genpath(fullfile(JCtrax_path,'pdollar_toolbox')));
addpath(fullfile(JCtrax_path,'misc'));
addpath(fullfile(JCtrax_path,'filehandling'));
jlabelpath = fileparts(which('JLabel'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));
%addpath(FlyBowlAnalysis_path);

%% parameters

rootdatadir = '/groups/branson/bransonlab/projects/olympiad/HackHitData';
perframe_params = struct;
perframe_params.fov = pi;
perframe_params.thetafil = [.0625,.25,.375,.25,.0635];
perframe_params.nbodylengths_near = 2.5;
perframe_params.max_dnose2ell_anglerange = 127;

experiment_names = {
  'GMR_65C11_AE_01_TrpA_Rig2Plate14BowlD_20110407T083514'
  'pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110707T085804'
  'GMR_71G01_AE_01_TrpA_Rig2Plate14BowlC_20110707T154934'
  'GMR_26F09_AE_01_TrpA_Rig2Plate17BowlD_20110817T133217'
  'pBDPGAL4U_TrpA_Rig2Plate14BowlD_20110615T164545'
  'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928'
  'GMR_53B02_AE_01_TrpA_Rig1Plate15BowlC_20110930T140347'
  'GMR_94B10_AE_01_TrpA_Rig1Plate15BowlC_20111007T155325'
  };

perframefns = {'dnose2ell_angle_30tomin30'
  'dnose2ell_angle_min20to20'
  'dnose2ell_angle_min30to30'
  'dnose2tail'
  'nflies_close'
  'angleonclosestfly'};

%% 

parfor i = 1:numel(experiment_names),
  expdir = fullfile(rootdatadir,experiment_names{i});
  fprintf('%s...\n',experiment_names{i});
  trx = Trx('trxfilestr','registered_trx.mat','perframe_params',perframe_params);
  trx.AddExpDir(expdir);
  for j = 1:numel(perframefns),
    fprintf('  %s...\n',perframefns{j});
    trx.(perframefns{j}); %#ok<VUNUS>
  end
end