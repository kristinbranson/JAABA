%% script demoing how to run JLabel

clear classes;

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

  case 'bransonk-lw2',

    JCtrax_path = 'C:\Code\JCtrax';
    %FlyBowlAnalysis_path = 'C:\Code\FlyBowlAnalysis';
    %rootdatadir = 'C:\Code\Jdetect\larva\fly_data\TrainingData';
    %settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
    configfilename = 'params\JLabelParams.xml';

  case 'bransonk-desktop',

    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    %FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    %rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/larva/fly_data/TrainingData/';
    %settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    configfilename = 'params/JLabelParams_bransonk-desktop.xml';
    
  case 'kabram-ws.janelia.priv',

    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    %FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    %rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/larva/fly_data/TrainingData/';
    %settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    configfilename = 'params/JLabelParamsMayank.xml';
  
  case 'robiea-ws'
    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    %FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    %rootdatadir = '/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data';
    %settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    configfilename = '/groups/branson/home/robiea/Projects_data/JLabel/JLabelParams.xml';
    
  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    %FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
    %rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/larva/fly_data/TrainingData/';
    %settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
    configfilename = 'params/JLabelParams_bransonk-desktop.xml';
    
end

addpath(genpath(fullfile(JCtrax_path,'pdollar_toolbox')));
addpath(fullfile(JCtrax_path,'misc'));
addpath(fullfile(JCtrax_path,'filehandling'));
jlabelpath = fileparts(which('JLabel'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));
%addpath(FlyBowlAnalysis_path);

%%
try 
  matlabpool(8);
catch
end
%% start JLabel

hfig = JLabel('configfilename',configfilename); handles = guidata(hfig);

