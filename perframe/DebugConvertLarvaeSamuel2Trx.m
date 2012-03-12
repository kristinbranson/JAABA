%% set up path

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-lw1',

    JCtrax_path = 'E:\Code\JCtrax';
    rootoutputdir = 'E:\Data\eric';
    rootdatadir = fullfile(rootoutputdir,'rawdata');

  case 'bransonk-lw2',

    JCtrax_path = 'C:\Code\JCtrax';
    magatanalyzer_path = 'C:\Code\MaGATAnalyzer-Matlab-Analysis';
    rootoutputdir = 'C:\Data\larvae_bruno';
    rootdatadir = fullfile(rootoutputdir,'rawdata');

  case 'bransonk-desktop',

    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    rootoutputdir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/data/eric';
    rootdatadir = fullfile(rootoutputdir,'rawdata');

  case 'kabram-ws.janelia.priv',

    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';

  case 'robiea-ws'
    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    
  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    rootoutputdir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/data/eric';
    rootdatadir = fullfile(rootoutputdir,'rawdata');

end

addpath(genpath(magatanalyzer_path));
addpath(fullfile(JCtrax_path,'misc'));
addpath(fullfile(JCtrax_path,'filehandling'));
jlabelpath = fileparts(which('JLabel'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));


%% parameters

experiment_name = '20120206_100300';
inexpdir = fullfile(rootdatadir,experiment_name);

%% 

[outexpdir] = ConvertLarvaeSamuel2Trx(inexpdir,rootoutputdir)