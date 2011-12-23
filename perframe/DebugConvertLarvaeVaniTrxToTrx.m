%% set up path

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
  
  case 'bransonk-desktop',

    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    configfilename = 'params/JLabelParams_bransonk-desktop.xml';
    rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/data/larvae_vani';
    
  case 'bransonk-lw2',

    JCtrax_path = 'C:\Code\JCtrax';
    configfilename = 'params\JLabelParamsKristinChase.xml';
    rootdatadir = 'C:\Data\larvae_vani';

  otherwise
    
    warning('Unknown computer %s. Paths may not be setup correctly',computername);
    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    configfilename = 'params/JLabelParams_bransonk-desktop.xml';
    rootdatadir = '/groups/branson/home/bransonk/behavioranalysis/data/larvae_vani';

end

addpath(genpath(fullfile(JCtrax_path,'pdollar_toolbox')));
addpath(fullfile(JCtrax_path,'misc'));
addpath(fullfile(JCtrax_path,'filehandling'));
jlabelpath = fileparts(which('JLabel'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));


%% parameters

genotype_name = 'genotypetest';
experiment_name = '20111213-144536';

%%

expdir = fullfile(rootdatadir,genotype_name,experiment_name);

ConvertLarvaeVaniTrxToTrx(expdir);
