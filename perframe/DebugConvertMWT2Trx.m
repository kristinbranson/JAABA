%% set up path

[~,computername] = system('hostname');
computername = strtrim(computername);

switch computername,
 
  case 'bransonk-lw2',

    JCtrax_path = 'C:\Code\JCtrax';
    rootdatadir = 'C:\Data\larvae_mwt';
    rootoutputdir = rootdatadir;
    
  case 'bransonk-desktop',

    JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
    rootdatadir = '/groups/branson/bransonlab/projects/larva_olympiad/data/tihana/rawdata';
    rootoutputdir = rootdatadir;
    
end

addpath(fullfile(JCtrax_path,'misc'));
addpath(fullfile(JCtrax_path,'filehandling'));
jlabelpath = fileparts(which('JLabel'));
addpath(fullfile(jlabelpath,'compute_perframe_features'));


%% parameters

experiment_name = '20120305_081215';
inexpdir = fullfile(rootdatadir,experiment_name);

%% 

[outexpdir] = ConvertMWT2Trx(inexpdir,rootoutputdir);