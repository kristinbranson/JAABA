%% set up paths

function success = ScriptTestPrepareJAABAData

success = true;

SetUpJAABAPath;
jaabapath = fileparts(which('StartJAABA'));
magatanalyzerdir = fullfile(jaabapath,'..','..','MAGATAnalyzer-Matlab-Analysis');

if ~exist(magatanalyzerdir,'dir'),
  magatanalyzerdir = uigetdir('.','Select MAGATAnalyzer-Matlab-Analysis directory');
  if ~ischar(magatanalyzerdir),
    return;
  end
end
addpath(magatanalyzerdir);
addpath(fullfile(magatanalyzerdir,'utility functions'));

if ispc,
  rootdir = 'C:\Data\JAABA\sampledata';
else
  rootdir = '/groups/branson/bransonlab/projects/JAABA/sampledata';
end
if ~exist(rootdir,'dir'),
  rootdir = uigetdir('.','Select root JAABA/sampledata directory');
  if ~ischar(rootdir),
    return;
  end
end  

%% parameters

datatypes = {
  'Ctrax'
  'CtraxPlusWings'
  'SimpleTwoFlies'
  'LarvaeRiveraAlba'
  'Qtrax'
  'MAGATAnalyzer'
  'MWT'
  'LarvaeReid'
  'LarvaeLouis'};

%% test

for i = 1:numel(datatypes),
  
  fprintf('Testing PrepareJAABAData for %s...\n',datatypes{i});
  success = success & testPrepareJAABAData(datatypes{i},rootdir);
  
end