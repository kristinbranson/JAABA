%% set up paths

setuppaths;
%rootdatadir = '/groups/branson/bransonlab/projects/JAABA/data/headfixedmice_adam/predta_20130520_39';
linuxroot = '/groups/branson/bransonlab';
dosoftlink = true;
doforce = false;
frontside = true;

%% pre-processing

[expdirlistfile,expdirs,expdirsCreated,expdirs_linux] = ...
  setUpDirUsingCluster({},'linuxroot',linuxroot,'doforce',doforce,'frontside',true,'dosoftlink',dosoftlink);

%% try running the function

genAllFeatures(expdirs{1},'doforce',doforce,'frontside',frontside);

%% command run on cluster

fprintf('/groups/branson/bransonlab/share/qsub_genAllFeatures.pl %s\n',expdirlistfile);

%% check that it worked

[success,problemdirs,msgs,expdirs] = ...
  CheckSetUpDir('expdirlist',expdirlistfile,...
  'rootdir','/groups/branson/bransonlab',...
  'linuxroot',linuxroot);
