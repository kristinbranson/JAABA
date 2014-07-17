SetUpJAABAPath;

datadir = '/groups/branson/bransonlab/projects/JAABA/data/treadmill_hannah';
expname = '9cones9cyliUSphere_DBDWTBM18_8d_lConst43_heavyB_1min_withHeader';
intrxfile = fullfile(datadir,[expname,'.txt']);
arenafile = fullfile(datadir,'9cones9cyli_Coords.txt');
expdir = fullfile(datadir,[expname,'_resample']);
resamplerate = 10; % 10 Hz
roiradius_mm = 4; % radius of cones & cylinders is 4 mm

[success,msg] =  ConvertTreadmill2JAABA('intrxfile',intrxfile,'arenafile',arenafile,'expdir',expdir,...
  'resamplerate',resamplerate,'roiradius_mm',roiradius_mm);