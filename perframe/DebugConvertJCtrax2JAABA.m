%% set up path

SetUpJAABAPath;

%% parameters

inmoviefile = '/groups/branson/bransonlab/projects/JAABA/sampledata/eric_jctrax/EH100706_GMR_71G01_AE_G-TRP-29_p1/EH100706_GMR_71G01_AE_G-TRP-29_p1_JAABA/movie.mp4';
intrxfile = '/groups/branson/bransonlab/projects/JAABA/sampledata/eric_jctrax/EH100706_GMR_71G01_AE_G-TRP-29_p1/EH100706_GMR_71G01_AE_G-TRP-29_p1-track.mat';
infofile = '/groups/branson/bransonlab/projects/JAABA/sampledata/eric_jctrax/EH100706_GMR_71G01_AE_G-TRP-29_p1/EH100706_GMR_71G01_AE_G-TRP-29_p1_info.mat';
outexpdir = '/groups/branson/bransonlab/projects/JAABA/sampledata/eric_jctrax/EH100706_GMR_71G01_AE_G-TRP-29_p1/Exp20131113T180359';
moviefilestr = 'movie.mp4';
trxfilestr = 'trx.mat';
perframedirstr = 'perframe';
dosoftlink = true;

%% run

[success,msg] = ConvertJCtrax2JAABA('inmoviefile',inmoviefile,...
  'intrxfile',intrxfile,...
  'infofile',infofile,...
  'expdir',outexpdir,...
  'moviefilestr',moviefilestr,...    
  'trxfilestr',trxfilestr,...
  'perframedirstr',perframedirstr,...    
  'dosoftlink',dosoftlink);
