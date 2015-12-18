%% set up path

addpath ../misc;
addpath ../filehandling;

expdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlCtrax/20120220_non_olympiad_azanchir_mating_galit_CS_20120211/results/mated20120519/EXT_CantonS_1220002_None_Rig2Plate17BowlC_20120519T133539';

%% parameters

scorefns = {
  'scoresCrabwalk.mat'
  'scoresBackup.mat'
  'scores_Walk.mat'
  'scores_Stops.mat'
  'scores_pivot_head.mat'
  'scores_pivot_tail.mat'
  'scores_Righting.mat'
  'scores_AttemptedCopulation.mat'
  'scores_Jump.mat'
  'scores_Chasev7.mat'
  'scores_wingflick.mat'
  'scoresWingGrooming.mat'
  'scoresTouch.mat'
  'scores_WingExtension.mat'
  };

%%

MakeJAABAResultsMovie_PositiveClips(expdir,scorefns{i})