function success=testJAABADetect()

% This tests compares scores that were generated earlier using JAABADetect
%%
success = false;
jabFileName='/groups/branson/home/kabram/bransonlab/projects/JAABA/test_data/test_windowdata.jab';

oldScores = load('/groups/branson/home/kabram/bransonlab/projects/JAABA/test_data/GMR_71G01_AE_01_TrpA_Rig1Plate15BowlA_20120316T144027/scores_test_gui_old.mat');

JAABADetect('/groups/branson/home/kabram/bransonlab/projects/JAABA/test_data/GMR_71G01_AE_01_TrpA_Rig1Plate15BowlA_20120316T144027',....
  'jabfiles',{jabFileName});

newScores = load('/groups/branson/home/kabram/bransonlab/projects/JAABA/test_data/GMR_71G01_AE_01_TrpA_Rig1Plate15BowlA_20120316T144027/scores_test_gui.mat');

if ~isequaln(oldScores,newScores),
  error('Score from JAABADetect are different');
end

success = true;