function success=testAddSecondExperimentToNewWingedFlyJab()

testDataDirName='/groups/branson/bransonlab/projects/JAABA/test_data/flies_wing';
%testDataDirName='/Users/taylora/jaaba/sampledata_copy/flies_wing';
nameOfExpToAdd1='pBDPGAL4U_TrpA_Rig2Plate17BowlA_20110929T143440';
nameOfExpToAdd2='GMR_71G01_AE_01_TrpA_Rig2Plate17BowlD_20110921T084818';

% Make a Macguffin, from which we'll open a new file in JLabelData
macguffin=Macguffin('flies_wing');
macguffin.setMainBehaviorName('may22ing');
macguffin.setScoreFileName('scores_may22ing.mat');
macguffin.setTrxFileName('wingtracking_results.mat');
macguffin.setMovieFileName('movie.ufmf');

% Create the JLabelData object
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()), ...
                'isInteractive',false);
data.perframeGenerate=true;  % tell JLabelData to generate any missing perframe files

% Create a new file in the JLabelData object
data.newJabFile(macguffin);

% Add a new exp dir
[success,msg]=data.AddExpDir(fullfile(testDataDirName,nameOfExpToAdd1));
if ~success ,
  error(msg);
end

% Add a second exp dir
[success,msg]=data.AddExpDir(fullfile(testDataDirName,nameOfExpToAdd2));
if ~success ,
  error(msg);
end

% Clean up
data.closeJabFile();
data=[];  %#ok
success=true;

end

