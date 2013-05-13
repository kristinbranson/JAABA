function success=testConvertClassifierWithScoreFeature()

testDataDirName='/groups/branson/bransonlab/projects/JAABA/test_data';
%testDataDirName='/Users/taylora/jaaba/test_data';

projectFileName=fullfile(testDataDirName,'metafooing_project.mat');
classifierFileName=fullfile(testDataDirName,'metafooing_labels_and_classifier.mat');
gtExpDirNames={};               
jabFileName='/tmp/metafooing_labels_and_classifier.jab';

nameOfExpInMetafooingClassifier='GMR_71G01_AE_01_TrpA_Rig2Plate14BowlD_20110707T154929';
          
% Make sure the fooing score are _not_ already in the perframe dir
scorePerframeFileNameAbs= ...
  sprintf('%s/%s/perframe/scores_fooing.mat', ...
          testDataDirName, ...
          nameOfExpInMetafooingClassifier);
if exist(scorePerframeFileNameAbs,'file') ,
  cmd=sprintf('rm "%s"',scorePerframeFileNameAbs);
  system(cmd);          
end

% Also deleted the converted .jab file if it exsits
if exist(jabFileName,'file') ,
  cmd=sprintf('rm "%s"',jabFileName);
  system(cmd);          
end

% Convert the old-style classifier
gum=GrandlyUnifyModel();
gum.setProjectFileName(projectFileName);
gum.setClassifierFileName(classifierFileName);
for i=1:length(gtExpDirNames)
  gum.addGTExpDirName(gtExpDirNames{i})
end
gum.convert(jabFileName);
gum=[];  %#ok

% Open the converted file
gtMode=false;
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()));
data.openJabFile(jabFileName,gtMode);

% Display the JLabelData object
data

% Do a few spot-checks
labels=data.labels;
if length(labels(1).t0s{1}) ~= 2 ,
  data.closeJabFile();
  error('.jab file has the wrong number of labels');
end
if abs(labels(1).timestamp{1}(2)-735364.459128824) > 1e-6 ,
  data.closeJabFile();
  error('the one timestamp I checked is wrong');
end

% Clean-up
data.closeJabFile();
success=true;

end
