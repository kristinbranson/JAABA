function success=testJabCreationWithScoreFeature()

testDataDirName='/groups/branson/bransonlab/projects/JAABA/test_data';
%testDataDirName='/Users/taylora/jaaba/test_data';
nameOfExpToAdd='GMR_71G01_AE_01_TrpA_Rig2Plate14BowlD_20110707T154929';

% % Copy the fooing scores into the exp dir
% cmd=sprintf('cp "%s/scores_fooing_GMR_71G01_AE_01_TrpA_Rig1Plate15BowlA_20120316T144027.mat" "%s/%s/scores_fooing.mat"', ...
%             testDataDirName, ...
%             testDataDirName, ...
%             nameOfExpToAdd);
% system(cmd);          
          
% Make sure the fooing score are _not_ already in the perframe dir
scorePerframeFileNameAbs= ...
  sprintf('%s/%s/perframe/scores_fooing.mat', ...
          testDataDirName, ...
          nameOfExpToAdd);
if exist(scorePerframeFileNameAbs,'file') ,
  cmd=sprintf('rm -f "%s"',scorePerframeFileNameAbs);
  system(cmd);          
end

% Make a Macguffin, from which we'll open a new file in JLabelData
macguffin=Macguffin('flies');
macguffin.setMainBehaviorName('metafooing');
macguffin.setScoreFileName('scores_metafooing.mat');
macguffin.setTrxFileName('registered_trx.mat');
macguffin.setMovieFileName('movie.ufmf');

% Create the JLabelData object
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()), ...
                'isInteractive',false);
              
% Create a new file in the JLabelData object
data.newJabFile(macguffin);

% Add the fooing score feature
scoreFeatureJabFileNameAbs=fullfile(testDataDirName,'fooing_labels_and_classifier.jab');
fooingJLabelData=JLabelData('setstatusfn',@(~)(nop()), ...
                            'clearstatusfn',@()(nop()));
fooingJLabelData.openJabFile(scoreFeatureJabFileNameAbs,false);
timeStamp=fooingJLabelData.classifierTS;
scoreFeatureName=baseNameFromFileName(fooingJLabelData.scorefilename);
fooingJLabelData.closeJabFile();
fooingJLabelData=[];  %#ok
data.setScoreFeatures({scoreFeatureJabFileNameAbs},timeStamp,{scoreFeatureName});

% Add a new exp dir
nameOfExpDirToAdd='/groups/branson/bransonlab/projects/JAABA/test_data/GMR_71G01_AE_01_TrpA_Rig2Plate14BowlD_20110707T154929';
%nameOfExpDirToAdd='/Users/taylora/jaaba/test_data/GMR_71G01_AE_01_TrpA_Rig2Plate14BowlD_20110707T154929';
data.SetGenerateMissingFiles();  % tell JLabelData to generate any missing perframe files
[success,msg]=data.AddExpDir(nameOfExpDirToAdd);
if ~success ,
  error(msg);
end

% Check the status table
perframeDirExists=data.FileExists('perframedir');
if perframeDirExists ,
  fprintf('The JLabelData object says that the perframe directory exists, as it should.\n');
else
  error('The JLabelData object says that the perframe directory does not exist, which is not as it should be.');
end

% Add the score feature to the current vocabulary
featureLexicon = data.featureLexicon;
scoreFeatures = data.scoreFeatures;
toBeCalculatedPFNames=data.allperframefns;
windowFeatureParams = data.GetPerframeParams();
fv= ...
  FeatureVocabularyForSelectFeatures(featureLexicon, ...
                                     scoreFeatures, ...
                                     toBeCalculatedPFNames, ...
                                     windowFeatureParams);
fv.setPFToWFAmount('scores_fooing','normal');
fv.setPFEnablement('scores_fooing',true);
windowFeaturesParams=fv.getInJLabelDataFormat();
data.setWindowFeaturesParams(windowFeaturesParams);

% Change the current target to the first fly of the new experiment
iExp=whichstr(nameOfExpDirToAdd,data.expdirs);
if isempty(iExp)
  error('Added experiment seems to be missing.');
elseif length(iExp)>1
  error('Added experiment seems to be present more than once in the JLabelData object.');
end
iFly=1;
data.setCurrentTarget(iExp,iFly);

% Add a label
ts1=5000+(0:50);  % frame indices to set
iBehavior1=whichstri('metafooing',data.labelnames);
isImportant1=true;
data.SetLabel(iExp,iFly,ts1,iBehavior1,isImportant1);

% Add a label
ts2=5100+(0:50);  % frame indices to set
iBehavior2=whichstri('none',data.labelnames);
isImportant2=true;
data.SetLabel(iExp,iFly,ts2,iBehavior2,isImportant2);

% Train
data.Train();

% The predictions should be correct for the labeled bouts
data.Predict(iExp,iFly,ts1(1),ts1(end));  % Make sure the prediction has been done
pred1=data.GetPredictedIdx(iExp,iFly,ts1(1),ts1(end));
if all(pred1.predictedidx==iBehavior1)
  fprintf('Predictions are correct for first labeled bout!\n');
else
  error('Predictions are wrong for first label');
end
data.Predict(iExp,iFly,ts2(1),ts2(end));  % Make sure the prediction has been done
pred2=data.GetPredictedIdx(iExp,iFly,ts2(1),ts2(end));
if all(pred2.predictedidx==iBehavior2)
  fprintf('Predictions are correct for second labeled bout!\n');
else  
  error('Predictions are wrong for second label');
end

% Clean up
data.closeJabFile();
data=[];  %#ok
success=true;

end
