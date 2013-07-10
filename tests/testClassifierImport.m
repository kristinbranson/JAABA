function success=testClassifierImport()

testDataDirName='/groups/branson/bransonlab/projects/JAABA/test_data';
jabFileName=fullfile(testDataDirName,'larping.jab');
nameOfJabFileWithClassifierToImport=fullfile(testDataDirName,'larping_different_features.jab');
nameOfExpDir='GMR_71G01_AE_01_TrpA_Rig1Plate15BowlA_20120316T144027';
gtMode=false;
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()), ...
                'isInteractive',true);
data.openJabFile(jabFileName,gtMode);

iExp=whichstr(nameOfExpDir,data.expnames);
if isempty(iExp)
  error('The experiment seems to be missing.');
elseif length(iExp)>1
  error('The experiment seems to be present more than once in the JLabelData object.');
end
iFly=1;
data.setCurrentTarget(iExp,iFly);

% The prediction should be larping for these frames
ts=3900:3940;
larpingBehavior=whichstri('larping',data.labelnames);
data.Predict(iExp,iFly,ts(1),ts(end));  % Make sure the prediction has been done
pred=data.GetPredictedIdx(iExp,iFly,ts(1),ts(end));
if all(pred.predictedidx==larpingBehavior)
  fprintf('Predictions are correct prior to importing classifier!\n');
else
  error('Predictions are wrong prior to importing classifier');
end

% Do the import
data.importClassifier(nameOfJabFileWithClassifierToImport);

% The prediction should be none for these frames
noneBehavior=whichstri('none',data.labelnames);
data.Predict(iExp,iFly,ts(1),ts(end));  % Make sure the prediction has been done
pred=data.GetPredictedIdx(iExp,iFly,ts(1),ts(end));
if all(pred.predictedidx==noneBehavior)
  fprintf('Predictions are correct after importing classifier!\n');
else
  error('Predictions are wrong after importing classifier');
end

data.closeJabFile();
data=[];  %#ok
success=true;

end
