function success=testClassifierImportFlyIntoWingedFly()

% perframeFileName='/Users/taylora/jaaba/sampledata_copy/larvae/OregonRA/perframe/area.mat';
% if exist(perframeFileName,'file') ,
%   system(sprintf('mv "%s" "%s.hidden"',perframeFileName,perframeFileName));
% end

testDataDirName='/groups/branson/bransonlab/projects/JAABA/test_data';
jabFileName=fullfile(testDataDirName,'abiding_winged_fly.jab');  % winged fly .jab file
nameOfJabFileWithClassifierToImport=fullfile(testDataDirName,'may22ing_fly.jab');  % fly .jab file
nameOfExpDir='pBDPGAL4U_TrpA_Rig2Plate17BowlA_20110929T143440';
gtMode=false;
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()), ...
                'isInteractive',false);
data.openJabFile(jabFileName,gtMode);
data.perframeGenerate=true;  % generate missing per-frame files as needed
data.perframeOverwrite=false;  % don't regenerate existing per-frame files

iExp=whichstr(nameOfExpDir,data.expnames);
if isempty(iExp)
  error('The experiment seems to be missing.');
elseif length(iExp)>1
  error('The experiment seems to be present more than once in the JLabelData object.');
end
iTarget=20;
data.setCurrentTarget(iExp,iTarget);

% The prediction should be abiding for these frames
ts=17100:17110;
abidingBehavior=whichstri('abiding',data.labelnames);
data.Predict(iExp,iTarget,ts(1),ts(end));  % Make sure the prediction has been done
pred=data.GetPredictedIdx(iExp,iTarget,ts(1),ts(end));
if all(pred.predictedidx==abidingBehavior)
  %fprintf('Predictions are correct prior to importing classifier!\n');
else
  error('Predictions are wrong prior to importing classifier');
end

% The prediction should be none for these frames
ts=17140:17150;
noneBehavior=whichstri('none',data.labelnames);
data.Predict(iExp,iTarget,ts(1),ts(end));  % Make sure the prediction has been done
pred=data.GetPredictedIdx(iExp,iTarget,ts(1),ts(end));
if all(pred.predictedidx==noneBehavior)
  fprintf('Predictions are correct prior to importing classifier!\n');
else
  error('Predictions are wrong prior to importing classifier');
end

% Do the import
data.importClassifier(nameOfJabFileWithClassifierToImport);

% The prediction should be none for these frames post-import
ts=17100:17110;
noneBehavior=whichstri('none',data.labelnames);
data.Predict(iExp,iTarget,ts(1),ts(end));  % Make sure the prediction has been done
pred=data.GetPredictedIdx(iExp,iTarget,ts(1),ts(end));
if all(pred.predictedidx==noneBehavior)
  %fprintf('Predictions are correct prior to importing classifier!\n');
else
  error('Predictions are wrong after importing classifier');
end

% The prediction should be abiding for these frames post-import
ts=17140:17150;
abidingBehavior=whichstri('abiding',data.labelnames);
data.Predict(iExp,iTarget,ts(1),ts(end));  % Make sure the prediction has been done
pred=data.GetPredictedIdx(iExp,iTarget,ts(1),ts(end));
if all(pred.predictedidx==abidingBehavior)
  fprintf('Predictions are correct after importing classifier!\n');
else
  error('Predictions are wrong after importing classifier');
end

data.closeJabFile();
data=[];  %#ok
success=true;

end

