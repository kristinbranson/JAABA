function success=testClassifierImportTransSpeciesReversed()

testDataDirName='/groups/branson/bransonlab/projects/JAABA/test_data';
jabFileName=fullfile(testDataDirName,'larping_different_features.jab');  % fly .jab file
nameOfJabFileWithClassifierToImport=fullfile(testDataDirName,'may19ing_with_classifier.jab');  % larvae .jab file
% jabFileName='~/jaaba/larping_different_features.jab';  % fly .jab file
% nameOfJabFileWithClassifierToImport='~/jaaba/may19ing_with_classifier.jab';  % larvae .jab file
%nameOfExpDir='OregonRA';
gtMode=false;
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()), ...
                'isInteractive',false);
data.openJabFile(jabFileName,gtMode);
data.perframeGenerate=true;  % generate missing per-frame files as needed
data.perframeOverwrite=false;  % don't regenerate existing per-frame files

% iExp=whichstr(nameOfExpDir,data.expnames);
% if isempty(iExp)
%   error('The experiment seems to be missing.');
% elseif length(iExp)>1
%   error('The experiment seems to be present more than once in the JLabelData object.');
% end
% iTarget=1;
% data.setCurrentTarget(iExp,iTarget);

% % The prediction should be larping for these frames
% ts=3900:3940;
% larpingBehavior=whichstri('larping',data.labelnames);
% pred=data.GetPredictedIdx(iExp,iTarget,ts(1),ts(end));
% if all(pred.predictedidx==larpingBehavior)
%   fprintf('Predictions are correct prior to importing classifier!\n');
% else
%   error('Predictions are wrong prior to importing classifier');
% end

% Do the import --- this should fail, but with a particular error
try
  data.importClassifier(nameOfJabFileWithClassifierToImport);
catch excp
  if isequal(excp.identifier,'JLabelData:unableToGeneratePerframeFile')
    fprintf('Import failed gracefully, with the expected exception.  (This is a good thing.)\n');
  else
    rethrow(excp);
  end
end

% % Do a prediction
% data.predictForCurrentTargetAndTimeSpan();

% % The prediction should be none for these frames
% noneBehavior=whichstri('none',data.labelnames);
% pred=data.GetPredictedIdx(iExp,iTarget,ts(1),ts(end));
% if all(pred.predictedidx==noneBehavior)
%   fprintf('Predictions are correct after importing classifier!\n');
% else
%   error('Predictions are wrong after importing classifier');
% end

data.closeJabFile();
data=[];  %#ok
success=true;

end

