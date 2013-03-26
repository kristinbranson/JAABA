function featureNamesNew=featureNamesNewFromOld(featureNamesOld, ...
                                                scoresAsInputOld, ...
                                                scoresAsInputNew)
% Figures out the feature names in the new feature lexicon, given the old
% feature names, the old scoresAsInput, and the new
% scoresAsInput.

% figure out which of the old were retained, which of the new are not
% present in the old
[kept,added]= ...
  setDifferencesScoresAsInput(scoresAsInputOld,scoresAsInputNew);

% Get a list of the deleted feature names 
scoreFeatureNamesOld={scoresAsInputOld(:).scorefilename};
scoreFeatureNamesDeleted=scoreFeatureNamesOld(~kept);

% Make a function that computes whether a feature name was kept
featureNameKept=@(featureName)(~ismember(featureName,scoreFeatureNamesDeleted));

% Get the feature names common to both old a new by filtering the old
% list
featureNamesCommon= ...
  cellFilter(featureNameKept,featureNamesOld);

% Get the score feature names that are truly novel
scoreBaseNameListNew={scoresAsInputNew(:).scorefilename};
featureNamesAdded=scoreBaseNameListNew(added);

% Merge the common and added lists
featureNamesNew = ...
  [featureNamesCommon;featureNamesAdded];

end
