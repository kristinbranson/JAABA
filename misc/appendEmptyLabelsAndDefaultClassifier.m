function everythingParams=appendEmptyLabelsAndDefaultClassifier(everythingParams,projectParams)

everythingParams.expdirs={};
everythingParams.labels=struct([]);
everythingParams.gtLabels=struct([]);
classifier=struct();
%animalType=projectParams.behaviors.type;
if isfield(projectParams.windowfeatures,'windowfeaturesparams')
  classifier.windowFeaturesParams=projectParams.windowfeatures.windowfeaturesparams;
else
  classifier.windowFeaturesParams=struct([]);    
end
classifier.type='boosting';  % e.g., 'boosting'
classifier.params=struct([]);
classifier.confThresholds=[0 0];
classifier.scoreNorm=[];
classifier.postProcessParams=struct([]);
classifier.trainingParams=struct('iter',100, ...
                                 'iter_updates',10, ...
                                 'numSample',2500, ...
                                 'numBins',30, ...
                                 'CVfolds',7, ...
                                 'baseClassifierTypes','Decision Stumps', ...
                                 'baseClassifierSelected',1);
classifier.timeStamp=0;
everythingParams.classifier=classifier;                            

end
