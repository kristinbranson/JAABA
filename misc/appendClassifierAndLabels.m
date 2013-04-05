function everythingParams=appendClassifierAndLabels(everythingParams,projectParams,classifierParams)

everythingParams.expdirs=classifierParams.expdirs;
everythingParams.labels=classifierParams.labels;
if isfield(classifierParams,'gt_labels')
  gtLabels=classifierParams.gt_labels;
else
  nExps=length(classifierParams.expdirs);
  gtLabels=struct();
  for i=1:nExps
    gtLabels(i).t0s={};
    gtLabels(i).t1s={};
    gtLabels(i).names={};
    gtLabels(i).flies=[];
    gtLabels(i).off=[];
    gtLabels(i).timestamp={};
    gtLabels(i).imp_t0s={};
    gtLabels(i).imp_t1s={};
  end  
end
everythingParams.gtLabels=gtLabels;

% copy the classifier params proper over
classifier=struct();
%classifier.animalType=projectParams.behaviors.type;
classifier.windowFeaturesParams=classifierParams.windowfeaturesparams;
if ~isempty(projectParams.behaviors.names)                             
  classifier.behaviorName=projectParams.behaviors.names{1};
else
  classifier.behaviorName='';
end
classifier.type=classifierParams.classifiertype;
classifier.params=classifierParams.classifier;
classifier.confThresholds=classifierParams.confThresholds;
classifier.scoreNorm=classifierParams.scoreNorm;
classifier.postProcessParams=classifierParams.postprocessparams;
classifier.trainingParams=classifierParams.classifier_params;
classifier.timeStamp=classifierParams.classifierTS;
everythingParams.classifier=classifier;                              

end  
