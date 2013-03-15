function everythingParams= ...
  everythingFromOldStyleProjectAndClassifier(projectParams, ...
                                             classifierParams)
                                           
everythingParams=struct();

% get the featureConfigParams from the relevant file
featureConfigFileNameRel=projectParams.file.featureconfigfile;
pathToMisc=fileparts(mfilename());
pathToJaaba=fileparts(pathToMisc);
featureConfigFileNameAbs= ...
  fullfile(pathToJaaba,'perframe',featureConfigFileNameRel);
featureConfigParams = ReadXMLParams(featureConfigFileNameAbs);

% copy the project params over
everythingParams.behaviors=projectParams.behaviors;
fileParams=projectParams.file;
% delete some field that we don't need/want
if isfield(fileParams,'featureconfigfile')
  fileParams=rmfield(fileParams,'featureconfigfile');
end
everythingParams.file=fileParams;
everythingParams.trxGraphicParams=projectParams.trx;
everythingParams.labelGraphicParams=projectParams.labels;

% copy the experiment dirs, labels over
if isempty(classifierParams)
  % if no classifierParams provided
  everythingParams.expdirs={};
  everythingParams.labels=struct([]);
  everythingParams.gtLabels=struct([]);
  classifier=struct();
  classifier.animalType=projectParams.behaviors.type;
  classifier.featureConfigParams=featureConfigParams;
  classifier.windowFeaturesParams=projectParams.windowfeatures.windowfeaturesparams;
  classifier.basicFeatureTable=projectParams.windowfeatures.basicFeatureTable;
  classifier.featureWindowSize=projectParams.windowfeatures.featureWindowSize;
  classifier.scoresAsInput=struct('classifierfile',{}, ...
                                  'ts',{}, ...
                                  'scorefilename',{});
  if ~isempty(projectParams.behaviors.names)                             
    classifier.behaviorName=projectParams.behaviors.names{1};
  else
    classifier.behaviorName='';
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
else
  % The usual case: caller provided a non-empty classifierParams
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
  classifier.animalType=projectParams.behaviors.type;
  classifier.featureConfigParams=featureConfigParams;
  classifier.windowFeaturesParams=classifierParams.windowfeaturesparams;
  classifier.scoresAsInput=classifierParams.scoresasinput;
  if ~isempty(projectParams.behaviors.names)                             
    classifier.behaviorName=projectParams.behaviors.names{1};
  else
    classifier.behaviorName='';
  end
  classifier.type=classifierParams.classifiertype;
  classifier.params=classifierParams.classifier;
  classifier.trainingParams=classifierParams.classifier_params;
  classifier.scoreNorm=classifierParams.scoreNorm;
  classifier.timeStamp=classifierParams.classifierTS;
  classifier.confThresholds=classifierParams.confThresholds;
  classifier.basicFeatureTable=classifierParams.basicFeatureTable;
  classifier.featureWindowSize=classifierParams.featureWindowSize;
  classifier.postProcessParams=classifierParams.postprocessparams;
  everythingParams.classifier=classifier;                              
end

% set the version number
everythingParams.ver='0.5.0';

end                                    
