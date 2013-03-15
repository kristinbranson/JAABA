function everythingParams= ...
  everythingFromOldStyleProjectAndClassifier(projectParams, ...
                                             classifierParams)
                                           
everythingParams=struct();

% get the featureConfigParams
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
  everythingParams.expdirs={};
  everythingParams.labels=struct([]);
  everythingParams.gtLabels=struct([]);
  classifier=struct();
  classifier.type='';
  classifier.featureConfigParams=featureConfigParams;
  classifier.windowFeaturesParams=struct([]);
  classifier.scoresAsInput=struct([]);
  everythingParams.classifier=classifier;                            
else
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
  classifier.featureConfigParams=featureConfigParams;
  classifier.windowFeaturesParams=classifierParams.windowfeaturesparams;
  classifier.scoresAsInput=classifierParams.scoresasinput;
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
