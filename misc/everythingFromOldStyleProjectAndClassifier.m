function everythingParams= ...
  everythingFromOldStyleProjectAndClassifier(projectParams, ...
                                             classifierParams)
                                           
everythingParams=struct();

% get the featureLexicon from the relevant file
featureLexiconFileNameRel=projectParams.file.featureconfigfile;
pathToMisc=fileparts(mfilename());
pathToJaaba=fileparts(pathToMisc);
featureLexiconFileNameAbs= ...
  fullfile(pathToJaaba,'perframe',featureLexiconFileNameRel);
featureLexicon = ReadXMLParams(featureLexiconFileNameAbs);

% Try to match up the relative feature lexicon file name with a
% feature lexicon name
[featureLexiconNameList, ...
 featureLexiconFileNameRelList, ...
 featureLexiconAnimalTypeList] = ...
  getFeatureLexiconListsFromXML();

isSameLexicon=strcmp(featureLexiconFileNameRel,featureLexiconFileNameRelList);
featureLexiconIndex=find(isSameLexicon);
if isempty(featureLexiconIndices)
  featureLexiconName='custom'
else
  featureLexiconIndex=featureLexiconIndices(1);
end

nLexicons=length(featureLexiconNameList);

for i=1:nLexicons
  
  if isequal(featureLexiconFileNameRel,featureLexiconFileNameRelList{i})
end

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
  %animalType=projectParams.behaviors.type;
  classifier.featureLexiconName='custom';
  classifier.featureLexicon=featureLexicon;
  if isfield(projectParams.windowfeatures,'windowfeaturesparams')
    classifier.windowFeaturesParams=projectParams.windowfeatures.windowfeaturesparams;
  else
    classifier.windowFeaturesParams=struct([]);    
  end
  if isfield(projectParams.windowfeatures,'basicFeatureTable')
    classifier.basicFeatureTable=projectParams.windowfeatures.basicFeatureTable;
  else
    classifier.basicFeatureTable={};    
  end
  if isfield(projectParams.windowfeatures,'featureWindowSize')
    classifier.featureWindowSize=projectParams.windowfeatures.featureWindowSize;
  else
    classifier.featureWindowSize=10;  % the default
  end    
  classifier.scoresAsInput=projectParams.scoresinput;
  %classifier.scoresAsInput=struct('classifierfile',{}, ...
  %                                'ts',{}, ...
  %                                'scorefilename',{});
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
  %classifier.animalType=projectParams.behaviors.type;
  classifier.featureLexiconName='custom';  
  classifier.featureLexicon=featureLexicon;
  classifier.windowFeaturesParams=classifierParams.windowfeaturesparams;
  %classifier.basicFeatureTable=classifierParams.basicFeatureTable;
  %classifier.featureWindowSize=classifierParams.featureWindowSize;
  classifier.scoresAsInput=projectParams.scoresinput;
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

% set the version number
everythingParams.ver='0.5.0';

end                                    
