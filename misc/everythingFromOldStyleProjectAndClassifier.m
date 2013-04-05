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
 featureLexiconFileNameRelList] = ...
  getFeatureLexiconListsFromXML();
isSameLexicon=strcmp(featureLexiconFileNameRel,featureLexiconFileNameRelList);
i=find(isSameLexicon);  
if isempty(i)
  featureLexiconName='custom';
  %featureLexiconAnimalType='';
else
  featureLexiconName=featureLexiconNameList{i};
  %featureLexiconAnimalType=featureLexiconAnimalTypeList{i};
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
everythingParams.featureLexiconName=featureLexiconName;
everythingParams.featureLexicon=featureLexicon;
everythingParams.scoresAsInput=projectParams.scoresinput;
if isempty(projectParams.behaviors.names)                             
  everythingParams.behaviorName='';
else
  everythingParams.behaviorName=projectParams.behaviors.names{1};
end

% copy the experiment dirs, labels over
if isempty(classifierParams)
  % if no classifierParams provided
  everythingParams=appendEmptyLabelsAndDefaultClassifier(everythingParams,projectParams);
else
  % The usual case: caller provided a non-empty classifierParams
  everythingParams=appendClassifierAndLabels(everythingParams,projectParams,classifierParams);
end

% set the version number
everythingParams.ver='0.5.0';

end                                    
