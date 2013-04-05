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
if isempty(i) ,
  featureLexiconName='custom';
  %featureLexiconAnimalType='';
else
  featureLexiconName=featureLexiconNameList{i};
  %featureLexiconAnimalType=featureLexiconAnimalTypeList{i};
end

% Convert featureparamlist, if present, to new format
% If not present, the sublexicon is identical to the lexicon
if isfield(projectParams,'featureparamlist') ,
  sublexiconPFNames=fieldnames(projectParams.featureparamlist);
else
  sublexiconPFNames=fieldnames(featureLexicon.perframe);
end

% copy the project params over
everythingParams.behaviors=projectParams.behaviors;
fileParams=projectParams.file;
% delete some field that we don't need/want
if isfield(fileParams,'featureconfigfile') ,
  fileParams=rmfield(fileParams,'featureconfigfile');
end
everythingParams.file=fileParams;
if isfield(projectParams,'trx')
  everythingParams.trxGraphicParams=projectParams.trx;
else
  everythingParams.trxGraphicParams=projectParams.plot.trx;
end
if isfield(projectParams,'labels')
  everythingParams.labelGraphicParams=projectParams.labels;
else
  everythingParams.labelGraphicParams=projectParams.plot.labels;
end
everythingParams.featureLexiconName=featureLexiconName;
everythingParams.featureLexicon=featureLexicon;
everythingParams.sublexiconPFNames=sublexiconPFNames;
everythingParams.scoresAsInput=projectParams.scoresinput;
if isempty(projectParams.behaviors.names) ,
  everythingParams.behaviorName='';
else
  everythingParams.behaviorName=projectParams.behaviors.names{1};
end

% copy the experiment dirs, labels over
if isempty(classifierParams) ,
  % if no classifierParams provided
  everythingParams=appendEmptyLabelsAndDefaultClassifier(everythingParams,projectParams);
else
  % The usual case: caller provided a non-empty classifierParams
  everythingParams=appendClassifierAndLabels(everythingParams,projectParams,classifierParams);
end

% set the version number
everythingParams.ver='0.5.0';

end
