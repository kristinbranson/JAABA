function ...
  everythingFileFromOldStyleProjectAndClassifierFiles(...
    everythingFileName, ...
    projectFileName, ...
    classifierFileName, ...
    gtExpDirNames)

% load the project params  
projectParams=load(projectFileName,'-mat');

% load the classifier params and (normal) labels
if isempty(classifierFileName)
  classifierParams=struct([]);
else
  classifierParams=load(classifierFileName,'-mat');
end

% load the gt labels
gtLabelFileNameLocal=projectParams.file.gt_labelfilename;
gtLabels=struct([]);
for i=1:length(gtExpDirNames)
  gtExpDirName=gtExpDirNames{i};
  gtLabelFilePathName=fullfile(gtExpDirName,gtLabelFileNameLocal);
  gtLabelsThis=load(gtLabelFilePathName,'-mat');
  if i==1 ,
    gtLabels=gtLabelsThis;
  else
    gtLabels(i)=gtLabelsThis;
  end
end

% convert to everything params
everythingParams= ...
  everythingFromOldStyleProjectAndClassifier(projectParams, ...
                                             classifierParams, ...
                                             gtExpDirNames, ...
                                             gtLabels);  %#ok
                                           
% save to a file                                           
save(everythingFileName,'-struct','everythingParams');  

end                                    

