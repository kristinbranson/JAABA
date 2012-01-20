function JLabelBatch(expName,classifierfilename,configfilename)

data = JLabelData(configfilename);
data.SetClassifierFileNameBatch(classifierfilename);

existingNdx = find(strcmp(expName,data.expdirs));
if ~isempty(existingNdx)
  data.PredictWholeMovie(existingNdx);
else
  data.AddExpDir(expName);
  ndx = find(strcmp(expName,data.expdirs));
  data.PredictWholeMovie(ndx);
end

scoreFileName = sprintf('scores_%s.mat',data.labelnames{1});
sfn = fullfile(data.rootoutputdir,expName,scoreFileName);
save(sfn,'classifierfilename','-append');