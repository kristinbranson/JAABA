function JLabelBatch(expName,classifierfilename,configfilename)

data = JLabelData(configfilename);
data.SetClassifierFileNameWoExp(classifierfilename);

existingNdx = find(strcmp(expName,data.expdirs));
if ~isempty(existingNdx)
  data.PredictWholeMovie(existingNdx);
  ndx = existingNdx;
else
  data.AddExpDir(expName);
  ndx = find(strcmp(expName,data.expdirs));
  data.PredictWholeMovie(ndx);
end

sfn = data.GetFile('scores',ndx);
save(sfn,'classifierfilename','-append');