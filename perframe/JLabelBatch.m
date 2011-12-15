function JLabelBatch(expName,classifierfilename,configfilename)

data = JLabelData(configfilename);
[success,msg] = data.SetClassifierFileName(classifierfilename);
if ~success,
  warning(msg);
end

existingNdx = find(strcmp(expName,data.expdirs));
if ~isempty(existingNdx)
  data.PredictWholeMovie(existingNdx);
else
  data.AddExpDir(expName);
  ndx = find(strcmp(expName,data.expdirs));
  data.PredictWholeMovie(ndx);
end
