function JLabelBatch(expName,classifierfilename,configfilename)

data = JLabelData(configfilename);
[success,msg] = data.SetClassifierFileName(classifierfilename);
if ~success,
  warning(msg);
end

data.AddExpDir(expName);
data.PredictWholeMovie(data.nexps); 
