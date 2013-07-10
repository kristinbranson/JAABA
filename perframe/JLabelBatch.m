function JLabelBatch(expName,classifierfilename,configfilename)

data = JLabelData(configfilename);
data.SetClassifierFileNameWoExp(classifierfilename);

existingNdx = find(strcmp(expName,data.expdirs));
if ~isempty(existingNdx)
  data.PredictWholeMovie(existingNdx);
  ndx = existingNdx;
else
  [success,msg] = data.AddExpDir(expName);
  if ~success, fprintf(msg), return, end
  
  if ~data.filesfixable,
    fprintf('Experiment %s is missing required files that cannot be generated within this interface.',expName); 
    return;
  end

  if data.filesfixable && ~data.allfilesexist,
    [success,msg] = data.GenerateMissingFiles(data.nexps,false);
    if ~success,
      fprintf('Error generating missing required files for experiment %s: %s. Removing...',expdir,msg);
      return;
    end
    
    [success,msg] = data.PreLoadPeriLabelWindowData();
    if ~success,
      fprintf('Error computing window data for experiment %s: %s. Removing...',expdir,msg);
      return;
    end
  end
  
  ndx = find(strcmp(expName,data.expdirs));
end

allScores = data.PredictWholeMovie(ndx);
sfn = data.GetFile('scores',ndx);
timestamp = now;
save(sfn,'allScores','timestamp','classifierfilename');