function data = LoadPerFrameData(obj,fn,n)

islabelfile = regexp(fn,'[lL]abels','once');
isscoresfile = regexp(fn,'[sS]cores','once');
if islabelfile,
  labelfilestr = [fn,'.mat'];
  filename = fullfile(obj.expdirs{n},labelfilestr);
  if ~exist(filename,'file'),
    scoresfn = regexprep(fn,'labels','scores','preservecase','once');
    scoresfilestr = [scoresfn,'.mat'];
    scoresfilename = fullfile(obj.expdirs{n},scoresfilestr);
    if ~exist(scoresfilename,'file'),
      error('Files %s and %s do not exist',filename,scoresfilename);
    end
    [scores,data] = obj.LoadScoresFromFile(scoresfilestr,n);
    obj.units.(scoresfn) = parseunits('unit');
    obj.SetPerFrameData(scoresfn,scores,n);
  else
    data = obj.LoadLabelsFromFile(labelfilestr,n);
  end
  obj.units.(fn) = parseunits('unit');
  obj.SetPerFrameData(fn,data,n);
elseif isscoresfile,
  scoresfilestr = [fn,'.mat'];
  filename = fullfile(obj.expdirs{n},scoresfilestr);
  if ~exist(filename,'file'),
    error('File %s does not exist',filename);
  end
  [data,labels] = obj.LoadScoresFromFile(scoresfilestr,n);
  obj.units.(fn) = parseunits('unit');
  obj.SetPerFrameData(fn,data,n);  
  labelfn = regexprep(fn,'scores','labels','preservecase','once');
  obj.SetPerFrameData(labelfn,labels,n);
else
  
  filename = obj.GetPerFrameFile(fn,n);
  if exist(filename,'file')
    x = load(filename,'data','units');
  else
    [x.data,x.units] = obj.ComputePerFrameData(fn,n);
  end
  obj.units.(fn) = x.units;
  obj.SetPerFrameData(fn,x.data,n);
  data = x.data;
end