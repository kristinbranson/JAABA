function data = LoadPerFrameData(obj,fn,n)

islabelfile = regexp(fn,'_labels$','once');
if islabelfile,
  labelfilestr = [fn,'.mat'];
  filename = fullfile(obj.expdirs{n},labelfilestr);
  if ~exist(filename,'file'),
    error('File %s does not exist',filename);
  end
  data = obj.LoadLabelsFromFile(labelfilestr,n);
  obj.units.(fn) = parseunits('unit');
  obj.SetPerFrameData(fn,data,n);
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