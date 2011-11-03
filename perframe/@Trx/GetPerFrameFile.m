function filename = GetPerFrameFile(obj,fn,n)

dirname = fullfile(obj.expdirs{n},obj.perframedir);
if ~exist(dirname,'file'),
  [success,msg,~] = mkdir(obj.expdirs{n},obj.perframedir);
  if ~success,
    error('Error creating per-frame directory: %s',msg);
  end
end
filename = fullfile(dirname,[fn,'.mat']);
