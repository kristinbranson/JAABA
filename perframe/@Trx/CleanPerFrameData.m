function CleanPerFrameData(obj,fns,ns)

if ~exist('ns','var'),
  ns = 1:obj.nexpdirs;
end
if ~exist('fns','var') || isempty(fns),
  fns = setdiff(Trx.PerFrameFieldNames(),Trx.TrajectoryFieldNames());
end
if ~iscell(fns),
  fns = {fns};
end
for i = 1:numel(fns),
  fn = fns{i};
  for n = ns,
    filename = obj.GetPerFrameFile(fn,n);
    if exist(filename,'file'),
      fprintf('Deleting per-frame data file %s\n',filename);
      delete(filename);
    end
    
    % clear from cache
    if isfield(obj.datacached{n},fn),
      obj.datacached{n} = rmfield(obj.datacached{n},fn);
    end
    
  end
end