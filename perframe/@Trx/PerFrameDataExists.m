function res = PerFrameDataExists(obj,fn,n)

res = 0;

if isfield(obj.datacached{n},fn)
  res = 1;
  return;
end

filename = obj.GetPerFrameFile(fn,n);
if exist(filename,'file'),
  res = 2;
end
