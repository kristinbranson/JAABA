function res = PerFrameDataExists(obj,fn,n)

res = 0;

if ismember(fn,obj.fnscached{n}),
  res = 1;
  return;
end

filename = obj.GetPerFrameFile(fn,n);
if exist(filename,'file'),
  res = 2;
end
