function [path,name,ext] = myfileparts(s)

path = '';
name = '';
if isempty(s),
  return;
end

if s(end) == '/' || ispc && s(end) == '\',
  s = s(1:end-1);
end

[path,name,ext] = fileparts_platind(s);
if nargout < 3,
  name = [name,ext];
end