% out = getstructarrayfield(s,field)
% s: array of structures
% field: name of a field of every s(i)
% out: concatenation of all s(i).(field)

function out = getstructarrayfield(s,field)

if length(s) == 1,
  out = s.(field);
  return;
end

out = cell(size(s));
[out{:}] = deal(s.(field));
out = cell2mat(out);