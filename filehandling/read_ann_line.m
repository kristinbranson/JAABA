function [s,param,value] = read_ann_line(fid)

s = fgetl(fid);
param = [];
value = [];

if s == -1
  return;
end
if strcmp(s,'end header'),
  return;
end

i = findstr(s,':');
if isempty(i),
  return;
end;
param = s(1:i-1);
value = s(i+1:end);

specialparams = {'background median','background mean',...
                 'background mad','background std',...
                 'hfnorm','roipolygons'};
isspecial = ismember(param,specialparams);

if isspecial,
  sz = str2num(value);
  value = fread(fid,sz,'uint8');
else
  tmp = str2num(value);
  if ~isempty(tmp),
    value = tmp;
  end;
end;

param = strrep(param,' ','_');