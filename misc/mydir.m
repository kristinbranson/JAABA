function [res,info] = mydir(dirname,varargin)

iscond = nargin > 1;
if iscond,
  [name,daterange,isdir] = myparse(varargin,'name','','daterange',[],'isdir',[]);
end

[p,f] = myfileparts(dirname);
if ~any(f == '*'),
  p = fullfile(p,f);
end
s = dir(dirname);
res = {};
info = [];
for i = 1:numel(s),
  if ismember(s(i).name,{'.','..'}),
    continue;
  end
  if iscond,
    if ~isempty(name),
      
      m = regexp(s(i).name,name,'once');
      if isempty(m),
        continue;
      end
      
    end
    
    if ~isempty(daterange),
      
      if s(i).datenum < daterange(1),
        continue;
      end
      if numel(daterange)>=2 && s(i).datenum > daterange(2),
        continue;
      end
      
    end
    
    if ~isempty(isdir),
      if s(i).isdir ~= isdir,
        continue;
      end
    end
    
  end
  res{end+1} = fullfile(p,s(i).name); %#ok<AGROW>
  infocurr = s(i);
  infocurr.path = p;
  info = structappend(info,infocurr);
end

