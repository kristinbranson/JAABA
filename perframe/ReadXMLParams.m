function params = ReadXMLParams(filename)

if ~exist(filename,'file'),
  error('File %s does not exist',filename);
end

DOMnode= xmlread(filename);
n = DOMnode.getDocumentElement();
params = parse(n);


function [out,outname] = parse(n)
    
out = struct;
outname = char(n.getNodeName());

% add attributes
a = n.getAttributes();
for i = 0:a.getLength()-1,
  name = char(a.item(i).getName());
  s = regexp(char(a.item(i).getValue()),',','split');
  if numel(s) == 1,
    s = s{1};
  end
  value = str2double(strtrim(s));
  if any(isnan(value)),
    value = s;
  end
  if isfield(out,name),
    warning('Overwriting field %s\n',name);
  end
  out.(name) = value;
end

cs = n.getChildNodes();
for i = 0:cs.getLength()-1,
  c = cs.item(i);
  if c.getNodeType() == 1,
    [in,name] = parse(c);
    if isfield(out,name),
      warning('Overwriting field %s\n',name);
    end
    out.(name) = in;
  end
end

