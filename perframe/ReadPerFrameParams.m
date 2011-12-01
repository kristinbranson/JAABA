function [params,cellparams] = ReadPerFrameParams(filename)

if ~exist(filename,'file'),
  error('File %s does not exist',filename);
end

DOMnode= xmlread(filename);
n = DOMnode.getDocumentElement();
params = parse(n);

% convert to cell params also
if nargout > 1,
  cellparams = struct;
  fns1 = fieldnames(params);
  for i1 = 1:numel(fns1),
    fn1 = fns1{i1};
    fns2 = fieldnames(params.(fn1));
    cellparams.(fn1) = {};
    feature_types = {};
    for i2 = 1:numel(fns2),
      fn2 = fns2{i2};
      if ~isstruct(params.(fn1).(fn2)),
        cellparams.(fn1)(end+1:end+2) = {fn2,params.(fn1).(fn2)};
      else
        cellparams.(fn1)(end+1:end+2) = {[fn2,'_params'],struct2paramscell(params.(fn1).(fn2))};
        feature_types{end+1} = fn2; %#ok<AGROW>
      end
    end
    cellparams.(fn1)(end+1:end+2) = {'feature_types',feature_types};
  end
end

function [out,outname] = parse(n)
    
out = struct;
outname = char(n.getNodeName());

% add attributes
a = n.getAttributes();
for i = 0:a.getLength()-1,
  name = char(a.item(i).getName());
  s = regexp(char(a.item(i).getValue()),',','split');
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
    
