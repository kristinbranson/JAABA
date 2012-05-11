function [params,cellparams,basicTable,windowSize] = ReadPerFrameParams(filename,featureconfigfile)

if ~exist(filename,'file'),
  error('File %s does not exist',filename);
end

DOMnode= xmlread(filename);
n = DOMnode.getDocumentElement();
%   basicData = n.
%   n = params;
v = parse(n);
if isfield(v,'basicParams')
    basicTable = convertToTable(v.basicParams,featureconfigfile);
else
    basicTable = convertToTable(struct(),featureconfigfile);
end
if isfield(v,'featureWindowSize')
    windowSize = v.featureWindowSize.size;
else
    windowSize = 10;
end
if isfield(v,'params')
    params = v.params;  
else
    params = v;
end
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
    

function outTable = convertToTable(inParams,featureconfigfile)
% allpf is required in case a new "type" of feature is defined.

settings = ReadXMLParams(featureconfigfile);
allpf = fieldnames(settings.perframe);
pftypeList = {};
for pfndx = 1:numel(allpf)
  curpf = allpf{pfndx};
  curtypes  = settings.perframe.(curpf).type; 
  if ischar(curtypes)
    curT = curtypes;
    if ~any(strcmp(pftypeList,curT))
      pftypeList{end+1} = curT;
    end
  else    
    for tndx = 1:numel(curtypes)
      curT = curtypes{tndx};
      if ~any(strcmp(pftypeList,curT))
        pftypeList{end+1} = curT;
      end
    end
  end
end

types = fieldnames(inParams);
outTable = cell(numel(pftypeList),3);
for ndx = 1:numel(pftypeList)
  matches = strcmp(types,pftypeList{ndx});
  if any(matches)
    mndx = find(matches,1);
    outTable{ndx,1} = types{mndx};
    outTable{ndx,2} = inParams.(types{mndx}).mode{1};
    outTable{ndx,3} = inParams.(types{mndx}).selection{1};  
  else
    outTable{ndx,1} = pftypeList{ndx};
    outTable{ndx,2} = 'None';
    outTable{ndx,3} = 'normal';  
  end
end

