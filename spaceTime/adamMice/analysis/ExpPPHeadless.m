function [data,stats] = ExpPPHeadless(expdirs,jabfiles,varargin)

[paramsfile,statframelimits] = myparse(varargin,'paramsfile','ExpPP.params.txt',...
  'statframelimits',[-inf,inf]);

data = ExpPP.loadexps(expdirs,jabfiles);
data = ExpPP.ensureCustomGroupInit(data);
params = ReadParams(paramsfile);

statframelimits(end) = min(statframelimits(end),max([data.trxt1]));
computeparams = {'statframelimits' statframelimits};

[data,stats] = ExpPP.computeStats(data,params{:},computeparams{:});

% possibly not the cleanest way to do this, copy-pasted
function params = ReadParams(paramsfile)

fh = fopen(paramsfile);
if isequal(fh,-1)
  error('ExpPPHeadless:errRead','Cannot read postprocessing parameter file.');
end
out = textscan(fh,'%s%s%s%s%s','delimiter',',','collectoutput',true);
fclose(fh);
out = out{1};
assert(isequal(out(1,:),{'name' 'type' 'nanok' 'default' 'desc'}));
out = out(2:end,:);
tbl = struct();
tbl.name = out(:,1);
tbl.type = out(:,2);
tbl.nanok = str2double(out(:,3));
tbl.default = str2double(out(:,4));
tbl.desc = out(:,5);
data = [tbl.name num2cell(tbl.default) tbl.desc];
params = struct();
for i = 1:size(data,1)
  prm = data{i,1};
  val = data{i,2};
  
  assert(strcmp(prm,tbl.name{i}));
  
  % check val
  if isnan(val)
    if tbl.nanok(i)
      % none
    else
      warningNoTrace('ExpPPControl:nanVal',...
        'Parameter ''%s'' cannot be NaN. Using default value of %d.\n',...
        prm,tbl.default(i));
      val = tbl.default(i);
    end
  else
    switch tbl.type{i}
      case 'nat'
        if floor(val)~=val
          warningNoTrace('ExpPPControl:fracVal',...
            'Parameter ''%s'' should be an integer. Rounding to %d.\n',prm,round(val));
          val = round(val);
        end
        if val<0
          warningNoTrace('ExpPPControl:negVal',...
            'Parameter ''%s'' cannot be negative. Using 0.\n',prm);
          val = 0;
        end
      otherwise
        assert(false,'All types currently nat');
    end
  end
  
  % assign
  subparams = regexp(prm,'#','split');
  switch numel(subparams)
    case 1
      params.(subparams{1}) = val;
    case 2
      params.(subparams{1}).(subparams{2}) = val;
    otherwise
      assert(false,'Unexpected parameter.');
  end
end

params = struct2paramscell(params);
