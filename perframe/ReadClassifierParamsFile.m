function classifierparams = ReadClassifierParamsFile(classifierparamsfiles,varargin)

[isrelativepath,fnsrelative] = myparse(varargin,...
  'isrelativepath',false,...
  'fnsrelative',{'featureparamfilename'});

if ~iscell(classifierparamsfiles)
  classifierparamsfiles = {classifierparamsfiles};
end
  
classifiermatfiles = {};
configfiles = {};
paths = {};

for filei = 1:numel(classifierparamsfiles),
  
  classifierparamsfile = classifierparamsfiles{filei};
  
  fid = fopen(classifierparamsfile,'r');
  
  pathcurr = fileparts(classifierparamsfile);
  
  while true,
    l = fgetl(fid);
    if ~ischar(l),
      break;
    end
    l = strtrim(l);
    if isempty(l), continue; end
    if l(1) == '#',
      continue;
    end
    ws = regexp(l,',','split');
    if isrelativepath,
      for j = 1:2,
        if ~isglobalpath(ws{j}),
          ws{j} = fullfile(pathcurr,ws{j});
        end
      end
    end
    classifiermatfiles{end+1} = ws{1}; %#ok<AGROW>
    configfiles{end+1} = ws{2}; %#ok<AGROW>
    paths{end+1} = pathcurr; %#ok<AGROW>
    
  end
  fclose(fid);
  
end

nbehaviors = numel(configfiles);

% names of scores files
classifierparams = [];
for i = 1:nbehaviors,
  [~,~,ext] = fileparts(configfiles{i});
  if strcmpi(ext,'.xml'),
    params = ReadXMLConfigParams(configfiles{i});
  else
    params = load(configfiles{i});
  end
  % obsolete: this was moved to ReadXMLConfigParams
%   if ~isfield(params.file,'scorefilename'),
%     if ~iscell(params.behaviors.names),
%       params.behaviors.names = {params.behaviors.names};
%     end
%     params.file.scorefilename = ['scores',sprintf('_%s',params.behaviors.names{:}),'.mat'];
%   end
  if isrelativepath,
    for j = 1:numel(fnsrelative),
      fn = fnsrelative{j};
      if isfield(params.file,fn),
        if ~isglobalpath(params.file.(fn)),
          params.file.(fn) = fullfile(paths{i},params.file.(fn));
        end
      end
    end
  end
  params.classifierfile = classifiermatfiles{i};
  params.configfile = configfiles{i};
  classifierparams = structappend(classifierparams,params);  
end