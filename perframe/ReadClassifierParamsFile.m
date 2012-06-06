function classifierparams = ReadClassifierParamsFile(classifierparamsfile)

fid = fopen(classifierparamsfile,'r');

classifiermatfiles = {};
configfiles = {};

while true,
  l = fgetl(fid);
  if ~ischar(l),
    break;
  end
  ws = regexp(l,',','split');
  classifiermatfiles{end+1} = ws{1}; %#ok<AGROW>
  configfiles{end+1} = ws{2}; %#ok<AGROW>
  
end
fclose(fid);

nbehaviors = numel(configfiles);

% names of scores files
classifierparams = [];
for i = 1:nbehaviors,
  params = ReadXMLParams(configfiles{i});
  if ~isfield(params.file,'scorefilename'),
    if ~iscell(params.behaviors.names),
      params.behaviors.names = {params.behaviors.names};
    end
    params.file.scorefilename = ['scores',sprintf('_%s',params.behaviors.names{:}),'.mat'];
  end
  params.classifierfile = classifiermatfiles{i};
  params.configfile = configfiles{i};
  classifierparams = structappend(classifierparams,params);  
end