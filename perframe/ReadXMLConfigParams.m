function params = ReadXMLConfigParams(filename)

params = ReadXMLParams(filename);

% backwards compatibility
if isfield(params,'file') && ~isfield(params.file,'scorefilename') && ...
    isfield(params,'behaviors') && isfield(params.behaviors,'names'),
  if iscell(params.behaviors.names),
    namescurr = params.behaviors.names;
  else
    namescurr = {params.behaviors.names};
  end
  params.file.scorefilename = ['scores_',namescurr{:},'.mat'];
  warning('scorefilename not set in config file. Setting to default value %s',params.file.scorefilename);
end
