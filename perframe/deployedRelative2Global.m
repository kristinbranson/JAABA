function glo = deployedRelative2Global(rel)

if isdeployed && ~ismac && ~isglobalpath(rel),
  path = pwd;
  glo = fullfile(path,rel);
else
  glo = rel;
end
