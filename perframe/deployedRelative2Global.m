function glo = deployedRelative2Global(rel)

if isdeployed && ~isunix && ~isglobalpath(rel),
  path = pwd;
  glo = fullfile(path,rel);
else
  glo = rel;
end
