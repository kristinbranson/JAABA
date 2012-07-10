function print_names = WindowFeatureParams2String(pff,wfs)

print_names = {};

if isempty(wfs),
  return;
end
iscellinput = iscell(wfs{1});
if ~iscellinput,
  wfs = {wfs};
end
n = numel(wfs);

print_names = cell(1,n);
for i = 1:n,
  print_names{i} = pff;
  for j = 1:numel(wfs{i}),
    if ischar(wfs{i}{j}),
      print_names{i} = [print_names{i},'_',wfs{i}{j}];
    else
      print_names{i} = [print_names{i},'_',num2str(wfs{i}{j})];
    end
  end
end

if ~iscellinput,
  print_names = print_names{1};
end