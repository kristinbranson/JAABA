function outfilenames = SaveFigLotsOfWays(hfig,basename)

formats = {'pdf','tiff','svg','ai'};
[path,basename] = myfileparts(basename);


outfilenames = {};
if ismember('svg',formats),
  filename = [basename,'.svg'];
  plot2svg(fullfile(path,filename),hfig,'png');
  outfilenames{end+1} = fullfile(path,filename);
end
formats = setdiff(formats,{'svg'});

if ismember('ai',formats),
  filename = [basename,'.ai'];
  saveas(hfig,fullfile(path,filename),'ai');
  outfilenames{end+1} = fullfile(path,filename);
end
formats = setdiff(formats,{'ai'});
  
for i = 1:numel(formats),
  filename = [basename,'.',formats{i}];
  savefig(filename,hfig,formats{i});
  if ~isempty(path),
    movefile(filename,fullfile(path,filename));
  end
  outfilenames{end+1} = fullfile(path,filename);
end