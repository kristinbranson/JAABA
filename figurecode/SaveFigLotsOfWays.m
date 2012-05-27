function outfilenames = SaveFigLotsOfWays(hfig,basename)

formats = {'pdf','tiff','svg','fig'};
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

if ismember('fig',formats),
  filename = [basename,'.fig'];
  saveas(hfig,fullfile(path,filename),'fig');
  outfilenames{end+1} = fullfile(path,filename);
end
formats = setdiff(formats,{'fig'});

for i = 1:numel(formats),
  filename = [basename,'.',formats{i}];
  try
    savefig(filename,hfig,formats{i});
    if ~isempty(path),
      movefile(filename,fullfile(path,filename));
    end
  catch ME,
    warning('Error while creating %s: %s',filename,getReport(ME));
    continue;
  end
  outfilenames{end+1} = fullfile(path,filename);
end