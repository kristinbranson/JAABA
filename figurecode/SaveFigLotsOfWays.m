function outfilenames = SaveFigLotsOfWays(hfig,basename)

if exist('savefig1','var'),
  savefigfun = @savefig1;
else
  savefigfun = @savefig;
end

% formats = {'pdf','tiff','svg','fig','eps'};
formats = {'pdf','tiff','png','fig','eps'};
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

tmpbasename = regexprep(basename,'[^a-zA-Z_0-9]','_');
for i = 1:numel(formats),
  filename = [basename,'.',formats{i}];
  tmpfilename = [tmpbasename,'.',formats{i}];
  try
    savefigfun(tmpfilename,hfig,formats{i});
    if ~isempty(path),
      movefile(tmpfilename,fullfile(path,filename));
    end
  catch ME,
    warning('Error while creating %s: %s',filename,getReport(ME));
    continue;
  end
  outfilenames{end+1} = fullfile(path,filename);
end