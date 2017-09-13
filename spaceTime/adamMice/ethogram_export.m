function ethogram_export(hFig0)

if ~exist('hFig0','var')
  hFig0 = gcf;
end

v = version('-release');
if ~strcmp(v,'2014b')
  error('ethogram_export:version','Export is intended for use with MATLAB R2014b.');
end

% copy the entire ethogram
hFig = copyobj(hFig0,0);

% delete stuff
TAGS2DELETE = {'udUserSelectionPatch' 'udCustomGroupMarker' 'udGroupStatsAx'};
for i = 1:numel(TAGS2DELETE)
  h = findall(hFig,'userdata',TAGS2DELETE{i});
  delete(h);
end

% massage x-axis, legend
XTICK = 500;
ax = findall(hFig,'type','axes');
assert(isscalar(ax),'Multiple axes found.');
set(ax,'XTick',XTICK);
leg = findall(hFig,'type','legend');
assert(isscalar(leg),'Multiple legends found.');
set(leg,'Location','NorthOutside')

% get a location/file
exportpath = ExpPP.loadConfigVal('exportpath');
if isempty(exportpath)
  exportpath = pwd;
end
filterspec = fullfile(exportpath,'*.fig');
[fname,pname] = uiputfile(filterspec,'Export data');
if isequal(fname,0) || isequal(pname,0)
  return;
end

[~,fbase] = fileparts(fname);
figname = fullfile(pname,fname);
pdfname = fullfile(pname,[fbase '.pdf']);
epsname = fullfile(pname,[fbase '.eps']);

% save fig, print pdf, eps
hgsave(hFig,figname);
hlpExport(hFig,pdfname,'pdf');
hlpExport(hFig,epsname,'epsc2');

ExpPP.saveConfigVal('exportpath',pname);

delete(hFig);

function hlpExport(hFig,fullfname,ext)
if exist(fullfname,'file')
  [fname,pname] = uiputfile(fullfname,['Export ' ext]);
  if isequal(fname,0) || isequal(pname,0)
    fullfname = [];
  else
    fullfname = fullfile(pname,fname);
  end
end
if ~isempty(fullfname)
  print(hFig,['-d' ext],'-noui',fullfname);  
end
