function [useclassifierfile,classifier] = checkForClassifierFile(jabfile,x)

if ~exist(jabfile,'file'),
  error('jabfile %s does not exist',jabfile);
end

isclassifierinjab = ~isempty(x.classifierStuff.params);

% check if the classifier file exists
[bdir,jabnn,jabext] = fileparts(jabfile);
classifiernn = [jabnn '_classifier.mat'];
jabnn = [jabnn,jabext];
classifierfile = fullfile(bdir,classifiernn);
isclassifierfile = exist(classifierfile,'file') && ~isempty(who('-file',classifierfile,'classifier'));

if ~isclassifierinjab && isclassifierfile,
  useclassifierfile = true;
  P = load(classifierfile);
  classifier = P.classifier;
  return;
end

if ~isclassifierfile,
  useclassifierfile = false;
  classifier = [];
  return;
end

classifierinfo = dir(classifierfile);
if x.classifierStuff.timeStamp < classifierinfo.datenum,
  useclassifierfile = true;
  P = load(classifierfile);
  classifier = P.classifier;
else
  useclassifierfile = false;
  classifier = [];
end
