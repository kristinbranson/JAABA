function [classifier,useclassifierfile,fixjabfile] = getClassifierFromJabFile(jabfile,varargin)

[useclassifierfile] = ...
  myparse(varargin,...
  'useclassifierfile',[]);

if ~exist(jabfile,'file'),
  error('jabfile %s does not exist',jabfile);
end

% check if the classifier file exists and is newer
[bdir,jabnn,jabext] = fileparts(jabfile);
classifiernn = [jabnn '_classifier.mat'];
jabnn = [jabnn,jabext];
classifierfile = fullfile(bdir,classifiernn);

% load in jab file
x = loadAnonymous(jabfile);

if isempty(useclassifierfile),
  
  if exist(classifierfile,'file') && ~isempty(who('-file',classifierfile,'classifier')),
    assert(~x.isMultiClassifier(),...
      'Unexpected external classifier file for multi-classifier jabfile.');
    
    jabinfo = dir(jabfile);
    classifierinfo = dir(classifierfile);
  
    % check to see if the classifier is set in the jab file
    isclassifierinjab = ~isempty(x.classifierStuff.params);
    
    if isclassifierinjab && jabinfo.datenum >= classifierinfo.datenum,
      useclassifierfile = false;
      fixjabfile = false;
    else      
    
      res = questdlg(sprintf('For jab file %s (last saved %s) classifier file %s (last saved %s) is newer. ',...
        jabnn,jabinfo.date,classifiernn,classifierinfo.date),...
        'Use classifier from jab file or classifier.mat file?',...
        'Fix jab file','Use jab file','Use classifier file, do not fix jab file',...
        'Fix jab file');
      
      if strcmpi(res,'Fix jab file'),
        fixjabfile = true;
        useclassifierfile = true;
      elseif strcmpi(res,'Use classifier file, do not fix jab file'),
        fixjabfile = false;
        useclassifierfile = true;
      else
        fixjabfile = false;
        useclassifierfile = false;
      end
    end
  else
    useclassifierfile = false;
    fixjabfile = false;
  end
else
  fixjabfile = useclassifierfile;
end

if useclassifierfile,      
  P = load(classifierfile);
  classifier = P.classifier;
  if fixjabfile,
    assert(~x.isMultiClassifier(),...
      'Unexpected external classifier file for multi-classifier jabfile.');
    try
      x.classifierStuff.params = classifier;
      [success,msg] = movefile(jabfile,[jabfile,'.old']);
      if ~success,
        warning('Could not rename %s to %s: %s',jabfile,[jabfile,'.old'],msg);
      end
      save(jabfile,'x');
    catch ME,
      warning('Could not fix jab file %s: %s',jabfile,getReport(ME));
    end
  end
else
  classifier = x.classifierStuff.params;
end

