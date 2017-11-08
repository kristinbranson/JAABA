function [outjabfile,outclassifierfile,injabfile,jabdata] = ExportClassifier(varargin)

jabdata = [];

persistent lastinjabfile;
if isempty(lastinjabfile),
  lastinjabfile = '';
end

[injabfile,outjabfile] = myparse(varargin,'injabfile','','outjabfile','');

if isempty(injabfile),
  
  [f,p] = uigetfile('*.jab','Select input jab file',lastinjabfile);
  if ~ischar(f),
    return;
  end
  
  injabfile = fullfile(p,f);
    
end

if ~exist(injabfile,'file'),
  error('Input jab file %s does not exist',injabfile);
end

% find corresponding classifier file
[p,n] = fileparts(injabfile);
inclassifierfile = fullfile(p,[n,'_classifier.mat']);

if ~exist(inclassifierfile,'file'),
  error('Classifier file %s corresponding to jab file %s does not exist',inclassifierfile,injabfile);
end

lastinjabfile = injabfile;
jabdata = load(injabfile,'-mat');

% remove experiments:

% labels
jabdata.x.labels = jabdata.x.labels(zeros(1,0));
jabdata.x.gtLabels = jabdata.x.gtLabels(zeros(1,0));

% experiment names
jabdata.x.expDirNames = jabdata.x.expDirNames(zeros(1,0));
jabdata.x.gtExpDirNames = jabdata.x.gtExpDirNames(zeros(1,0));

if isempty(outjabfile),

  inp = fileparts(injabfile);
  
  [f,p] = uiputfile('*.jab','Select output jab file to export the classifier to',inp);
  if ~ischar(f),
    return;
  end

  outjabfile = fullfile(p,f);
  
end

  
[p,n] = fileparts(outjabfile);
outclassifierfile = fullfile(p,[n,'_classifier.mat']);
[success,msg] = copyfile(inclassifierfile,outclassifierfile);
if ~success,
  error('Error copying classifier file %s to %s: %s',inclassifierfile,outclassifierfile,msg);
end
  
save(outjabfile,'-struct','jabdata');

fprintf('Output jab file: %s\n',outjabfile);
fprintf('Output classifier file: %s\n',outclassifierfile);
