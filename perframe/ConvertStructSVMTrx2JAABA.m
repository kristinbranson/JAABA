function ConvertStructSVMTrx2JAABA(trxfiles,labelfiles,expdir,varargin)

[trxfilestr,moviefilestr,perframedirstr] = myparse(varargin,...
  'trxfilestr','trx.mat','moviefilestr','movie.ufmf',...
  'perframedirstr','perframe');

if numel(trxfiles) ~= numel(labelfiles),
  error('Number of trx and label files do not match');
end

trx = [];
alldata = {};
allfieldnames = {};
allunits = {};
alllabels = [];
allbehaviors = {};

% loop through trxfiles and labelfiles, one for each fly
for i = 1:numel(trxfiles),
  
  trxfile = trxfiles{i};
  labelfile = labelfiles{i};
  
  % read the trx file
  [trxcurr,fieldnamescurr,datacurr,unitscurr] = ReadStructSVMTrxFile(trxfile);
  trx = structappend(trx,trxcurr);
  id = trxcurr.id;
  
  % check that it is from the same movie
  if i == 1,
    moviefile = trxcurr.moviename;
  else
    if ~strcmp(moviefile,trxcurr.moviename),
      error('Movies mismatch');
    end
  end

  % add on the perframe data
  [ism,idx] = ismember(fieldnamescurr,allfieldnames);
  if any(ism),
    for j = find(ism),
      alldata{idx(j)}{id} = datacurr{j};
    end
  end
  if any(~ism),
    for j = find(~ism),
      alldata{end+1} = {}; %#ok<*AGROW>
      alldata{end}{id} = datacurr{j};
      allunits{end+1} = unitscurr{j};
      allfieldnames{end+1} = fieldnamescurr{j};
    end
  end
  
  % read the labels    
  [labelscurr,behaviorscurr] = ReadStructSVMLabelFile(labelfile);
  
  % add them to alllabels
  [ism,idx] = ismember(behaviorscurr,allbehaviors);
  if any(ism),
    for j = find(ism),
      alllabels{idx(j)}.t0s{end+1} = labelscurr{j}.t0s;
      alllabels{idx(j)}.t1s{end+1} = labelscurr{j}.t1s;
      alllabels{idx(j)}.names{end+1} = labelscurr{j}.names;
      alllabels{idx(j)}.flies(:,end+1) = labelscurr{j}.flies;
      alllabels{idx(j)}.off(end+1) = labelscurr{j}.off;
      alllabels{idx(j)}.timestamp{end+1} = labelscurr{j}.timestamp;
      alllabels{idx(j)}.imp_t0s{end+1} = labelscurr{j}.imp_t0s;
      alllabels{idx(j)}.imp_t1s{end+1} = labelscurr{j}.imp_t1s;
    end
  end
  if any(~ism),
    for j = find(~ism),
      alllabels{end+1} = labelscurr{j};
      alllabels{end}.t0s = {alllabels{end}.t0s};
      alllabels{end}.t1s = {alllabels{end}.t1s};
      alllabels{end}.names = {alllabels{end}.names};
      alllabels{end}.timestamp = {alllabels{end}.timestamp};
      alllabels{end}.imp_t0s = {alllabels{end}.imp_t0s};
      alllabels{end}.imp_t1s = {alllabels{end}.imp_t1s};
      allbehaviors{end+1} = behaviorscurr{j};
    end
  end  

end

if ~exist(expdir,'dir'),
  mkdir(expdir);
end

% sort trx by id
[~,order] = sort([trx.id]);
trx = trx(order);

% create timestamps
timestamps = nan(1,max([trx.endframe]));
for i = 1:numel(trx),
  timestamps(trx(i).firstframe:trx(i).endframe) = trx(i).timestamps;
end

% save trx
trxfile = fullfile(expdir,trxfilestr);
save(trxfile,'trx','timestamps');

% create perframe dir
perframedir = fullfile(expdir,perframedirstr);
if ~exist(perframedir,'dir'),
  mkdir(perframedir);
end

% save perframe data
for i = 1:numel(alldata),
  data = alldata{i}; %#ok<NASGU>
  fn = allfieldnames{i};
  units = allunits{i}; %#ok<NASGU>
  filename = fullfile(perframedir,[fn,'.mat']);
  save(filename,'data','units');
end

% save labels
for i = 1:numel(alllabels),
  filestr = sprintf('labeled_%s.mat',allbehaviors{i});
  labels = alllabels{i}; %#ok<NASGU>
  save(fullfile(expdir,filestr),'-struct','labels');
end
  

if exist(moviefile,'file'),
  copyfile(moviefile,fullfile(expdir,moviefilestr));
end