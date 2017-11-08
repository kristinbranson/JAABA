%% set up paths

addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/misc;
addpath /groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/filehandling;

datafiles = {'m119nocnoraw.csv'
  'm119cnoraw.csv'};

timestamp = now;

%% read in data files

rawdata = [];
for i = 1:numel(datafiles),
  
  rawdatacurr = ReadRawDataFile(datafiles{i});
  rawdata = structappend(rawdata,rawdatacurr);
  
end

%% 

ntrials = numel(rawdata);

fnsspecial = {'auto_Grab_success','auto_Grab_successtype'};
fns = fieldnames(rawdata);
datafns = setdiff(fns(~cellfun(@isempty,regexp(fns,'^auto'))),fnsspecial);
nstats = numel(datafns);

% some behaviors don't have counts
tmp = regexp(datafns,'^auto_([^_]+)_0$','tokens','once');
behaviors = unique([tmp{:}]);
nbehaviors = numel(behaviors);
for i = 1:nbehaviors,
  j = 0;
  fnprev = '';
  fncount = sprintf('auto_%s_num',behaviors{i});
  if ismember(fncount,datafns),
    continue;
  end
  count = zeros(ntrials,1);
  while true,
    fn = sprintf('auto_%s_%d',behaviors{i},j);
    if ~ismember(fn,datafns),
      break;
    end
    tmp = ~isnan([rawdata.(fn)]);
    if j > 0,
      tmp = tmp & [rawdata.(fn)]~=[rawdata.(fnprev)];
    end
    count(tmp) = count(tmp) + 1;
    j = j+1;
    fnprev = fn;
  end
  for j = 1:ntrials,
    rawdata(j).(fncount) = count(j); %#ok<SAGROW>
  end
end

datafns = setdiff(fns(~cellfun(@isempty,regexp(fns,'^auto'))),fnsspecial);
nstats = numel(datafns);

% set to nan when del = 0
for i = find(~cellfun(@isempty,regexp(datafns,'^auto_del','once'))),
  fn = datafns{i};
  for j = 1:ntrials,
    if rawdata(j).(fn) == 0,
      rawdata(j).(fn) = nan;
    end
  end
end

X = nan(ntrials,nstats);
for i = 1:nstats,
  X(:,i) = [rawdata.(datafns{i})];
end

% find redundant stats
d = inf(nstats,nstats);
n = nan(nstats,1);
idxreplace = zeros(1,nstats);
for i = 1:nstats,
  idx = ~isnan(X(:,i));
  n(i) = nnz(idx);
  for j = 1:nstats,
    if i == j,
      continue;
    end
    d(i,j) = nnz(X(idx,j) ~= X(idx,i)) / nnz(idx);
    if idxreplace(j) > 0,
      continue;
    end
    if d(i,j) == 0,
      idxreplace(i) = j;
      idxreplace(idxreplace==i) = j;
      fprintf('Can compute %d: %s from %d: %s\n',i,datafns{i},idxreplace(i),datafns{idxreplace(i)});
    end
  end
end

% is there any information in the nan pattern?
idx = ~cellfun(@isempty,regexp(datafns,'num$','once'));
countfns = datafns(idx);

Z = nan(ntrials,numel(countfns)+2);
for i = 1:numel(countfns),
  Z(:,i) = [rawdata.(countfns{i})];
end
Z(:,end-1) = double([rawdata.success]);
[successtypes,~,Z(:,end)] = unique({rawdata.successtype});

mindiff = inf(1,nstats);
idxnanpattern = zeros(6,nstats);
for i = find(idxreplace>0),
  isdata = ~isnan(X(:,i));
  
  for j1 = 1:size(Z,2),
    
    for j2 = 1:size(Z,2),
      
      for k1 = min(Z(:,j1)):max(Z(:,j1)),

        if j1 == size(Z,2),
          isdata1 = Z(:,j1) == k1;
        else
          isdata1 = (Z(:,j1) >= k1);
        end
        
        for k2 = min(Z(:,j2)):max(Z(:,j2)),

          if j2 == size(Z,2),
            isdata2 = Z(:,j2) == k2;
          else
            isdata2 = (Z(:,j2) >= k2);
          end

          for l1 = 1:2,
            
            if l1 == 1,
              isdata1x = isdata1;
            else
              isdata1x = ~isdata1;
            end
            for l2 = 1:2,
              if l2 == 1,
                isdata2x = isdata2;
              else
                isdata2x = ~isdata2;
              end

              isdata3 = isdata1x & isdata2x;
              
              nmismatch = nnz(isdata ~= isdata3);
              if nmismatch < mindiff(i),
                idxnanpattern(:,i) = [j1,j2,k1,k2,l1,l2];
                mindiff(i) = nmismatch;
              end
            end
          end
        end
      end
    end
  end
end

fprintf('Max info in nan patterns = %d / %d\n',max(mindiff(idxreplace>0)),ntrials);

datafns = datafns(idxreplace==0);
nstats = numel(datafns);

X = nan(ntrials,nstats);
for i = 1:nstats,
  X(:,i) = [rawdata.(datafns{i})];
end

%% compare features and various statistics

% success or not
y = [rawdata.success];

hfig = CompareStatsAcrossGroups(X,y,datafns,{'Failure','Success'},'title','Success vs failure, all trials');

savefig(sprintf('SortedStatisticsPredictingSuccess_AllTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsPredictingSuccess_AllTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

% cno vs not cno
y = [rawdata.iscno];
hfig = CompareStatsAcrossGroups(X,y,datafns,{'No CNO','CNO'},'title','No CNO vs CNO, all trials');

savefig(sprintf('SortedStatisticsPredictingCNO_AllTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsPredictingCNO_AllTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

% success CNO vs success not CNO
idx = [rawdata.success] == 1;
y = [rawdata.iscno];
hfig = CompareStatsAcrossGroups(X(idx,:),y(idx),datafns,{'No CNO','CNO'},'title','No CNO vs CNO, success trials');

savefig(sprintf('SortedStatisticsPredictingCNO_SuccessTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsPredictingCNO_SuccessTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

% failure CNO vs failure not CNO
idx = [rawdata.success] == 0;
y = [rawdata.iscno];
hfig = CompareStatsAcrossGroups(X(idx,:),y(idx),datafns,{'No CNO','CNO'},'title','No CNO vs CNO, failure trials');

savefig(sprintf('SortedStatisticsPredictingCNO_FailureTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsPredictingCNO_FailureTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

% correlate trial to statistics
X1 = [X,double([rawdata.success])'];
datafns1 = [datafns,{'success'}];
y = [rawdata.trial];

hfig = MeasureStatDependencies(X1,y,datafns1,'Effect of trial number, all trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_AllTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_AllTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

idx = [rawdata.success]==1;
hfig = MeasureStatDependencies(X(idx,:),y(idx),datafns,'Effect of trial number, success trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_SuccessTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_SuccessTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

idx = [rawdata.success]==0;
hfig = MeasureStatDependencies(X(idx,:),y(idx),datafns,'Effect of trial number, failure trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_FailureTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_FailureTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

idx = [rawdata.iscno]==1;
hfig = MeasureStatDependencies(X1(idx,:),y(idx),datafns1,'Effect of trial number, CNO trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_CNOTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_CNOTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

idx = [rawdata.iscno]==0;
hfig = MeasureStatDependencies(X1(idx,:),y(idx),datafns1,'Effect of trial number, no CNO trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_noCNOTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_noCNOTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

%% how much does a single date contribute?

% correlate trial to statistics
X1 = [X,double([rawdata.success])'];
datafns1 = [datafns,{'success'}];
y = [rawdata.trial];

[uniquedates,~,dateidx] = unique({rawdata.date});

idx = [rawdata.iscno]==1;
[tmp,~,groupidx] = unique(dateidx(idx));
groupnames = uniquedates(tmp);

hfig = MeasureStatDependenciesPerGroup(X1(idx,:),y(idx),groupidx,datafns1,groupnames,'Effect of trial number per date, CNO trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_PerGroup_CNOTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_PerGroup_CNOTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');


idx = [rawdata.iscno]==0;
[tmp,~,groupidx] = unique(dateidx(idx));
groupnames = uniquedates(tmp);

hfig = MeasureStatDependenciesPerGroup(X1(idx,:),y(idx),groupidx,datafns1,groupnames,'Effect of trial number per date, no CNO trials');

savefig(sprintf('SortedStatisticsVsTrialNumber_PerGroup_noCNOTrials%s.pdf',datestr(timestamp,'yyyymmdd')),hfig,'pdf');
savefig(sprintf('SortedStatisticsVsTrialNumber_PerGroup_noCNOTrials%s.png',datestr(timestamp,'yyyymmdd')),hfig,'png');

%% make a table showing all statistics, all groups

ys = nan(0,numel(rawdata));
idxs = false(0,numel(rawdata));
ynames = {};
idxnames = {};

ys(end+1,:) = [rawdata.success];
ynames{end+1} = 'Success';
ys(end+1,:) = [rawdata.iscno];
ynames{end+1} = 'CNO';

idxs(end+1,:) = true(1,numel(rawdata));
idxnames{end+1} = 'All trials';
idxs(end+1,:) = [rawdata.success];
idxnames{end+1} = 'Success trials';
idxs(end+1,:) = ~[rawdata.success];
idxnames{end+1} = 'Failure trials';
idxs(end+1,:) = [rawdata.iscno];
idxnames{end+1} = 'CNO trials';
idxs(end+1,:) = ~[rawdata.iscno];
idxnames{end+1} = 'No CNO trials';

pairsignore = {'Success','Success trials'
  'Success','Failure trials'
  'CNO','CNO trials'
  'CNO','No CNO trials'};

hfig = CompareStatsAcrossMultipleGroups(X,ys,ynames,idxs,idxnames,datafns,pairsignore);
