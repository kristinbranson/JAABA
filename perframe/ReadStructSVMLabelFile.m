function [labels,behaviors] = ReadStructSVMLabelFile(labelfile)

fid = fopen(labelfile,'rb');
version = fread(fid,1,'double'); %#ok<NASGU>
moviename_length = fread(fid,1,'int');
moviename = fread(fid,moviename_length,'char=>char')'; %#ok<NASGU>
matname_length = fread(fid,1,'int');
matname = fread(fid,matname_length,'char=>char')'; %#ok<NASGU>
trxname_length = fread(fid,1,'int');
trxname = fread(fid,trxname_length,'char=>char')'; %#ok<NASGU>
nflies = fread(fid,1,'int'); %#ok<NASGU>
fly_ids = fread(fid,1,'int');
firstframe = fread(fid,1,'int');
firstframe = firstframe + 1;
lastframe = fread(fid,1,'int');
lastframe = lastframe + 1;
nbehaviors = fread(fid,1,'int');
nframes = lastframe-firstframe+1;

behaviornames = cell(1,nbehaviors);

for i = 1:nbehaviors,
  behaviorname_length = fread(fid,1,'int');
  behaviorname = fread(fid,behaviorname_length,'char=>char')';
  behaviornames{i} = behaviorname;
end
nbouts = fread(fid,1,'int');
data = fread(fid,nbouts*3,'int');
fclose(fid);

% assuming startframes and endframes are relative -- CHECK THIS
startframes = data(1:3:end-2);
endframes = data(2:3:end-1);
behavioridx = data(3:3:end);
startframes = startframes+1;
endframes = endframes+1;
behavioridx = behavioridx+1;

% initialize everything as None
labelidx0 = repmat(2,[1,nframes]);
isunknown = ismember(behaviornames,{'*Unknown*','*Tracker Failure*'});
idx_behavior = find(~isunknown);
idx = find(isunknown(behavioridx));
% set unknown labels
for i = idx(:)',
  labelidx0(startframes(i):endframes(i)-1) = 0;
end

labels = cell(1,numel(idx_behavior));
behaviors = cell(1,numel(idx_behavior));
timestamp = now;
for ii = 1:numel(idx_behavior),
  i = idx_behavior(ii);
  labelidx = labelidx0;
  idx = find(behavioridx == i);
  for j = idx(:)',
    labelidx(startframes(j):endframes(j)-1) = 1;
  end
  
  behaviors_curr = {behaviornames{i},'None'};
  labels_curr = struct;
  labels_curr.t0s = [];
  labels_curr.t1s = [];
  labels_curr.names = {};
  labels_curr.timestamp = [];
  labels_curr.flies = fly_ids;
  labels_curr.imp_t0s = [];
  labels_curr.imp_t1s = [];
  labels_curr.off = 1-firstframe;

  for j = 1:2,
    [i0s,i1s] = get_interval_ends(labelidx==j);        
    if ~isempty(i0s),
      n = numel(i0s);
      labels_curr.t0s(end+1:end+n) = i0s + firstframe - 1;
      labels_curr.t1s(end+1:end+n) = i1s + firstframe - 1;
      labels_curr.imp_t0s(end+1:end+n) = i0s + firstframe - 1;
      labels_curr.imp_t1s(end+1:end+n) = i1s + firstframe - 1;

      labels_curr.names(end+1:end+n) = repmat(behaviors_curr(j),[1,n]);
      labels_curr.timestamp(end+1:end+n) = timestamp;
    end
  end

  labels{ii} = labels_curr;
  behaviors{ii} = behaviornames{i};
  
end

