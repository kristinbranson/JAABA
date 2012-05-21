function [trx,alltimestamps] = MergeMWTData(varargin)

MINDT = .0001;

trxall = varargin(1:3:end-2);
fnsall = varargin(2:3:end-1);
fnscopyall = varargin(3:3:end);
nsets = numel(trxall);
if nsets < 1,
  error('Usage: MergeMWTData(trx1,fns1,trx2,fns2,...)');
end

ids1 = unique([trxall{1}.id]);
ntargets = numel(ids1);
for seti = 2:nsets,
  ids2 = [trxall{seti}.id];
  if ~isempty(setxor(ids1,ids2)),
    error('ids do not match for sets 1 and %d',seti);
  end
end
idsall = cell(1,nsets);
for seti = 1:nsets,
  idsall{seti} = [trxall{seti}.id];
end

idx = cell(1,nsets);
for seti = 1:nsets,
  [~,idx{seti}] = ismember(ids1,idsall{seti});
end

alltimestamps = [];
for seti = 1:nsets,
  alltimestamps = union(alltimestamps,[trxall{seti}.timestamp]);
end
dts = diff(alltimestamps);
if any(dts < MINDT),
  error('Rounding error check failed: there are timestamps < %f apart',MINDT);
end


% frames per target
framesall = cell(ntargets,nsets);
firstframes = nan(ntargets,nsets);
endframes = nan(ntargets,nsets);
for seti = 1:nsets,
  for i = 1:ntargets,
    j = idx{seti}(i);
    [~,framesall{i,seti}] = ismember(trxall{seti}(j).timestamp,alltimestamps);
    firstframes(i,seti) = min(framesall{i,seti});
    endframes(i,seti) = max(framesall{i,seti});
  end
end

trx = [];
for i = 1:ntargets,
  trk = struct;
  trk.id = ids1(i);
  % union all frames for which there is data
  trk.firstframe = min(firstframes(i,:));
  trk.endframe = max(endframes(i,:));
  trk.nframes = trk.endframe - trk.firstframe + 1;
  trk.off = 1 - trk.firstframe;
  trk.timestamp = alltimestamps(trk.firstframe:trk.endframe);
  trk.dt = diff(trk.timestamp);
  for seti = 1:nsets,
    j = idx{seti}(i);
    frames = framesall{i,seti};
    for k = 1:numel(fnsall{seti}),
      fn = fnsall{seti}{k};
      sz = size(trxall{seti}(j).(fn));
      if numel(frames) ~= sz(2),
        error('Number of frames does not match size of %s for input %d',fn,seti);
      end
      if sz(2) ~= trk.nframes,
        fprintf('Adding %d empty frames to id %d, input %d, field %s\n',...
          trk.nframes-sz(2),ids1(i),fn);
      end
      sz(2) = trk.nframes;
      if iscell(trxall{seti}(j).(fn)),
        tmp = cell(sz);
      else
        tmp = nan(sz);
      end
      tmp(:,frames-trk.firstframe+1,:) = trxall{seti}(j).(fn);
      if isfield(trk,fn),
        if iscell(trk.(fn)) || iscell(tmp) || ...
            ndim(trk.(fn)) ~= numel(sz) || ...
            ~all(size(trk.(fn))==sz),
          warning('Field %s in multiple inputs and cannot be merged, skipping input %d',fn,seti);
        else
          idxdouble = ~isnan(trk.(fn)) & ~isnan(tmp);
          if any(trk.(fn)(idxdouble) ~= tmp(idxdouble)),
            warning('Field %s is in multiple inputs and does not match, skipping mismatches from input %d',fn,seti);
          end
          idxadd = isnan(trk.(fn)) & ~isnan(tmp);
          trk.(fn)(idxadd) = tmp(idxadd);
        end
      else
        trk.(fn) = tmp;
      end
    end
    for k = 1:numel(fnscopyall{seti}),
      fn = fnscopyall{seti}{k};
      trk.(fn) = trxall{seti}.(fn);
    end
  end
  trx = structappend(trx,trk);
end

% check for fields that got lost
fnsmerge = fieldnames(trx);
for seti = 1:nsets,
  fns = fieldnames(trxall{seti});
  missingfns = setdiff(fns,fnsmerge);
  if ~isempty(missingfns),
    fprintf('The following fields from input %d were not merged:',seti);
    fprintf(' %s',missingfns{:});
    fprintf('\n');
  end
end