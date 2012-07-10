function params = WindowFeatureName2Params(wfscurr)

STATIDX = 3;
TRANSIDX = 5;
RADIUSIDX = 7;
OFFSETIDX = 9;
RESTIDX = 10;

restdict.change = {'change_window_radius','change_window_radii'};

stats = cellfun(@(x)x{STATIDX},wfscurr,'UniformOutput',false);
trans = cellfun(@(x)x{TRANSIDX},wfscurr,'UniformOutput',false);
radii = cellfun(@(x)x{RADIUSIDX},wfscurr);
offsets = cellfun(@(x)x{OFFSETIDX},wfscurr);
rest = cellfun(@(x)x(RESTIDX:end),wfscurr,'UniformOutput',false);

[stats,~,statidx] = unique(stats);

params = {};

for i = 1:numel(stats),
  idxcurr = statidx == i;
  transcurr = unique(trans(idxcurr));
  windowscurr = unique([radii(idxcurr)',offsets(idxcurr)'],'rows');
  restcurr = cat(1,rest{idxcurr});
  paramscurr = {'windows',windowscurr,'trans_types',transcurr};
  for j = 1:2:size(restcurr,2),
    if isfield(restdict,stats{i}),
      k = find(strcmp(restcurr{i,j},restdict.(stats{i})(:,1)));
      if ~isempty(k),
        restcurr{1,j} = restdict.(stats{i}){k,2};
      end
    end
    paramscurr{end+1} = restcurr{1,j};
    if strcmp(stats{i},'harmonic') && strcmp(restcurr{1,j},'num_harmonic'),
      paramscurr{end+1} = max([restcurr{:,j+1}]);
    else
      paramscurr{end+1} = unique([restcurr{:,j+1}]);
    end
  end
  params(end+1:end+2) = {[stats{i},'_params'],paramscurr};
end
params(end+1:end+2) = {'feature_types',stats};
