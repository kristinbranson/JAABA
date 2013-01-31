function checkTransTypes(clfile)
Q = load(clfile);
params = Q.windowfeaturesparams;
[params ok]= convertTransTypes2Cell(params);
if ok,
  cellparams = convertParams2CellParams(params);
  if ~isequal(cellparams,Q.windowfeaturescellparams),
    fprintf('Cell Params don'' match\n');
  end
end

function [params ok]= convertTransTypes2Cell(params)
  % Convert the trans_types field into cell type
  ok = true;
  if ~isstruct(params),  return; end
  fnames = fieldnames(params);
  for ndx = 1:numel(fnames)
    if isstruct(params.(fnames{ndx}))
      [params.(fnames{ndx}) ok1]= convertTransTypes2Cell(params.(fnames{ndx}));
      ok = ok & ok1;
    end
  end
  if isfield(params,'trans_types')&& ~iscell(params.trans_types)
    fprintf('Not all params had trans_types as cell\n');
    ok = false;
    params.trans_types = {params.trans_types};
  end

  function cellparams = convertParams2CellParams(params)
    cellparams = struct;
    fns1 = fieldnames(params);
    for i1 = 1:numel(fns1),
      fn1 = fns1{i1};
      fns2 = fieldnames(params.(fn1));
      cellparams.(fn1) = {};
      feature_types = {};
      for i2 = 1:numel(fns2),
        fn2 = fns2{i2};
        if ~isstruct(params.(fn1).(fn2)),
          cellparams.(fn1)(end+1:end+2) = {fn2,params.(fn1).(fn2)};
        else
          cellparams.(fn1)(end+1:end+2) = {[fn2,'_params'],struct2paramscell(params.(fn1).(fn2))};
          feature_types{end+1} = fn2; %#ok<AGROW>
        end
      end
      cellparams.(fn1)(end+1:end+2) = {'feature_types',feature_types};
    end
    

