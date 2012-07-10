function v = WindowFeatureNameCompare(f1,f2)

v = [];

if isempty(f1) || isempty(f2),
  return;
end

iscell1 = iscell(f1{1});
iscell2 = iscell(f2{1});

if iscell1 && iscell2,
  if numel(f1) ~= numel(f2),
    error('Inputs f1 and f2 do not match');
  else
    v = false(size(f1));
    for i = 1:numel(f1),
      if numel(f1{i}) == numel(f2{i}),
        vcurr = true;
        for j = 1:numel(f1{i}),
          if strcmp(class(f1{i}{j}),class(f2{i}{j})),
            if isnumeric(f1{i}{j}),
              if ~all(f1{i}{j} == f2{i}{j}),
                vcurr = false;
                break;
              end
            else
              if ~strcmp(f1{i}{j},f2{i}{j}),
                vcurr = false;
                break;
              end
            end
          end
        end
        v(i) = vcurr;
      end
    end
  end
  return;
end

if iscell2,
  v = WindowFeatureNameCompare(f2,f1);
  return;
end

if iscell1,
  v = false(size(f1));
  for i = 1:numel(f1),
    if numel(f1{i}) == numel(f2),
      vcurr = true;
      for j = 1:numel(f1{i}),
        if strcmp(class(f1{i}{j}),class(f2{j})),
          if isnumeric(f1{i}{j}),
            if ~all(f1{i}{j} == f2{j}),
              vcurr = false;
              break;
            end
          else
            if ~strcmp(f1{i}{j},f2{j}),
              vcurr = false;
              break;
            end
          end
        end
      end
      v(i) = vcurr;
    end
  end
  return;
end


v = true;
for j = 1:numel(f1),
  if strcmp(class(f1{j}),class(f2{j})),
    if isnumeric(f1{j}),
      if ~all(f1{j} == f2{j}),
        v = false;
        break;
      end
    else
      if ~strcmp(f1{j},f2),
        v = false;
        break;
      end
    end
  end
end
