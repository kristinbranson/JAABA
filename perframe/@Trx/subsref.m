function [varargout] = subsref(obj,s)

n = numel(s);
if n == 0,
  varargout{1} = obj;
  return;
end

% which flies are referenced?
i = 1;
switch s(i).type,
  case '()',

    switch numel(s(i).subs),
      
      % one argument means fly index
      case 1,
        if ischar(s(i).subs{1}) && strcmp(s(i).subs{1},':'),
          idx = 1:obj.nflies;
        elseif isnumeric(s(i).subs{1}),
          idx = s(i).subs{1};
          if any(idx < 1 | idx > obj.nflies | idx ~= round(idx)),
            error('Indices must be whole numbers between 1 and %d',obj.nflies);
          end
        else
          error('Bad indexing argument');
        end
        
      % two arguments means experiment + fly
      case 2,
        
        % experiment indices
        j = 1;
        if ischar(s(i).subs{j}) && strcmp(s(i).subs{j},':'),
          exps = 1:obj.nexpdirs;
        elseif ischar(s(i).subs{j}) && ismember(s(i).subs{j},obj.expdirs),
          exps = find(strcmp(s(i).subs{j},obj.expdirs),1);
        elseif iscell(s(i).subs{j}),
          [didfind,exps] = ismember(s(i).subs{j},obj.expdirs);
          if any(~didfind),
            error(['Bad experiment indexing arguments: ',sprintf('%s ',s(i).subs{j}{~didfind})]);
          end
        elseif isnumeric(s(i).subs{j}),
          exps = s(i).subs{j};
          if any(exps < 1 | exps > obj.nexpdirs | exps ~= round(exps)),
            error('Experiment indices must be whole numbers between 1 and nexpdirs = %d',obj.nexpdirs);
          end
        else
          error('Bad experiment indexing argument');
        end
        
        % fly indices
        j = 2;
        if ischar(s(i).subs{j}) && strcmp(s(i).subs{j},':'),
          idx = [obj.exp2flies{exps}];
        elseif isnumeric(s(i).subs{j}),
          idx = cell2mat(cellfun(@(x)x(s(i).subs{j}),obj.exp2flies(exps),'UniformOutput',false));
        else
          error('Bad fly indexing argument');
        end
        
      otherwise
        error('Number of indexing arguments must be <= 2');
    end

    % go on to next argument
    i = 2;
    
  case '.',
    
    if ismember(s(i).subs,methods(obj)) || ismember(s(i).subs,properties(obj)),
      [varargout{1:nargout}] = builtin('subsref',obj,s);
      return;
    end
    
    idx = 1:obj.nflies;
    
    % use this argument
    i = 1;
    
  otherwise
    error('{} referencing not supported');
end

if n < i,
  % TODO: maybe implement this?
  error('Field name must be specified');
end

% get field
if ~strcmpi(s(i).type,'.'),
  error('Field name must be specified');
end
fn = s(i).subs;
i = i + 1;

% get indices
if n < i,
  t = ':';
else
  if ~strcmp(s(i).type,'()'),
    error('Frame indices must be indexed by ()');
  end
  if numel(s(i).subs) > 1,
    error('Only one frame index allowed');
  end
  t = s(i).subs{1};
end

if strcmp(fn,'nframes'),
  res = {obj.nframes(idx)};
elseif strcmp(fn,'firstframe'),
  res = {obj.firstframes(idx)};
elseif strcmp(fn,'endframe'),
  res = {obj.endframes(idx)};
elseif strcmp(fn,'off'),
  res = {1-obj.firstframes(idx)};
else
  res = cell(1,numel(idx));
  for j = 1:numel(idx),
    x = obj.GetPerFrameData(fn,idx(j));
    if ischar(t),
      res{j} = x;
    else
      res{j} = x(t);
    end
  end
end
if nargout < numel(res),
  if numel(res) == 1,
    varargout{1:nargout} = res{1};
  else
    varargout{1:nargout} = res;
  end
else
  varargout = res;
end

