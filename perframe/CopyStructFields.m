function a = CopyStructFields(b,a)

fns = fieldnames(b);
for i = 1:numel(fns),
  fn = fns{i};
  if isfield(a,fn) && isstruct(a.(fn)) && isstruct(b.(fn)),
    a.(fn) = CopyStructFields(b.(fn),a.(fn));
  else
    if isfield(a,fn) && isstruct(a.(fn)) ~= isstruct(b.(fn)),
      warning('CopyStructFields: isstruct does not match for %s, overwriting',fn);
    end
    a.(fn) = b.(fn);
  end
end