function metadata = GetExpTypeFromDir(expdir)

m = regexp(expdir,'(M\d+)[^/\\]*[/\\]([^/\\]+)[/\\](\d+)(CNO)?[/\\]([^/\\]+_v)(\d+)$','once','tokens');
if iscell(expdir),
  assert(~any(cellfun(@isempty,m)));
  mouse = cellfun(@(x) x{1},m,'Uni',false);
  session = cellfun(@(x) x{2},m,'Uni',false);
  day = cellfun(@(x) x{3},m,'Uni',false);
  iscno = ~cellfun(@(x) isempty(x{4}),m);
  trial = cellfun(@(x) str2double(x{6}),m);
  exp = cellfun(@(x) [x{5},x{6}],m,'Uni',false);
  id = cell(1,numel(expdir));
  for i = 1:numel(expdir),
    if iscno(i),
      id{i} = [mouse{i}(2:end),'#',day{i},'CNO','#',exp{i}];
    else
      id{i} = [mouse{i}(2:end),'#',day{i},'#',exp{i}];
    end
  end
  metadata = struct('mouse',mouse,'session',session,'day',day,...
    'iscno',num2cell(iscno),'trial',num2cell(trial),'exp',exp,...
    'expdir',expdir,'id',id);
else
  assert(~isempty(m));
  mouse = m{1};
  session = m{2};
  day = m{3};
  iscno = ~isempty(m{4});
  trial = str2double(m{6});
  exp = [m{5},m{6}];
  if iscno,
    id = [mouse(2:end),'#',day,'CNO#',exp];
  else
    id = [mouse(2:end),'#',day,'#',exp];
  end
  metadata = struct('mouse',mouse,'session',session,'day',day,'iscno',iscno,...
    'trial',trial,'exp',exp,'expdir',expdir,'id',id);
end