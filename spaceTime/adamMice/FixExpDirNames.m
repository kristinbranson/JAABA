function outx = FixExpDirNames(injabfile,outjabfile)

outx = [];

if nargin < 1,
  [f,p] = uigetfile('*.jab','Select jab file to convert');
  if ~ischar(f),
    return;
  end
  injabfile = fullfile(p,f);
end

Q = load(injabfile,'-mat');
persistent lastdir;

if isempty(lastdir),
  lastdir = '';
end

% find common root directories in input
inrootdirs = {};
for i = 1:numel(Q.x.expDirNames),
  s1 = regexp(Q.x.expDirNames{i},'[/\\]','split');
  for j = i+1:numel(Q.x.expDirNames),
    s2 = regexp(Q.x.expDirNames{j},'[/\\]','split');
    k = min(numel(s1),numel(s2));
    k = find(cellfun(@(x,y) strcmp(x,y),s1(1:k),s2(1:k)),1,'last');
    if ~isempty(k),
      inrootdir = fullfile(s1{1:k});
      if isempty(s1{1}),
        inrootdir = ['/',inrootdir]; %#ok<AGROW>
      end
        
      inrootdirs = union(inrootdirs,{inrootdir});
    end
  end
end

fprintf('Found %d common root directories in the input:\n',numel(inrootdirs));
fprintf('%s\n',inrootdirs{:});

conversions = cell(0,2);
isconverted = false(1,numel(Q.x.expDirNames));
outexpdirs = Q.x.expDirNames;
while true,

  if ~isempty(conversions),
    fprintf('So far, the following conversions have been created:\n');
    for i = 1:size(conversions,1),
      fprintf('IN = %s --> OUT = %s\n',conversions{i,1},conversions{i,2});
    end
    fprintf('These rules apply to %d / %d input experiment directories\n',nnz(isconverted),numel(isconverted));
  end
  fprintf('Type the input root directory you would like to replace.\n');
  fprintf('Hit enter with an empty line if you are done.\n');

  res = input('Input directory: ','s');
  if isempty(res),
    break;
  end
  inrootdir = res;

  if ~ismember(res,inrootdirs),
    res = input(sprintf('%s is not in the list of input root directories auto-discovered. Are you sure you want to use it? [Y/n]: ',inrootdir),'s');
    if ~isempty(res) && strcmpi(res(1),'n'),
      continue;
    end
  end
  
  s1 = regexp(inrootdir,'[/\\]','split');

  doesmatch = false(1,numel(Q.x.expDirNames));
  
  for i = find(~isconverted),
    s2 = regexp(Q.x.expDirNames{i},'[/\\]','split');
    doesmatch(i) = numel(s2) >= numel(s1) && all(strcmp(s1,s2(1:numel(s1))));
  end
  fprintf('Input root directory %s matches %d / %d of the remaining experiments:\n',inrootdir,nnz(doesmatch),nnz(~isconverted));
  fprintf('%s\n',Q.x.expDirNames{doesmatch});
  
  res = input('Convert all of these? [Y/n]: ','s');
  if ~isempty(res) && strcmpi(res(1),'n'),
    continue;
  end
  
  fprintf('Select output root directory corresponding to %s...\n',inrootdir);
  outrootdir = uigetdir(lastdir,inrootdir);  
  
  if ~ischar(outrootdir),
    continue;
  end
  
  lastdir = fileparts(outrootdir);

  newoutexpdirs = outexpdirs;
  for i = find(doesmatch),
    s2 = regexp(Q.x.expDirNames{i},'[/\\]','split');
    newoutexpdir = outrootdir;
    didmatch = true;
    for j = numel(s1)+1:numel(s2),
      d = dir(newoutexpdir);
      d = d([d.isdir]);
      didmatch = false;
      for k = 1:numel(d),
        didmatch = lenientMatch(d(k).name,s2{j});
        if didmatch,
          break;
        end
      end
      if ~didmatch,
        fprintf('Output directory %s does not exist, not matching\n',fullfile(outrootdir,s2{numel(s1)+1:end}));
        break;
      end
      newoutexpdir = fullfile(newoutexpdir,d(k).name);
    end
    if didmatch,
      newoutexpdirs{i} = newoutexpdir;
    else
      doesmatch(i) = false;
    end
  end
  fprintf('The following conversions will be made:\n');
  for i = find(doesmatch),
    fprintf('IN = %s -> OUT = %s\n',Q.x.expDirNames{i},newoutexpdirs{i});
  end
  res = input('Proceed with conversion? [Y/n]: ','s');
    
  if ~isempty(res) && strcmpi(res(1),'n'),
    continue;
  end

  outexpdirs = newoutexpdirs;
  isconverted = isconverted | doesmatch;

  if all(isconverted),
    break;
  end
  
end

fprintf('%d / %d directories were converted.\n',nnz(isconverted),numel(isconverted));

outx = Q.x;
outx.expDirNames = outexpdirs;

res = input('Save the new jab file with these conversions? [Y/n]: ','s');
if  ~isempty(res) && strcmpi(res(1),'n'),
  return;
end

if nargin < 2,
  [p,f,e] = fileparts(injabfile);
  outjabfile = fullfile(p,[f,'_converted',e]);
  [f,p] = uiputfile('*.jab','Choose output jab file',outjabfile);
  if ~ischar(f),
    return;
  end
  outjabfile = fullfile(p,f);
end

x = outx;
save(outjabfile,'x','-mat');
