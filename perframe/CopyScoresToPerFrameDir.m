function CopyScoresToPerFrameDir(expdir,varargin)

[scoresfilestrs,perframedirstr,doforce] = myparse(varargin,'scoresfilestrs',{},...
  'perframedirstr','perframe','doforce',false);

if isempty(scoresfilestrs),
  scoresfilestrs = dir(fullfile(expdir,'scores*mat'));
  scoresfilestrs = {scoresfilestrs.name};
end

for i = 1:numel(scoresfilestrs),
  infile = fullfile(expdir,scoresfilestrs{i});
  outfile = fullfile(expdir,perframedirstr,scoresfilestrs{i});
  if ~exist(infile,'file'),
    fprintf('File %s does not exist, skipping\n',infile);
    continue;
  end
  if exist(outfile,'file') && ~doforce,
    res = questdlg(sprintf('File %s exists, overwrite?',outfile));
    if strcmpi(res,'No'),
      continue;
    elseif strcmpi(res,'Cancel'),
      break;
    end
  end
  scores = load(infile);
  nflies = numel(scores.allScores.scores);
  perframe = struct;
  perframe.data = cell(1,nflies);
  for j = 1:nflies,
    perframe.data{j} = scores.allScores.scores{j}(scores.allScores.tStart(j):scores.allScores.tEnd(j));
  end
  perframe.units = parseunits('unit');
  save(outfile,'-struct','perframe');
  
  %  copyfile(infile,outfile);
end