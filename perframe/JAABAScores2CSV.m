function JAABAScores2CSV(scorefile,outfile)

% function JAABAScores2CSV(scorefile,outfile)
% Convert JAABA Scores to comma separated values.
% The format of CSV file is Animal Num, Frame Num, Score Value
% (unnormalized).


if ~exist(scorefile,'file'),
  error('Scorefile %s doesn''t exist\n',scorefile);
end
Q = load(scorefile);

fid = fopen(outfile,'w');

for ndx = 1:numel(Q.allScores.scores)
  for t = Q.allScores.tStart(ndx):Q.allScores.tEnd(ndx)
    fprintf(fid,'%d,%d,%.4f\n',ndx,t,Q.allScores.scores{ndx}(t));
  end
end
fclose(fid);
