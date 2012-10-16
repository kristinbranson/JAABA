function JAABAScores2CSV(infile,outfile)

if ~exist(infile,'file'),
  error('Scores file %s does not exist',infile);
end

Q = load(infile);

fid = fopen(outfile,'w');
if ~fid,
  error('Cannot open %s for writing',outfile);
end

if isfield(Q.allScores,'postprocessed')
  ispp = true;
else
  ispp = false;
end
  
fprintf(fid,'Fly Number, Frame Number, Scores');
if ispp
  fprintf(fid,',Postprocessed Prediction\n');
else
  fprintf(fid,',Prediction\n');
end

for fly = 1:numel(Q.allScores.scores)
  for fno = Q.allScores.tStart(fly):Q.allScores.tEnd(fly)
    if ispp
      pred = sign(Q.allScores.postprocessed{fly}(fno)-0.5);
    else
      pred = sign(Q.allScores.scores{fly}(fno));
    end
    
    fprintf(fid,'%d,%d,%.4f,%d\n',fly,fno,Q.allScores.scores{fly}(fno),pred);
    
  end
  
end

fclose(fid);