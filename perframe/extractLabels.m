function extractLabels(jabfile,outfile)

J = load(jabfile,'-mat');
oo = fopen(outfile,'w');

for exp = 1:numel(J.x.labels)
  for ndx = 1:numel(J.x.labels(exp).t0s)
    for bndx = 1:numel(J.x.labels(exp).t0s{ndx})
      fprintf(oo,'%s %d %d %d %s\n',J.x.expDirNames{exp},...
        J.x.labels(exp).flies(ndx),...
        J.x.labels(exp).t0s{ndx}(bndx),...
        J.x.labels(exp).t1s{ndx}(bndx)-1,...
        J.x.labels(exp).names{ndx}{bndx}...
        );
    end
    
  end
  
end
fclose(oo);