function copyData(rootdir,outdir)


dd = dir(rootdir);
for ndx =1:numel(dd)
  if strcmp(dd(ndx).name(1),'.'), continue; end
  curd = fullfile(rootdir,dd(ndx).name);
  if isdir(curd)
    if ~exist(fullfile(outdir,dd(ndx).name),'dir')
      mkdir(fullfile(outdir,dd(ndx).name));
    end
      copyData(curd,fullfile(outdir,dd(ndx).name));
  else
    if ~exist(fullfile(outdir,dd(ndx).name),'file')
      infile = regexprep(curd,' ','\\ ');
      outfile = regexprep(fullfile(outdir,dd(ndx).name),' ','\\ ');
      cmd = sprintf('ln -s %s %s',infile,outfile);
      unix(cmd);
    end
  end
  
end