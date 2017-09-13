jabfile = '/groups/branson/bransonlab/projects/JAABA/data/headfixedmice_adam/annotateddata/Jabfiles/Liftm76_pastonly.jab';
[success,isproblem,msgs,expdirs] = checkJabFile(jabfile,'frontside',true);

for i = find(isproblem),
  fprintf('\n%s:\n',expdirs{i});
  fprintf('  %s\n',msgs{i}{:});
end

for i = find(isproblem),
  expdir = expdirs{i};
  genAllFeatures(expdir,'doforce',true,'dorecurse',false,'frontside',true);
end