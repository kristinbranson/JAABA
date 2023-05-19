function jabFileToClassifier(injabfile,outmatfile)

jd = loadAnonymous(injabfile);
multiclassifier_stuff(numel(jd.classifierStuff)) = struct; 
mult_beh = struct('Lift','',  'Handopen','', 'Grab', '','Supinate','', 'Chew' ,'','Atmouth','');
beh = fieldnames(mult_beh);
for j=1:numel(jd.classifierStuff)
  inClassifierParams = jd.classifierStuff(j).params;
  fns = fieldnames(inClassifierParams);
  for i = 1:numel(fns)
    fn = fns{i};
    multiclassifier_stuff(j).(fn) = reshape([inClassifierParams.(fn)],size(inClassifierParams));
  end
  mult_beh.(beh{j}) = multiclassifier_stuff(j);  
end
save(outmatfile,'-struct', 'mult_beh','-v7.3');