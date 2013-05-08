function success=testConvertClassifierWithoutLabels()

projectFileName='/groups/branson/bransonlab/projects/JAABA/test_data/fooing_project.mat';
classifierFileName='/groups/branson/bransonlab/projects/JAABA/test_data/fooing_labels_and_classifier_without_labels.mat';
gtExpDirNames={'/groups/branson/bransonlab/projects/JAABA/test_data/GMR_82E08_AE_01_TrpA_Rig1Plate10BowlA_20110323T092430'}';               
jabFileName='/tmp/fooing_labels_and_classifier_without_labels.jab';

gum=GrandlyUnifyModel();
gum.setProjectFileName(projectFileName);
gum.setClassifierFileName(classifierFileName);
for i=1:length(gtExpDirNames)
  gum.addGTExpDirName(gtExpDirNames{i})
end
gum.convert(jabFileName);
gum=[];  %#ok

gtMode=false;
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()));
data.openJabFile(jabFileName,gtMode);

labels=data.labels;
if length(labels)~=2 ,
  data.closeJabFile();
  error('.jab file has the wrong number of labels');
end
if abs(labels(1).timestamp{1}(2)-735244.723198031)>1e-6
  data.closeJabFile();
  error('The one timestamp I checked is wrong');
end

data.closeJabFile();

success=true;

end
