function success=testJAABAGui()

% This tests to make sure that JAABA gui can label and train a classifier.
% Would be good to include setting features as part of it.
%%
success = false;
jabFileName='/groups/branson/home/kabram/bransonlab/projects/JAABA/test_data/testGui.jab';
expdir = '/groups/branson/home/kabram/bransonlab/projects/JAABA/test_data/GMR_71G01_AE_01_TrpA_Rig1Plate15BowlA_20120316T144027';
%% Open JAABA

SetUpJAABAPath;

try
  c=parcluster;
  if (c.NumWorkers>2) && (matlabpool('size')<1)
    matlabpool('open',c.NumWorkers-1);  % BJA: must save one for frame cache thread
  end
end
% Start JAABA.
JLabel('jabfile',jabFileName);

jObject = findall(0,'type','figure','name','JAABA');
if numel(jObject)~=1, 
  error('Could not find the jaaba window');
end
jhandles = guidata(jObject);

%%

jhandles.data.AddExpDir(expdir);
JLabel('modifyFilesDone',jObject,true);

%% Label Some data

JLabel('SetCurrentFrame',jhandles,1,4000,jObject,true);

nbouts = 100;
blen = 30;

frameStart = 4000;
for ndx = 1:nbouts
  jhandles = JLabel('SetLabelPlot',jhandles,frameStart+1,frameStart+blen,1,1);
  jhandles = JLabel('SetLabelPlot',jhandles,frameStart+blen+1,frameStart+2*blen,2,1);
  frameStart = frameStart + 2*blen;
end

JLabel('UpdatePlots',jhandles,...
              'refreshim',false, ...
              'refreshflies',true, ...
              'refreshtrx',false, ...
              'refreshlabels',true,...
              'refresh_timeline_manual',true,...
              'refresh_timeline_auto',false,...
              'refresh_timeline_suggest',false,...
              'refresh_timeline_error',true,...
              'refresh_timeline_xlim',false,...
              'refresh_timeline_hcurr',false,...
              'refresh_timeline_props',false,...
              'refresh_timeline_selection',false,...
              'refresh_curr_prop',false);
 
guidata(jObject,jhandles);            

%% Train

JLabel('pushbutton_train_Callback',jhandles.pushbutton_train,[],jhandles);

%% Save classifier.

jhandles.data.saveJabFile(tempname);


%% Quit

JLabel('figure_JLabel_CloseRequestFcn',jObject,[],jhandles);

matlabpool('close');

success = true;


