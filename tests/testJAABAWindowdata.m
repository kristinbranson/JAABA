function success=testJAABAWindowdata()

% This tests to make sure that loading cached windowdata and computing
% windowdata again learn the same classifier.

% AL20150602. This test doesn't work when I run it on the current codebase;
% I think it's because window feature computation is no longer
% deterministic, even with the setting of RNG state below.
% ComputeWindowFeatures() is called within a parfor and has a call to
% randsample.

%%
success = false;
jabFileName='/groups/branson/home/kabram/bransonlab/projects/JAABA/test_data/test_windowdata.jab';
gtMode = false;
data=JLabelData('setstatusfn',@(str)(fprintf('%s\n',str)), ...
                'clearstatusfn',@()(nop()));
data.openJabFile(jabFileName,gtMode);

if matlabpool('size')>0, matlabpool('close'); end
oldcl = data.classifier;

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

data.Train;
newcl1 = data.classifier;

data.setWindowFeaturesParams(data.windowfeaturesparams);

s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);

data.Train;
newcl2 = data.classifier;

if ~isequal(oldcl,newcl1) || ~isequal(newcl1,newcl2),
  error('Classifiers trained are not same');
  
end
success = true;
