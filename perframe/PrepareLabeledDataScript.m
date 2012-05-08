% script for preparing labeled behavior data

% set up path
if ispc,
  JCtrax_path = 'C:\Code\JCtrax';
  FlyBowlAnalysis_path = 'C:\Code\FlyBowlAnalysis';
  rootdatadir = 'C:\Code\Jdetect\larva\fly_data\TrainingData';
  settingsdir = 'C:\Code\FlyBowlAnalysis\settings';
else
  JCtrax_path = '/groups/branson/home/bransonk/tracking/code/JCtrax';
  FlyBowlAnalysis_path = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis';
  settingsdir = '/groups/branson/bransonlab/projects/olympiad/FlyBowlAnalysis/settings';
end
addpath(fullfile(JCtrax_path,'misc'));
addpath(fullfile(JCtrax_path,'filehandling'));
addpath(FlyBowlAnalysis_path);

%% load test data

ctraxfilestr = 'fixed_ctrax_results.mat';
trxfilestr = 'registered_trx.mat';
registrationdatafilestr = 'registrationdata.mat';
analysis_protocol = '20110407';
perframefns = {'dtheta','velmag','dist2wall','dnose2ell','dcenter','dell2nose'};

labelmatnames = {
  'labeledsharpturns_movie01_pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110202T105734_fly1_11.mat'
  'labeledsharpturns_movie01_pBDPGAL4U_TrpA_Rig1Plate10BowlA_20110202T105734_fly2_08.mat'
  'labeledsharpturns_movie02_pBDPGAL4U_TrpA_Rig2Plate14BowlA_20110202T110111_fly1_01.mat'
  'labeledsharpturns_movie02_pBDPGAL4U_TrpA_Rig2Plate14BowlA_20110202T110111_fly2_15.mat'
  'labeledsharpturns_movie03_pBDPGAL4U_TrpA_Rig2Plate14BowlB_20110202T110116_fly1_12.mat'
  'labeledsharpturns_movie03_pBDPGAL4U_TrpA_Rig2Plate14BowlB_20110202T110116_fly2_03.mat'
  'labeledsharpturns_movie07_GMR_14G03_AE_01_TrpA_Rig2Plate14BowlA_20110202T102516_fly2_17.mat'
  'labeledsharpturns_movie08_GMR_14H07_AE_01_TrpA_Rig1Plate10BowlB_20110202T084231_fly1_09.mat'
  'labeledsharpturns_movie09_GMR_16E02_AE_01_TrpA_Rig1Plate10BowlC_20110202T140143_fly1_09.mat'
  };

for expdiri = 1:numel(labelmatnames),

  labeldata = load(fullfile(rootdatadir,labelmatnames{expdiri}));
  
  % replace root dir
  [pathstr,moviefilestr] = myfileparts(labeldata.moviename);
  %moviefilestr = [moviefilestr,moviefileext];
  [pathstr,expdir] = myfileparts(pathstr);
  expdir = fullfile(rootdatadir,expdir);
  moviename = fullfile(expdir,moviefilestr);
  
  % load trx
  ctraxfile = fullfile(expdir,ctraxfilestr);
  load(ctraxfile,'trx');
  
  % apply latest spatial registration
  rd = load(fullfile(expdir,registrationdatafilestr));
  rd = detectRegistrationMarks('registrationData',rd);
  
  trx = ApplyRegistration(trx,rd);
  
  % save to fixed_registered_trx mat file
  trxfile = fullfile(expdir,trxfilestr);
  save(trxfile,'trx');
  
  % compute some per-frame features
  if expdiri > 1,
    FlyBowlComputePerFrameFeatures(expdir,'analysis_protocol',analysis_protocol,...
      'settingsdir',settingsdir,'perframefns',perframefns);
  end
  
  clear trx;
  
  fly = labeldata.fly;
  
  trx = Trx('analysis_protocol',analysis_protocol,...
    'settingsdir',settingsdir);
  trx.AddExpDir(expdir);

  % start of labeled sequence, relative to video
  video_labelstart = labeldata.t0tolabelcurr;
  video_labelend = labeldata.t1tolabelcurr;
  % start of labeled sequence, relative to trajectory indices
  trk_labelstart = video_labelstart + trx(fly).off;
  trk_labelend = video_labelend + trx(fly).off;
  
  % load all the per-frame data for fly
  perframedata = cell(1,numel(perframefns));
  for fni = 1:numel(perframefns);
    perframefn = perframefns{fni};
    perframedata{fni} = trx(fly).(perframefn);
  end

  % create 0/1/nan per-frame labels
  n = video_labelend-video_labelstart+1;
  label = zeros(1,n-1);
  for labeli = 1:numel(labeldata.segstarts),
    if strcmpi(labeldata.labels{labeli},'unknown'),
      val = nan;
    else
      val = 1;
    end
    label(labeldata.segstarts(labeli)-video_labelstart+1:labeldata.segends(labeli)-video_labelstart+1) = val;    
  end
  behaviorname = setdiff(labeldata.labels,{'unknown'});
  
  hold off;
  plot(1:n-1,abs(perframedata{1}(trk_labelstart:trk_labelend-1)),'k.-');
  hold on;
  plot(find(label==1),abs(perframedata{1}(trk_labelstart-1+find(label==1))),'r.');
  plot(find(isnan(label)),abs(perframedata{1}(trk_labelstart-1+find(isnan(label)))),'b.');
  axisalmosttight;
  drawnow;

  [~,savename] = myfileparts(labelmatnames{expdiri});
  savename = fullfile(rootdatadir,['perframe_',savename]);
  save(savename,'perframefns','perframedata','label',...
    'video_labelstart','video_labelend','trk_labelstart','trk_labelend',...
    'behaviorname','expdir','moviename','ctraxfile','trxfile');

end