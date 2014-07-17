%% parameters

expdirs = {
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_26E12_AE_01_TrpA_Rig1Plate15BowlA_20111202T160314'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_64D11_AE_01_TrpA_Rig1Plate10BowlB_20110401T144120'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_21G01_AE_01_TrpA_Rig1Plate15BowlD_20110812T111328'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_55C09_AE_01_TrpA_Rig2Plate14BowlC_20110330T151605'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_69A12_AE_01_TrpA_Rig1Plate10BowlD_20110408T110922'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_68B12_AE_01_TrpA_Rig1Plate15BowlD_20120208T161323'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_92B03_AE_01_TrpA_Rig1Plate15BowlA_20111104T095956'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_20C03_AE_01_TrpA_Rig1Plate10BowlA_20110720T133402'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_16B01_AE_01_TrpA_Rig1Plate15BowlC_20120307T153820'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/GMR_50B07_AE_01_TrpA_Rig2Plate14BowlA_20110331T092823'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_67B01_AE_01_TrpA_Rig2Plate17BowlB_20111014T091451'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_22E02_AE_01_TrpA_Rig2Plate17BowlD_20111201T093552'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_13A08_AE_01_TrpA_Rig1Plate15BowlC_20111111T104418'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_73E12_AE_01_TrpA_Rig1Plate15BowlD_20120517T140351'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_52H10_AE_01_TrpA_Rig2Plate17BowlA_20110914T145037'
  '/groups/branson/home/robiea/Projects_data/JAABA/GroundTruth/specific_behaviors/GMR_17F12_AE_01_TrpA_Rig1Plate15BowlD_20120411T092053'
  };

scorefns = {'scores_Chasev7','scores_Righting','scores_Jump'};
behaviornames = {'Chase','Righting','Jump'};
rootdatadir = '/groups/branson/bransonlab/mayank/myFlyBowl';
nexps = numel(expdirs);
experiment_names = cell(1,nexps);
nexpsbad = 10;
for i = 1:nexps,
  [~,experiment_names{i}] = fileparts(expdirs{i});
  if i <= nexpsbad,
    expdirs{i} = fullfile(rootdatadir,experiment_names{i});
  end
end

nbehaviors = numel(scorefns);

T0 = 6000;
nframesplot = 5000;
buffert = 200;
T1 = T0+nframesplot-1;

behaviorcolors = lines(3);
yrad = .45;


%% 

% reset random number generator

% stream = RandStream.getGlobalStream;
% reset(stream);


% load the data

labels = nan(2,nexps,nframesplot,nbehaviors);
fliesselected = nan(2,nexps,nbehaviors);
for i = 1:nexps,
  for j = 1:nbehaviors,

    % load the trx
    trxfile = fullfile(expdirs{i},'registered_trx.mat');
    load(trxfile,'trx');
    
    % find one male and one female that is alive the entire interval
    nflies = numel(trx);
    isallowedmale = false(1,nflies);
    isallowedfemale = false(1,nflies);
    
    t0 = T0 + min([trx.firstframe]) - 1;
    t1 = t0+nframesplot-1;
    
    for fly = 1:nflies,
      if trx(fly).firstframe > t0-buffert || trx(fly).endframe < t1 + buffert,
        continue;
      end
      i0 = t0+trx(fly).off;
      i1 = t1+trx(fly).off;
      if all(strcmpi(trx(fly).sex(i0:i1),'M')),
        isallowedmale(fly) = true;
      elseif all(strcmpi(trx(fly).sex(i0:i1),'F')),
        isallowedfemale(fly) = true;
      end
    end
    flymale = randsample(find(isallowedmale),1);
    flyfemale = randsample(find(isallowedfemale),1);
    
    % load the scores
    scorefile = fullfile(expdirs{i},scorefns{j});
    scoredata = load(scorefile);
    labels(1,i,:,j) = scoredata.allScores.scores{flymale}(t0:t1) > 0;
    labels(2,i,:,j) = scoredata.allScores.scores{flyfemale}(t0:t1) > 0; 
    fliesselected(1,i,j) = flymale;
    fliesselected(2,i,j) = flyfemale;
    
  end
end

%% plot

hfig = 1;
figure(hfig);
clf;

labels1 = reshape(labels,[nexps*2,nframesplot,nbehaviors]);

h = cell(2*nexps,nbehaviors);
hlegend = nan(1,nbehaviors);
isfirst = true;
for j = 1:nbehaviors,
  for i = 1:nexps*2,
    [starts,ends] = get_interval_ends(labels1(i,:,j));
    ends = ends-1;
    nbouts = numel(starts);
    if nbouts == 0,
      continue;
    end
    x = [starts;starts;ends;ends;starts]+T0-1;
    y = i+yrad*[-1;1;1;-1;-1];
    h{i,j} = nan(1,nbouts);
    for k = 1:nbouts,
      h{i,j}(k) = patch(x(:,k),y,behaviorcolors(j,:),...
        'EdgeColor','none','MarkerFaceColor',behaviorcolors(j,:));
      if isfirst,
        isfirst = false;
        hold on;
      end
      if isnan(hlegend(j)),
        hlegend(j) = h{i,j}(1);
      end
    end
  end
end

line_names = cell(1,nexps);
for i = 1:nexps,
  [~,expname] = fileparts(expdirs{i});
  m = regexp(expname,'^GMR_(.*)_AE_01_TrpA_Rig.*$','tokens','once');
  line_names{i} = m{1};
end

legend(hlegend,behaviornames);
set(gca,'YTick',.5+1:2:2*nexps,'YTickLabel',line_names,'XLim',[T0-10,T1+10],'YLim',[0,2*nexps+1],'YDir','reverse','Box','off','TickDir','out');

%SaveFigLotsOfWays(hfig,'../figures/GMROut/Ethogram2')