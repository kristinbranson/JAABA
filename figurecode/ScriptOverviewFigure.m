%% script for making plots used to demonstrate trx, per-frame, window features

%% set up path

addpath ../perframe;
addpath ../perframe/compute_perframe_features;
addpath ../misc;
addpath ../filehandling;
outfigdir = 'C:\Code\Jdetect\figures\OverviewOut';

%% parameters

rootdatadir = 'C:\Data\JAABA\FlyBowl';
experiment_name = 'pBDPGAL4U_TrpA_Rig1Plate15BowlB_20110922T145928';
fly = 4;
t = 6823;

%% load the data

trx = Trx('trxfilestr','registered_trx.mat','moviefilestr','movie.ufmf','perframedir','perframe');
expdir = fullfile(rootdatadir,experiment_name);
trx.AddExpDir(expdir);
moviename = trx.movienames{1};
[readframe,nframes,fid,headerinfo] = get_readframe_fcn(moviename);

%% show the trajectory

tradius = 50;

hfig = 1;
figure(hfig);
set(hfig,'Units','pixels','Position',[10 10 1026 964]);
clf;

plot(trx(fly).x,trx(fly).y,'-','color',[.4,.4,.4]);
hold on;
idx = t-tradius+trx(fly).off:t+tradius+trx(fly).off;
plot(trx(fly).x(idx),trx(fly).y(idx),'-','linewidth',3,'color','k');
axis equal ij;
axisalmosttight;
set(gca,'Box','off');

SaveFigLotsOfWays(hfig,fullfile(outfigdir,sprintf('ExampleTrajectory_%s_fly%02d',experiment_name,fly)));

%% choose per-frame features

T0 = min(cellfun(@(x)x(1),trx.timestamps));

files = dir(fullfile(trx.expdirs{1},'perframe','*.mat'));
perframefns = cellfun(@(x) x(1:end-4),{files.name},'UniformOutput',false);
ignorefns = {'a','b','x','y','theta','sex','timestamps','dt','x_mm','y_mm','xnose_mm','ynose_mm','arena_angle','arena_r','signdtheta'};
perframefns = setdiff(perframefns,ignorefns);
idx = ~cellfun(@isempty,regexp(perframefns,'^closestfly_.*$'));
perframefns(idx) = [];

%% plot per-frame features

naxr = 5;
naxc = 3;
nax = naxr*naxc;

t0 = t-tradius;
t1 = t+tradius;
i0 = t0+trx(fly).off;
i1 = t1+trx(fly).off;

hfig = 2;
fignum = 1;
for i = 1:nax:numel(perframefns),
  figure(hfig);
  set(hfig,'Units','pixels','Position',[20 20 2243 1030]);
  clf;
  hax = createsubplots(naxr,naxc,[.025,.025;.025,.01]);
  for j = i:min(numel(perframefns),i+nax-1),
    %axes(hax(j-i+1));
    fn = perframefns{j};
    plot(hax(j-i+1),trx(fly).timestamps(i0:i1)-T0,trx(fly).(fn)(i0:i1),'k-');
    axisalmosttight([],hax(j-i+1));
    set(hax(j-i+1),'XLim',[trx(fly).timestamps(i0-1)-T0,trx(fly).timestamps(i1-1)-T0],'Box','off');
    if mod(j-i,naxr) ~= naxr-1,
      set(hax(j-i+1),'XTickLabel',{});
    end
    ylabel(hax(j-i+1),perframefns{j},'Interpreter','none');
  end
  SaveFigLotsOfWays(hfig,fullfile(outfigdir,sprintf('ExamplePerFrameFeatures_%s_fly%02d_%d',experiment_name,fly,fignum)));
  hfig = hfig + 1;  
  fignum = fignum+1;
end

%% compute window features

featureparamsfilename = '../perframe/params/WindowFeatures_Touch.xml';
featureconfigfilename = '../perframe/params/featureConfig.xml';

[windowfeaturesparams,windowfeaturescellparams,basicFeatureTable,featureWindowSize] = ...
  ReadPerFrameParams(featureparamsfilename,featureconfigfilename); 

perframefns1 = fieldnames(windowfeaturescellparams);

x_curr_all = cell(1,numel(perframefns1));
feature_names_all = cell(1,numel(perframefns1));


for j = 1:numel(perframefns1),
  fn = perframefns1{j};
        
  % get per-frame data
  perframedata = trx(fly).(fn);
  
  i11 = min(i1,numel(perframedata));
  [x_curr,feature_names_curr] = ...
    ComputeWindowFeatures(perframedata,windowfeaturescellparams.(fn){:},'t0',i0,'t1',i11);
  if any(imag(x_curr(:)))
    fprintf('Feature values are complex, check input\n');
  end
  
  if i11 < i1,
    x_curr(:,end+1:end+i1-i11) = nan;
  end
        
  x_curr_all{j} = x_curr;
  feature_names_all{j} = feature_names_curr;
        
end

%% plot window features

rng(0);

hfig = 100;

fignum = 0;
nwindowfeatures = naxc;
for i = 1:numel(perframefns1),
  
  pfn = perframefns1{i};
  
  if mod(i,naxr) == 1,
    hfig = hfig+1;
    fignum = fignum+1;
    figure(hfig);
    set(hfig,'Units','pixels','Position',[98 143 1804 963]);
    clf;
    hax = reshape(createsubplots(naxr,naxc,[.03,.03;.04,.04]),[naxr,naxc]);
    r = 1;
  else
    r = r+1;
  end
  
  nw = numel(feature_names_all{i});
  js = sort(randsample(nw,nwindowfeatures))';
  
  for jj = 1:nwindowfeatures,
    j = js(jj);
    tmp = feature_names_all{i}{j}(5:end);
    tmp(2:2:end) = cellfun(@num2str,tmp(2:2:end),'UniformOutput',false);
    wfn = [sprintf('%s_%s',feature_names_all{i}{j}{2},feature_names_all{i}{j}{4}),...
      sprintf('_%s%s',tmp{:})];
    plot(hax(r,jj),trx(fly).timestamps(i0:i1)-T0,x_curr_all{i}(j,:),'k-');
    axisalmosttight([],hax(r,jj));
    set(hax(r,jj),'XLim',[trx(fly).timestamps(i0-1)-T0,trx(fly).timestamps(i1-1)-T0],'Box','off');
    if r ~= naxr,
      set(hax(r,jj),'XTickLabel',{});
    end
    %if jj == 1,
    ylabel(hax(r,jj),perframefns{i},'Interpreter','none');
    xlabel(hax(r,jj),wfn,'Interpreter','none');
    %else
    %  ylabel(hax(r,jj),wfn,'Interpreter','none');
    %end
  end
  if mod(i,naxr) == 0,
    SaveFigLotsOfWays(hfig,fullfile(outfigdir,sprintf('ExampleWindowFeatures_%s_fly%02d_%02d',experiment_name,fly,fignum)));
  end
end

%%

hfig = 200;
plotall_perframefns = {'velmag_ctr','dnose2ell'};
for fnii = 1:numel(plotall_perframefns),
  fignum = 1;
  pfn = plotall_perframefns{fnii};
  fni = find(strcmp(perframefns1,pfn),1);
  nw = numel(feature_names_all{fni});
  for i = 1:nax:nw,

    figure(hfig);
    set(hfig,'Units','pixels','Position',[98 143 1804 963]);
    clf;
    hax = createsubplots(naxr,naxc,[.03,.03;.04,.04]);
    for j = i:min(nw,i+nax-1),
      %axes(hax(j-i+1));
      tmp = feature_names_all{fni}{j}(5:end);
      tmp(2:2:end) = cellfun(@num2str,tmp(2:2:end),'UniformOutput',false);
      wfn = [sprintf('%s_%s',feature_names_all{fni}{j}{2},feature_names_all{fni}{j}{4}),...
        sprintf('_%s%s',tmp{:})];
      plot(hax(j-i+1),trx(fly).timestamps(i0:i1)-T0,x_curr_all{fni}(j,:),'k-');
      axisalmosttight([],hax(j-i+1));
      set(hax(j-i+1),'XLim',[trx(fly).timestamps(i0-1)-T0,trx(fly).timestamps(i1-1)-T0],'Box','off');
      if mod(j-i,naxr) ~= naxr-1,
        set(hax(j-i+1),'XTickLabel',{});
      end
      %ylabel(hax(j-i+1),pfn,'Interpreter','none');
      xlabel(hax(j-i+1),wfn,'Interpreter','none');
    end
    pause(1);
    SaveFigLotsOfWays(hfig,fullfile(outfigdir,sprintf('ExampleWindowFeatures_%s_%s_fly%02d_%d',pfn,experiment_name,fly,fignum)));
    hfig = hfig + 1;
    fignum = fignum+1;
  end
end