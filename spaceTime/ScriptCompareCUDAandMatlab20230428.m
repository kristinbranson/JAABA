expdir = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe/test_hoghof/M274_20180814_v002';
originalexpdir = '/groups/branson/bransonlab/hantman_data/jab_experiments/M274Vglue2_Gtacr2_TH/20180814/M274_20180814_v002';
cudafeatfile = fullfile(expdir,'cuda_features_M274/features_cuda_averaged.mat');
cudafeat = load(cudafeatfile);
cudafeatfile0 = fullfile(expdir,'cuda_features_M274/hoghof.h5');
datasets = h5info(cudafeatfile0).Datasets;
cudafeat0 = struct;
for i = 1:numel(datasets),
  n = datasets(i).Name;
  cudafeat0.(n) = h5read(cudafeatfile0,['/',n]);
end
%cudafeatdirnew = '/groups/branson/bransonlab/patilr/realtime_classifier_test/cuda_features_M274/avgwithinplugin';
cudafeatdirnew = '/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe/test_hoghof/cudafeat';
cudafeatnew = struct;
cudafeatnew.front = csvread(fullfile(cudafeatdirnew,'hoghof_avg_front_biasjaaba.csv'));
cudafeatnew.side = csvread(fullfile(cudafeatdirnew,'hoghof_avg_side_biasjaaba.csv'));

featorder.side = readtable('/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe/test_hoghof/translatedindexes_mat2C_side.csv');
featorder.front = readtable('/groups/branson/home/bransonk/behavioranalysis/code/Jdetect/Jdetect/perframe/test_hoghof/translatedindexes_mat2C_front.csv');

matlabfeat = load(fullfile(expdir,'features.mat'));
matlabfeat0 = load(fullfile(originalexpdir,'features.mat'));


%%
prct = 99.9;
n = size(matlabfeat.curFeatures,2);
x = matlabfeat.curFeatures(:,1:n/2);
scale_hof = prctile(x(:),prct);
x = matlabfeat.curFeatures(:,n/2+1:end);
scale_hog = prctile(x(:),prct);
scale = zeros(1,n) + scale_hof;
scale(n/2+1:end) = scale_hog;

%%

nfront = size(cudafeatnew.front,2);
nside = size(cudafeatnew.side,2);
cudafeatnew.combined = [cudafeatnew.side(:,nside/2+1:end),...
  cudafeatnew.front(:,nside/2+1:end),...
  cudafeatnew.side(:,1:nside/2),...
  cudafeatnew.front(:,1:nside/2)];
cudafeatnew.reordered = nan(size(matlabfeat.curFeatures));
for i = 1:numel(featorder.side.Matlab),
  mati = featorder.side.Matlab(i);
  cudai = featorder.side.Cuda(i);
  cudafeatnew.reordered(:,mati) = cudafeatnew.side(:,cudai);
end
for i = 1:numel(featorder.front.Matlab),
  mati = featorder.front.Matlab(i);
  cudai = featorder.front.Cuda(i);
  cudafeatnew.reordered(:,mati) = cudafeatnew.front(:,cudai);
end

%%

winsize = 5;
cudafeat0.combined = combine_hog_hof_features(cudafeat0,true);
nframes = size(cudafeat0.combined,1);
fil = ones(winsize,1);
n = imfilter(ones(nframes,1),fil,0,'same');
cudafeat0.smoothed = imfilter(cudafeat0.combined,fil,0,'same')./n;

%%

%scale = 1;
figure(1);
clf;
hax = createsubplots(3,1,.05);

clim = [0,prctile([matlabfeat.curFeatures(:);cudafeat.curFeatures(:)],99)];

axi = 1;

axes(hax(axi));
imagesc((matlabfeat.curFeatures./scale)',clim);
title('matlab (computed 20230502)');
colorbar;

axi = axi + 1;
axes(hax(axi));
imagesc((cudafeatnew.reordered./scale)',clim);
title('cuda (computed 20230502, reordered)');
colorbar;

idx = any(~isnan(cudafeatnew.reordered),1);
axi = axi+1;
axes(hax(axi));
tmp = abs(cudafeatnew.reordered-matlabfeat.curFeatures)./scale;
diff_clim = prctile(tmp(:),99);
imagesc(((cudafeatnew.reordered(:,idx)-matlabfeat.curFeatures(:,idx))./scale(idx))',diff_clim*[-1,1]);
title('cuda - matlab');
colorbar;

colormap jet;
xlabel('Frame');
ylabel('Feature number');



%% 

nfeat = size(cudafeat.curFeatures,2);
nfeatpertype = nfeat/2;
dfeat = sum(abs(cudafeat.curFeatures-matlabfeat.curFeatures),1);
flowfeatnums = zeros(1,2);
flowfeatnums(1) = argmax(dfeat(1:nfeatpertype));
dfeatmedian = median(dfeat(1:nfeatpertype));
flowfeatnums(2) = argmin(abs(dfeat(1:nfeatpertype)-dfeatmedian));

hogfeatnums = zeros(1,2);
hogfeatnums(1) = argmax(dfeat(nfeatpertype+1:end));
dfeatmedian = median(dfeat(nfeatpertype+1:end));
hogfeatnums(2) = argmin(abs(dfeat(nfeatpertype+1:end)-dfeatmedian));

figure(2);
clf;

for i = 1:2,
  subplot(4,1,i);
  hmatlab = plot(matlabfeat.curFeatures(:,flowfeatnums(i)));
  hold on;
  hcuda = plot(cudafeat.curFeatures(:,flowfeatnums(i)));
  hmatlab0 = plot(matlabfeat0.curFeatures(:,flowfeatnums(i)));
  legend([hmatlab,hcuda,hmatlab0],{'new matlab','cuda','old matlab'});
  title(sprintf('Flow feature %d',flowfeatnums(i)));
end

for i = 1:2,
  subplot(4,1,i+2);
  hmatlab = plot(matlabfeat.curFeatures(:,hogfeatnums(i)));
  hold on;
  hcuda = plot(cudafeat.curFeatures(:,hogfeatnums(i)));
  hmatlab0 = plot(matlabfeat0.curFeatures(:,hogfeatnums(i)));
  legend([hmatlab,hcuda,hmatlab0],{'new matlab','cuda','old matlab'});
  title(sprintf('HOG feature %d',hogfeatnums(i)));
end

