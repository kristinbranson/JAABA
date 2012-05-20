% data,labels are the training data and labels. 
% info has the information about these examples (expi,flies,t)
% newdata,newlabels are the data and labels for examples not part of training data. 
% newinfo has information about these examples.
% expdirs have the expdirs, trx has the loaded tracks for that experiment, and moviefilename is 'movie.ufmf';

%%
moviefilename = 'movie.ufmf';

clparams = struct('iter',100,'iter_updates',10,...
      'numSample',2500,'numBins',30,'CVfolds',7,...
      'baseClassifierTypes',{'Decision Stumps'},'baseClassifierSelected',1);

binVals = findThresholds(data,clparams);
bins = findThresholdBins(data,binVals,clparams);
[bagModels distmat] = doBagging(data,labels,[],binVals,bins,clparams);
newdistmat = findBagCoords(newdata,bagModels);

alldistmat = [distmat; newdistmat];

numex = size(alldistmat,1);
xdistmat = zeros(numex);
for ndx = 1:numex
  for indx = 1:numex
    xdistmat(indx,ndx) = nanmean(abs(alldistmat(indx,:)-alldistmat(ndx,:)));
  end
end

Z = mdscale(xdistmat,2,'Start','random');


%% Scatter plot the embedded points.

alllabels = [labels(:)' (newlabels(:)'+2)];

C = zeros(numel(alllabels,3));
for ndx = 1:numel(labels)
  switch labels(ndx)
    case 1
      C(ndx,:) = [1 0 0];
    case 2
      C(ndx,:) = [0 0 1];
      
    case 3
      C(ndx,:) = [1 0.6 0];
    case 4 
      C(ndx,:) = [0 0.6 1];
  end
end
% Scale by mfac so that the image patches are more separate
mfac = 6000;
figure; f = scatter(Z(:,1)*mfac,Z(:,2)*mfac,40,C,'filled');

%% Load images and tracks for training data.
im = {}; traj = {};
for ndx = 1:numel(info.expi)
  expi = info.expi(ndx);
  flies = info.flies(ndx);
  t = info.t(ndx);
  mname = fullfile(expdirs{expi},moviefilename);
  [im{ndx},traj{ndx}] = vid2im(mname,flies,t,trx{expi});
end

%% Load images and tracks for test data.
newim = {}; newtraj = {};
for ndx = 1:numel(newinfo.expi)
  expi = newinfo.expi(ndx);
  flies = newinfo.flies(ndx);
  t = newinfo.t(ndx);
  mname = fullfile(expdirs{expi},moviefilename);
  [newim{ndx},newtraj{ndx}] = vid2im(mname,flies,t,trx{expi});
end

%% Show the training images.

hold on;
curZ = Z(1:size(distmat,1),:);
for ndx = 1:numel(im)
  ii = uint8(repmat(im{ndx},[1 1 3]));
  newX = curZ(ndx,1)*mfac; newY = curZ(ndx,2)*mfac;
  image(newX,newY,ii);
  plot(traj{ndx}.Y+newX,traj{ndx}.X+newY)
end
axis equal;

%% Show the "out of training" images
hold on;
curZ = Z( (size(distmat,1)+1):end,:);
for ndx = 1:numel(newim)
  ii = uint8(repmat(newim{ndx},[1 1 3]));
  newX = curZ(ndx,1)*mfac; newY = curZ(ndx,2)*mfac;
  image(newX,newY,ii);
  plot(newtraj{ndx}.Y+newX,newtraj{ndx}.X+newY)
end
axis equal;

