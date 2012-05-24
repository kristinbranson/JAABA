% data,labels are the training data and labels. 
% info has the information about these examples (expi,flies,t)
% newdata,newlabels are the data and labels for examples not part of training data. 
% newinfo has information about these examples.
% expdirs have the expdirs, trx has the loaded tracks for that experiment, and moviefilename is 'movie.ufmf';

%% Bouts.
bouts = struct('ndx',[],'label',[],'timestamp',[]);
for expNdx = 1:obj.nexps
  for flyNdx = 1:obj.nflies_per_exp(expNdx)
    curLabels = obj.GetLabels(expNdx,flyNdx);
    for boutNum = 1:numel(curLabels.t0s)
      idx =  obj.FlyNdx(expNdx,flyNdx) & ...
        obj.windowdata.t >= curLabels.t0s(boutNum) & ...
        obj.windowdata.t < curLabels.t1s(boutNum);
      bouts.ndx(end+1,:) = obj.FlyNdx(expNdx,flyNdx) & ...
        obj.windowdata.t >= curLabels.t0s(boutNum) & ...
        obj.windowdata.t < curLabels.t1s(boutNum);
      bouts.label(end+1) = find(strcmp(obj.labelnames,curLabels.names{boutNum}));
      bouts.timestamp(end+1) = curLabels.timestamp(boutNum);
    end
    
  end
end

zz = bouts.timestamp;

moviefilename = 'movie.ufmf';

clparams = struct('iter',100,'iter_updates',10,...
      'numSample',2500,'numBins',30,'CVfolds',7,...
      'baseClassifierTypes',{'Decision Stumps'},'baseClassifierSelected',1);

allbagModels = {}; celldistmat = {};

%% To show the images.. load the trx.
for ndx = 1:obj.nexps
  trxfile = fullfile(obj.expdirs{ndx},obj.GetFileName('trx'));
  trx{ndx} = load_tracks(trxfile);
end
expdirs = obj.expdirs;
%%

for outqndx = 1:numel(Q.traintimes_chaseactive)-1
%%
  qq(1) = Q.traintimes_chaseactive(outqndx);
  if outqndx+1>numel(Q.traintimes_chaseactive)
    qq(2) = Q.traintimes_chaseactive(end);
  else
    qq(2) = Q.traintimes_chaseactive(outqndx+1);
  end

%%

oldbouts = find(zz<qq(1));
oldidx = false(size(bouts.ndx(1,:)));
oldidxall = false(size(bouts.ndx(1,:)));
for ndx = oldbouts(:)'
  nn = find(bouts.ndx(ndx,:),1);
  xx = find(bouts.ndx(ndx,:),1,'last');
  ii = round( (nn+xx)/2);
  oldidx(ii) = true;
  oldidxall = oldidxall | bouts.ndx(ndx,:);
end

oldndxintoall = find(oldidx(oldidxall));
data = obj.windowdata.X(oldidxall,:);
labels = obj.windowdata.labelidx_cur(oldidxall);
info = struct;
info.expi = obj.windowdata.exp(oldidxall);
info.flies = obj.windowdata.flies(oldidxall);
info.t = obj.windowdata.t(oldidxall);


newbouts = find(zz>=qq(1) & zz<qq(2) );
newidx = false(size(bouts.ndx(1,:)));
newidxall = false(size(bouts.ndx(1,:)));
for ndx = newbouts(:)'
  nn = find(bouts.ndx(ndx,:),1);
  xx = find(bouts.ndx(ndx,:),1,'last');
  ii = round( (nn+xx)/2);
  newidx(ii) = true;
  newidxall = newidxall | bouts.ndx(ndx,:);
end

newdata = obj.windowdata.X(newidxall,:);
newlabels = obj.windowdata.labelidx_cur(newidxall);
newinfo = struct;
newinfo.expi = obj.windowdata.exp(newidxall);
newinfo.flies = obj.windowdata.flies(newidxall);
newinfo.t = obj.windowdata.t(newidxall);
fprintf('Training,pos:%d,neg:%d Test:pos%d,neg:%d\n',nnz(labels==1),...
  nnz(labels==2),nnz(newlabels==1),nnz(newlabels==2));


%%
%{
binVals = findThresholds(data,clparams);
bins = findThresholdBins(data,binVals);
[bagModels distmat] = doBagging(data,labels,[],binVals,bins,clparams);
allbagModels{outqndx} = bagModels;
celldistmat{outqndx} = distmat;
%}

%%
bagModels = allbagModels{outqndx};
distmat = celldistmat{outqndx};

%%
newdistmat = findBagCoords(newdata,bagModels);

alldistmat = [distmat; newdistmat];
alllabels = [labels(:)' (newlabels(:)'+2)];

numex = size(alldistmat,1);
xdistmat = zeros(numex);
for ndx = 1:numex
    tt = repmat(alldistmat(ndx,:),numex,1);
    xdistmat(:,ndx) = nanmean(abs(tt-alldistmat),2);
end

wt = getWeights(sign([labels(:)' newlabels(:)']-1.5));
wt = sqrt(wt*wt');
% wt = ones(size(xdistmat));
Z = mdscale(xdistmat,2,'Start','random','Weights',wt);
allnlabels = [labels(:)' newlabels(:)'];
R(1) = sum(Z(1:numel(labels),1).*sign(labels-1.5));
R(2) = sum(Z(1:numel(labels),2).*sign(labels-1.5));
Z = Z(:,1:2);
theta = -atan2(R(2),R(1));
rotmat = [cos(theta) sin(theta); -sin(theta) cos(theta)];
Zr = Z*rotmat;
Zold = Z;
Z = Zr;
Z(:,1) = (Z(:,1)-min(Z(:,1)))/(max(Z(:,1))-min(Z(:,1)));
Z(:,2) = (Z(:,2)-min(Z(:,2)))/(max(Z(:,2))-min(Z(:,2)));


%%
%{
end
save embedBagsChase allbagModels celldistmat;
%}
%% Scatter plot the embedded points.

wt = getWeights(sign(labels-1.5));
[~,bb] = loglossboostLearnRandomFeatures(data,sign(labels-1.5),100,wt,binVals,bins,clparams);
allpred = myBoostClassify([data;newdata],bb);
tr = prctile(abs(allpred),80);
allnlabels = [labels(:)' newlabels(:)'];

%%
C = zeros(numel(alllabels),3);
% S = max(10,1./(1+exp(allpred.*(sign(allnlabels'-1.5))))*100);
S = [];
for ndx = 1:numel(alllabels)
  switch alllabels(ndx)
    case 1
      C(ndx,:) = [1 0 0];
    case 2
      C(ndx,:) = [0 0 1];
      
    case 3
      C(ndx,:) = [1 0.6 0];
    case 4 
      C(ndx,:) = [0 0.6 1];
  end
  
%  S(ndx) = 1/(min(1,abs(allpred(ndx)/tr)))*10;
  
%   if sign(allpred(ndx))== sign(allnlabels(ndx)-1.5)
%     S(ndx) = 20;
%   else
%     S(ndx) = 100;
%   end
end
% Scale by mfac so that the image patches are more separate

%%
%{
mfac = 6000;
figure; f = scatter(Z(:,1)*mfac,Z(:,2)*mfac,S,C,'filled');

%% Connect frames from same bouts

hold on;
curZ = Z(1:numel(labels),:);
for ndx = 1:numel(bouts.timestamp)
  
  if bouts.timestamp(ndx)>qq(1); continue; end;
  curidx = bouts.ndx(ndx,:);
  lndx = find(curidx(oldidxall));
  if bouts.label(ndx)==1
    col = [1 0 0];
  else
    col = [0 0 1];
  end
  plot(curZ(lndx,1)*mfac,curZ(lndx,2)*mfac,'Color',col);
end

curZ = Z(numel(labels)+1:end,:);
for ndx = 1:numel(bouts.timestamp)
  if bouts.timestamp(ndx)>qq(2) || bouts.timestamp(ndx)<qq(1); continue; end;
  curidx = bouts.ndx(ndx,:);
  lndx = find(curidx(newidxall));
  if bouts.label(ndx)==1
    col = [1 0.6 0];
  else
    col = [0 0.6 1];
  end
  plot(curZ(lndx,1)*mfac,curZ(lndx,2)*mfac,'Color',col);
end

%% To show the images.. load the trx.

for ndx = 1:obj.nexps
  trxfile = fullfile(obj.expdirs{ndx},obj.GetFileName('trx'));
  trx{ndx} = load_tracks(trxfile);
end

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

%% Show the most negative images on out of training bouts

hold on;

curZ = Z(numel(labels)+1:end,:);
for ndx = 1:numel(bouts.timestamp)
  if bouts.timestamp(ndx)>qq(2) || bouts.timestamp(ndx)<qq(1); continue; end;
  curidx = bouts.ndx(ndx,:);
  lndx = find(curidx(newidxall));
  curScores = myBoostClassify(newdata(lndx,:),bb);
  [~,minndx] = min(curScores*sign(bouts.label(ndx)-1.5));
  curlndx = lndx(minndx);
  ii = uint8(repmat(newim{curlndx},[1 1 3]));
  newX = curZ(curlndx,1)*mfac; newY = curZ(curlndx,2)*mfac;
  image(newX,newY,ii);
  plot(newtraj{curlndx}.Y+newX,newtraj{curlndx}.X+newY)
end
axis equal;

%% Show the most negative images on  training bouts

hold on;

curZ = Z(1:numel(labels),:);
for ndx = 1:numel(bouts.timestamp)
  if bouts.timestamp(ndx)>qq(1), continue; end;
  curidx = bouts.ndx(ndx,:);
  lndx = find(curidx(oldidxall));
  curScores = myBoostClassify(data(lndx,:),bb);
  [~,minndx] = min(curScores*sign(bouts.label(ndx)-1.5));
  curlndx = lndx(minndx);
  ii = uint8(repmat(im{curlndx},[1 1 3]));
  newX = curZ(curlndx,1)*mfac; newY = curZ(curlndx,2)*mfac;
  image(newX,newY,ii);
  plot(traj{curlndx}.Y+newX,traj{curlndx}.X+newY)
end
axis equal;

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
%
%%
mfac = 1200; sfac = 30;
clear Zm;
Zm(:,1) = Z(:,1)-min(Z(:,1));
Zm(:,2) = Z(:,2)-min(Z(:,2));

zszx = max(Zm(:,1))*mfac/sfac;
zszy = max(Zm(:,2))*mfac/sfac;
bimg = ones(round(zszy)+1,round(zszx)+1,3);

counts = zeros(round(zszy)+1,round(zszx)+1);
for ndx= 1:numel(labels)
  curx = round( Zm(ndx,1) *mfac/sfac)+1;
  cury = round( Zm(ndx,2) *mfac/sfac)+1;
  bimg(cury,curx,:) = reshape(C(ndx,:),[1 1 3]);
end

bbimg = imfilter(bimg,fspecial('gaussian',2*3*8,8),'symmetric');
hsvbimg = rgb2hsv(bbimg);
satimg = hsvbimg(:,:,2);
satimg = satimg-mean(satimg(:))+0.1;
satimg(satimg>1) = 1;
hsvbimg(:,:,2) = satimg;
bbimg = hsv2rgb(hsvbimg);
fbimg = imresize(bbimg,sfac);


fbimg(fbimg(:)>1) = 1;


figure; image(fbimg);
%}

%%

mfac = 500; sfac = 3; 

clear Zm;
Zm(:,1) = Z(:,1)-min(Z(:,1))+0.15;
Zm(:,2) = Z(:,2)-min(Z(:,2))+0.15;

zszx = (max(Zm(:,1))+0.15)*mfac/sfac;
zszy = (max(Zm(:,2))+0.15)*mfac/sfac;

numpos = nnz(labels==1);
numneg = nnz(labels==2);
ratio = numneg/numpos;
rad = sqrt(ratio);
totArea = 4*rad.^2*numpos + numneg;
bfac = sqrt(8000/totArea);

bimgy = round(zszy)+round(bfac+rad)+sfac;
bimgx = round(zszx)+round(bfac+rad)+sfac;
bimg = ones(bimgy,bimgx,3);

for ndx= 1:numel(labels)
  curx = round( Zm(ndx,1) *mfac/sfac)+1;
  cury = round( Zm(ndx,2) *mfac/sfac)+1;
 
  bimg(cury,curx,:) = reshape(C(ndx,:),[1 1 3]);
  if labels(ndx) ==1
    bboxy = cury+round(-rad*bfac:rad*bfac);
    bboxy(bboxy<1) = [];
    bboxx = curx+round(-rad*bfac:rad*bfac);
    bboxx(bboxx<1) = [];
    bimg(bboxy,bboxx,1) = C(ndx,1);
    bimg(bboxy,bboxx,2) = C(ndx,2);
    bimg(bboxy,bboxx,3) = C(ndx,3);
  else
    bboxy = cury+round(-bfac:bfac);
    bboxy(bboxy<1) = [];
    bboxx = curx+round(-bfac:bfac);
    bboxx(bboxx<1) = [];
    bimg(bboxy,bboxx,1) = C(ndx,1);
    bimg(bboxy,bboxx,2) = C(ndx,2);
    bimg(bboxy,bboxx,3) = C(ndx,3);

  end
  
end
aoccupied = bimg(:,:,2)~=1;
fprintf('Area occupied for %d:%d\n ',outqndx,nnz(aoccupied));
bbimg = imfilter(bimg,fspecial('gaussian',2*3*32,32),'symmetric');
fbimg = imresize(bbimg,sfac);


fbimg(fbimg(:)>1) = 1;
%%

hfig = figure; image(fbimg);

%
hold on;
% scatter(Zm(numel(labels)+1:end,1)*mfac,Zm(numel(labels)+1:end,2)*mfac,...
%   S(numel(labels)+1:end),C(numel(labels)+1:end,:),'filled');
% hold on;

occupied = false(size(fbimg));
rem = 55;

curZ = Zm(numel(labels)+1:end,:);
px = []; py = []; 
bx_p = []; by_p = [];
bx_n = []; by_n = [];
bx_np = []; by_np = [];
bx_nn = []; by_nn = [];
phandle = [];
count = 1;
for ndx = 1:numel(bouts.timestamp)
  if bouts.timestamp(ndx)>qq(2) || bouts.timestamp(ndx)<qq(1); continue; end;
  curidx = bouts.ndx(ndx,:);
  lndx = find(curidx(newidxall));
  curScores = myBoostClassify(newdata(lndx,:),bb);
  [~,minndx] = min(curScores*sign(bouts.label(ndx)-1.5));
  curlndx = lndx(minndx);
  expi = newinfo.expi(curlndx);
  flies = newinfo.flies(curlndx);
  t = newinfo.t(curlndx);
  mname = fullfile(expdirs{expi},moviefilename);
  [im,traj] = vid2im(mname,flies,t,trx{expi});
  
  ii = uint8(repmat(flipud(im),[1 1 3]));
  ii(:,[1:rem end-rem-1:end],:)=[];
  halfsz = round(size(ii,1)/2);
  bottomout = size(ii,1)-halfsz-30+1;
  ii(halfsz+30:end,:,:) = [];
  
  alpha = (ii(:,:,2)<170);
  newX = curZ(curlndx,1)*mfac-size(ii,1); 
  newY = curZ(curlndx,2)*mfac-size(ii,2);

  testX = floor(newX:newX+size(ii,1)); testX(testX<1) = []; testX(testX>size(occupied,1))=[];
  testY = floor(newY:newY+size(ii,1)); testY(testY<1) = []; testY(testY>size(occupied,2))=[];
  oo = occupied(testX,testY);
  if(nnz(oo(:))/numel(oo)>1), continue; end;

  if bouts.label(ndx) == 1
    bx_np = [bx_np nan newX newX newX+size(ii,2) newX+size(ii,2) newX];
    by_np = [by_np nan newY newY+size(ii,1) newY+size(ii,1) newY newY];
  else
    bx_nn = [bx_nn nan newX newX newX+size(ii,2) newX+size(ii,2) newX];
    by_nn = [by_nn nan newY newY+size(ii,1) newY+size(ii,1) newY newY];
  end
  

  occupied(testX,testY) = true;
  
%   if bouts.label(ndx) == 1
%     plot(newX,newY,'Marker','*','MarkerEdgeColor',[1 0.6 0],'MarkerSize',3);
%   else
%     plot(newX,newY,'Marker','*','MarkerEdgeColor',[0 0.6 1],'MarkerSize',3);
%   end
%   
%   hfignew = figure; 
%   image(1,1,ii,'AlphaData',alpha);
%   hold on; axis image; axis off; axis tight;
%   plot(-traj.Y+hszx-rem,-traj.X+hszx,'LineStyle','none','Marker','.','MarkerSize',1.5,'MarkerEdgeColor','k');
%   fname = sprintf('/groups/branson/home/kabram/paperFigs/classifier_chase_new_plotted_%d_individual_%d',outqndx,count);
%   count = count+1;
%   SaveFigLotsOfWays(hfignew,fname);
%   close(hfignew);

  image(newX,newY,ii,'AlphaData',alpha);
  hszx = size(ii,1) + bottomout;
  px = [px nan -traj.Y+newX+hszx-rem];
  py = [py nan -traj.X+newY+hszx];

end

curZ = Zm(1:numel(labels),:);
for ndx = 1:numel(bouts.timestamp)
  if bouts.timestamp(ndx)>qq(1), continue; end;
  curidx = bouts.ndx(ndx,:);
  lndx = find(curidx(oldidxall));
  curlndx = round(median(lndx));
  expi = info.expi(curlndx);
  flies = info.flies(curlndx);
  t = info.t(curlndx);
  mname = fullfile(expdirs{expi},moviefilename);
  [im,traj] = vid2im(mname,flies,t,trx{expi});

  ii = uint8(repmat(flipud(im),[1 1 3]));
  ii(:,[1:rem end-rem-1:end],:)=[];
  halfsz = round(size(ii,1)/2);
  bottomout = size(ii,1)-halfsz-30+1;
  ii(halfsz+30:end,:,:) = [];
  alpha = (ii(:,:,2)<170);
  newX = curZ(curlndx,1)*mfac - size(ii,1)/2; 
  newY = curZ(curlndx,2)*mfac - size(ii,2)/2;

  testX = floor(newX:newX+size(ii,1)); testX(testX<1) = []; testX(testX>size(occupied,1))=[];
  testY = floor(newY:newY+size(ii,1)); testY(testY<1) = []; testY(testY>size(occupied,2))=[];
  oo = occupied(testX,testY);
  if(nnz(oo(:))/numel(oo)>0.4), continue; end;

  if bouts.label(ndx) == 1
    bx_p = [bx_p nan newX newX newX+size(ii,2) newX+size(ii,2) newX];
    by_p = [by_p nan newY newY+size(ii,1) newY+size(ii,1) newY newY];
  else
    bx_n = [bx_n nan newX newX newX+size(ii,2) newX+size(ii,2) newX];
    by_n = [by_n nan newY newY+size(ii,1) newY+size(ii,1) newY newY];
  end

  occupied(testX,testY) = true;
  image(newX,newY,ii,'AlphaData',alpha);
  hszx = size(ii,1) + bottomout; 
  px = [px nan -traj.Y+newX+hszx-rem];
  py = [py nan -traj.X+newY+hszx];
end

phandle = [phandle plot(bx_np,by_np,'-','Color',[1 0.6 0])];
phandle = [phandle plot(bx_nn,by_nn,'-','Color',[0 0.6 1])];
phandle = [phandle plot(bx_p,by_p,'-','Color',[1 0 0])];
phandle = [phandle plot(bx_n,by_n,'-','Color',[0 0 1])];

for pndx = 1:numel(phandle);
  uistack(phandle(pndx),'bottom');
  uistack(phandle(pndx),'up');
end

plot(px,py,'LineStyle','none','Marker','.','MarkerSize',1.5,'MarkerEdgeColor','k');
plot(px+0.5,py,'LineStyle','none','Marker','.','MarkerSize',1.5,'MarkerEdgeColor','k');
plot(px,py+0.5,'LineStyle','none','Marker','.','MarkerSize',1.5,'MarkerEdgeColor','k');
plot(px+0.5,py+0.5,'LineStyle','none','Marker','.','MarkerSize',1.5,'MarkerEdgeColor','k');
axis equal; axis tight; axis off;
% title(sprintf('Train:%d',outqndx));
fname = sprintf('/groups/branson/home/kabram/paperFigs/classifier_chase_new_images_%d',outqndx);
SaveFigLotsOfWays(hfig,fname);
end