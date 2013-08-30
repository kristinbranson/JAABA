function [binidx, nbins, featurenames, featureboundaries, featurecenters] = compute_spacetime_mask(meana, meanb)

boxwidth2 = round(meanb*8);
boxheight2 = round(meana*4);

nbinsr = 4;
nbinstheta = 8;  grow_theta=1.25;

for j=1:2
  as = linspace(meana,boxheight2,nbinsr);
  bs = linspace(meanb,boxwidth2,nbinsr);
  if(j==2)
    as = as - meana/2;
    bs = bs - meanb/2;
  end
  [xgrid,ygrid] = meshgrid(-boxwidth2:boxwidth2,-boxheight2:boxheight2);
  rbin = zeros(size(xgrid));
  for i = 1:nbinsr,
    idx = rbin == 0 & xgrid.^2/bs(i)^2 + ygrid.^2/as(i)^2 <= 1;
    rbin(idx) = i;
  end

  thetagrid = atan2(ygrid,xgrid);
  thetagrid2 = thetagrid+pi/2;
  idx = thetagrid2>pi;
  thetagrid2(idx) = thetagrid2(idx)-2*pi;
  splits=cumsum(grow_theta.^(0:(nbinstheta/2-1)));
  if(j==2)
    splits=splits/splits(end)*(pi+0.0001);
    thetabin=arrayfun(@(x) find((splits-x)>0,1,'first'), abs(thetagrid2));
    thetabin=thetabin.*((thetagrid2>0)*2-1);
    thetabin(thetabin==-4)=5;
    thetabin(thetabin==-3)=6;
    thetabin(thetabin==-2)=7;
    thetabin(thetabin==-1)=8;
  else
    tmp=splits(end);
    splits=conv([0 splits],[1 1],'valid')/2;
    splits=[splits/tmp*pi inf];
    thetabin=arrayfun(@(x) find((splits-x)>0,1,'first'), abs(thetagrid2));
    thetabin=thetabin.*((thetagrid2>0)*2-1);
    min(min(abs(thetabin)));
    thetabin(thetabin==-ans)=ans;
    max(max(abs(thetabin)));
    thetabin(thetabin==-ans)=ans;
    thetabin(thetabin==-4)=6;
    thetabin(thetabin==-3)=7;
    thetabin(thetabin==-2)=8;
  end

  binidx0 = zeros(size(xgrid));
  binidx0(rbin>0) = sub2ind([nbinsr+1,nbinstheta],rbin(rbin>0),thetabin(rbin>0));
  if(j==1)
    binidx0(rbin==1) = 1;
  else
    binidx0(rbin==1) = 0;
  end
  [binis,~,binidx0] = unique(binidx0);
  binidx0(binidx0==find(binis==0)) = 0;
  binidx0 = binidx0-1;
  binis(binis==0) = [];
  binidx0 = reshape(binidx0,size(xgrid));
  nbins0 = numel(binis);
  binidx{j} = binidx0;
  nbins{j} = nbins0;

  featurenames{j} = cell(1,nbins0);
  for i = 1:nbins0,
    [r,t] = ind2sub([nbinsr+1,nbinstheta],binis(i));
    r = r-1;
    
    %featurenames{i} = sprintf('theta%d_r%d',round((thetalims(1)+dtheta*(t-.5))*180/pi),r);
    featurenames{j}{i} = sprintf('t%d_r%d',t,r);
    if(j==2)
      featurenames{j}{i} = [featurenames{j}{i} '_overlap'];
    end
  end
  if(j==1)
    featurenames{j}{1} = 't0_r0';
  end

  for i = 1:nbins0,
    [featureboundaries{j}{i} featurecenters{j}{i}]=compute_feature_boundaries_and_centers(binidx0,i);
  end
  
  if(j==1)
    binidx{3}=binidx{1};
    binidx{3}((binidx{3}==1)|(binidx{3}==-1))=0;
    nbins{3}=1;
    for r=1:3
      for theta=[1 5]
        idx=find(strcmp(featurenames{j},['t' num2str(theta) '_r' num2str(r)]));
        binidx{3}(binidx{3}==idx)=0;
      end
      for theta=[2 8; 3 7; 4 6]'
        idx=find(strcmp(featurenames{j},['t' num2str(theta(1)) '_r' num2str(r)]));
        idx2=find(strcmp(featurenames{j},['t' num2str(theta(2)) '_r' num2str(r)]));
        nbins{3}=nbins{3}+1;
        binidx{3}(binidx{3}==idx )=-nbins{3};
        binidx{3}(binidx{3}==idx2)=-nbins{3};
        featurenames{3}{nbins{3}-1}=['t' num2str(theta(1)) num2str(theta(2)) '_r' num2str(r) '_symmetric'];
        [featureboundaries{3}{nbins{3}-1} featurecenters{3}{nbins{3}-1}]=...
              compute_feature_boundaries_and_centers(binidx{3},-nbins{3});
      end
    end
    binidx{3}=abs(binidx{3})-1;
    nbins{3}=nbins{3}-1;

    binidx{5}=binidx{1};
    binidx{5}((binidx{5}==1)|(binidx{5}==-1))=0;
    nbins{5}=1;
    for theta=1:8
      idx=find(strcmp(featurenames{j},['t' num2str(theta) '_r1']));
      idx2=find(strcmp(featurenames{j},['t' num2str(theta) '_r2']));
      idx3=find(strcmp(featurenames{j},['t' num2str(theta) '_r3']));
      nbins{5}=nbins{5}+1;
      binidx{5}(binidx{5}==idx )=-nbins{5};
      binidx{5}(binidx{5}==idx2)=-nbins{5};
      binidx{5}(binidx{5}==idx3)=-nbins{5};
      featurenames{5}{nbins{5}-1}=['t' num2str(theta) '_r0_symmetric'];
      [featureboundaries{5}{nbins{5}-1} featurecenters{5}{nbins{5}-1}]=...
            compute_feature_boundaries_and_centers(binidx{5},-nbins{5});
    end
    binidx{5}=abs(binidx{5})-1;
    nbins{5}=nbins{5}-1;
  end

  if(j==2)
    binidx{4}=binidx{2};
    binidx{4}(binidx{4}==-1)=0;
    nbins{4}=1;
    for r=1:3
      for theta=[1 8; 2 7; 3 6; 4 5]'
        idx=find(strcmp(featurenames{j},['t' num2str(theta(1)) '_r' num2str(r) '_overlap']));
        idx2=find(strcmp(featurenames{j},['t' num2str(theta(2)) '_r' num2str(r) '_overlap']));
        nbins{4}=nbins{4}+1;
        binidx{4}(binidx{4}==idx )=-nbins{4};
        binidx{4}(binidx{4}==idx2)=-nbins{4};
        featurenames{4}{nbins{4}-1}=['t' num2str(theta(1)) num2str(theta(2)) '_r' num2str(r) '_symmetric'];
        [featureboundaries{4}{nbins{4}-1} featurecenters{4}{nbins{4}-1}]=...
              compute_feature_boundaries_and_centers(binidx{4},-nbins{4});
      end
    end
    binidx{4}=abs(binidx{4})-1;
    nbins{4}=nbins{4}-1;
  end
end


function [boundary center]=compute_feature_boundaries_and_centers(binidx,idx)

% col = min(find(sum(binidx==idx)));
% row = min(find(binidx(:,col)==idx));
% boundary = bwboundaries(binidx==idx,[row col],'N');  
boundary = bwboundaries(binidx==idx);  
% x = mean(find(sum(binidx==idx,1)));
% y = mean(find(sum(binidx==idx,2)));
x = mean(boundary{1}(:,2));
y = mean(boundary{1}(:,1));
center = [x y];
