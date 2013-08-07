function ret_val = compute_spacetime_gradient(imnorm,imnorm_last,binidx,nbins,featurenames,dt)

imdiff=abs(imnorm-imnorm_last);
meanval = accumarray(binidx(binidx>0),imdiff(binidx>0),[nbins,1],@mean);
%   meanim = zeros(size(imnorm));
%   for i = 1:nbins,
%     meanim(binidx==i) = meanval(i);
%   end
%   subplot(1,4,4);
%   imagesc(meanim,[0,1]);
%   axis image;

%data{i1}(framei-firstframes(i1)+1,:) = meanval./dt0{fly1}(framei-firstframes(i1)+1);
ret_val = meanval./dt;

for theta=[2 8; 3 7; 4 6]'
  for r=1:3
    idx=find(strcmp(featurenames,['theta' num2str(theta(1)) '_r' num2str(r)]));
    idx2=find(strcmp(featurenames,['theta' num2str(theta(2)) '_r' num2str(r)]));
    ret_val(end+1) = (ret_val(idx) + ret_val(idx2))/2;
  end
end

for theta=[4 2; 5 1; 6 8]'
  for r=1:3
    idx=find(strcmp(featurenames,['theta' num2str(theta(1)) '_r' num2str(r)]));
    idx2=find(strcmp(featurenames,['theta' num2str(theta(2)) '_r' num2str(r)]));
    ret_val(end+1) = (ret_val(idx) + ret_val(idx2))/2;
  end
end

for theta=1:8
  idx=find(strcmp(featurenames,['theta' num2str(theta) '_r1']));
  idx2=find(strcmp(featurenames,['theta' num2str(theta) '_r2']));
  idx3=find(strcmp(featurenames,['theta' num2str(theta) '_r3']));
  ret_val(end+1) = (ret_val(idx) + ret_val(idx2) + ret_val(idx3))/3;
end