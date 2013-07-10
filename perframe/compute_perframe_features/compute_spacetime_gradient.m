function ret_val = compute_spacetime_gradient(imnorm,imnorm_last,binidx,nbins,dt)

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
