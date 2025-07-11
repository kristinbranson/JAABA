function [mdist] = malah_dist(x,y,a,b,theta,xLD,yLD )

%compute mahalanobis distance between an apt landmark of toucher fly1 and
%the ctrax ellipse of the touchee fly2

%compute covariance mat of fly2
N = numel(theta);
cost = cos(theta);
sint = sin(theta);

R(1,1,:) = cost;
R(1,2,:) = -sint;
R(2,1,:) = sint;
R(2,2,:) = cost;

D(1,1,:) = a.^2;
D(1,2,:) = zeros(1,N);
D(2,1,:) = zeros(1,N);
D(2,2,:) = b.^2;

Sigma = zeros(2,2,N);
for i= 1:N
    Sigma(:,:,i) = R(:,:,i) * D(:,:,i) * R(:,:,i)';
end



% offset of leg point (LD - ctr) 
xoff = xLD'-x;
yoff = yLD'-y;
mdist = zeros(N,1);
for i = 1:N
mdist(i) = sqrt([xoff(i),yoff(i)]*inv(Sigma(:,:,i)) * [xoff(i),yoff(i)]');
end 
