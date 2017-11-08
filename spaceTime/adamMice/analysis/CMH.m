function [chiMH,p] = CMH(m)
% Cochran-Mantel-Haenszel test for repeated tests of independence
% Ref: www.biostathandbook.com/cmh.html
% 
% Sample data from Ref:
% 
% m = [];
% m(:,:,1) = [56 69;40 77];
% m(:,:,2) = [61 257;57 301];
% m(:,:,3) = [73 65;71 79];
% m(:,:,4) = [71 48;55 48];
% [chiMH,p] = CMH(m)
%
% m = [];
% m(:,:,1) = [2 46;11 41];
% m(:,:,2) = [4 67;12 60];
% m(:,:,3) = [1 86;4  76];
% m(:,:,4) = [1 37;6  32];
% m(:,:,5) = [2 92;1  93];
% [chiMH,p] = CMH(m)

assert(isnumeric(m) && size(m,1)==2 && size(m,2)==2);

a = m(1,1,:);
b = m(1,2,:);
c = m(2,1,:);
d = m(2,2,:);
n = a+b+c+d;
a_exp = (a+b).*(a+c)./n;
a_var = (a+b).*(a+c).*(b+d).*(c+d)./(n.^3-n.^2);

num = (abs(sum(a-a_exp))-0.5)^2; % includes continuity correction
den = sum(a_var);
chiMH = num/den;

p = 1-chi2cdf(chiMH,1);
