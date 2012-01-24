function [rfrac,r,isonfly] = center_of_rotation(trk,singleframe,debug)

if ~exist('singleframe','var')
  singleframe = false;
end

cost = cos(trk.theta);
sint = sin(trk.theta);
if singleframe,
  dcos = diff(cost);
  dsin = diff(sint);
  dx = diff(trk.x);
  dy = diff(trk.y);
else
  dcos = [(cost(2)-cost(1))*2,cost(3:end)-cost(1:end-2),(cost(end)-cost(end-1))*2];
  dsin = [(sint(2)-sint(1))*2,sint(3:end)-sint(1:end-2),(sint(end)-sint(end-1))*2];
  dx = [(trk.x(2)-trk.x(1))*2,trk.x(3:end)-trk.x(1:end-2),(trk.x(end)-trk.x(end-1))*2];
  dy = [(trk.y(2)-trk.y(1))*2,trk.y(3:end)-trk.y(1:end-2),(trk.y(end)-trk.y(end-1))*2];
end
r = -( dx.*dcos + dy.*dsin )./ (dcos.^2 + dsin.^2);
if singleframe,
  a = (trk.a(1:end-1) + trk.a(2:end))/2;
  rfrac = r./a/2;
else
  rfrac = r./trk.a/2;
end
% try the end points if rfrac > 1, < -1
badidx = abs(rfrac) > 1;
isonfly = ~badidx;
% if any(badidx),
%   errplus = (dx(badidx) + 2*trk.a(badidx).*dcos(badidx)).^2 + ...
%     (dy(badidx) + 2*trk.a(badidx).*dsin(badidx)).^2;
%   errminus = (dx(badidx) - 2*trk.a(badidx).*dcos(badidx)).^2 + ...
%     (dy(badidx) - 2*trk.a(badidx).*dsin(badidx)).^2;
%   newrfrac = ones(1,length(errplus));
%   newrfrac(errplus>errminus) = -1;  
%   rfrac(badidx) = newrfrac;
%   r(badidx) = newrfrac.*trk.a(badidx)*2;
% end

if exist('debug','var') && debug,
  if singleframe,
    n = trk.nframes-1;
  else
    n = trk.nframes;
  end
  bestremp = zeros(1,n);
  for i = 1:n,
    rtry = linspace(-1,1,101).*trk.a(i)*2;
    err = (dx(i) + rtry.*dcos(i)).^2 + (dy(i) + rtry.*dsin(i)).^2;
    bestremp(i) = rtry(argmin(err));
  end
  figure(1234);
  clf;
  plot(r,bestremp,'.');
  xlabel('analytically found best r');
  ylabel('empirically found best r');
  axis equal;
  hold on;
  tmp1 = min(min(r),min(bestremp));
  tmp2 = max(max(r),max(bestremp));
  plot([tmp1,tmp2],[tmp1,tmp2],'r');
end