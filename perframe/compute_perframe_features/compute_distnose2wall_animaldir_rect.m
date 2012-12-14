% distance to wall
function [data,units] = compute_distnose2wall_animaldir_rect(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  
  xnose = trx(fly).x_mm + 2*trx(fly).a_mm.*cos(trx(fly).theta_mm);
  ynose = trx(fly).y_mm + 2*trx(fly).a_mm.*sin(trx(fly).theta_mm);
  nose = [xnose;ynose];
  body = [trx(fly).x_mm;trx(fly).y_mm ];
  tl = [trx.landmark_params{n}.tl_x(fly),trx.landmark_params{n}.tl_y(fly)];
  bl = [trx.landmark_params{n}.bl_x(fly),trx.landmark_params{n}.bl_y(fly)];
  tr = [trx.landmark_params{n}.tr_x(fly),trx.landmark_params{n}.tr_y(fly)];
  br = [trx.landmark_params{n}.br_x(fly),trx.landmark_params{n}.br_y(fly)];

  
  dtop = getDist(nose,body,tl,tr);
  dleft = getDist(nose,body,bl,tl);
  dright = getDist(nose,body,tr,br);
  dbottom = getDist(nose,body,br,bl);

  data{i} = min([dtop;dleft; dright; dbottom],[],1);

end
units = parseunits('mm');


function d = getDist(a,b,u,v)

x1=a(1,:);
x2=b(1,:);
x3=u(1);
x4=v(1);
y1=a(2,:);
y2=b(2,:);
y3=u(2);
y4=v(2);

ua = ((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3))./((y4-y3)*(x2-x1)-(x4-x3)*(y2-y1));
x = x1 + ua.*(x2 - x1);
y = y1 + ua.*(y2 - y1);

side = sign((x4-x3)*(y1-y3)-(y4-y3)*(x1-x3));
d = sqrt( (x-x1).^2 + (y-y1).^2);
d = side.*d;

% sign of ua tell us if the intersection point was away from the
% the body or not. ua>0 means that the intersection point was towards the
% tail rather than towards the head.
% side tells us whether the nose was beyond the wall or not.
% ua>0 & side>0 means that nose was inside the wall and the intersection
% point was behind the mice. In this case set the dist to 0.
d(ua>0 & side>0) = inf;
