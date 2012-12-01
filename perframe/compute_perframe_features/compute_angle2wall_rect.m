% distance to wall
function [data,units] = compute_angle2wall_rect(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  x = trx(fly).x_mm;
  y = trx(fly).y_mm;
  
  walla = zeros(1,4);
  walla(1) = atan2(trx.landmark_params{n}.tr_y(fly)-trx.landmark_params{n}.tl_y(fly),trx.landmark_params{n}.tr_x(fly)-trx.landmark_params{n}.tl_x(fly));
  walla(2) = atan2(trx.landmark_params{n}.br_y(fly)-trx.landmark_params{n}.tr_y(fly),trx.landmark_params{n}.br_x(fly)-trx.landmark_params{n}.tr_x(fly));
  walla(3) = atan2(trx.landmark_params{n}.bl_y(fly)-trx.landmark_params{n}.br_y(fly),trx.landmark_params{n}.bl_x(fly)-trx.landmark_params{n}.br_x(fly));
  walla(4) = atan2(trx.landmark_params{n}.tl_y(fly)-trx.landmark_params{n}.bl_y(fly),trx.landmark_params{n}.tl_x(fly)-trx.landmark_params{n}.bl_x(fly));
  
  
  [dtop] = getDist(trx.landmark_params{n}.tl_x(fly),trx.landmark_params{n}.tl_y(fly), trx.landmark_params{n}.tr_x(fly),trx.landmark_params{n}.tr_y(fly),...
    x,y);
  [dright]= getDist(trx.landmark_params{n}.tr_x(fly),trx.landmark_params{n}.tr_y(fly), trx.landmark_params{n}.br_x(fly),trx.landmark_params{n}.br_y(fly),...
    x,y);
  [dbottom] = getDist(trx.landmark_params{n}.bl_x(fly),trx.landmark_params{n}.bl_y(fly), trx.landmark_params{n}.br_x(fly),trx.landmark_params{n}.br_y(fly),...
    x,y);
  [dleft] = getDist(trx.landmark_params{n}.tl_x(fly),trx.landmark_params{n}.tl_y(fly), trx.landmark_params{n}.bl_x(fly),trx.landmark_params{n}.bl_y(fly),...
    x,y);
  
  [~,closest] = min([dtop;dright; dbottom; dleft],[],1);
  
  closest_walla = walla(closest);
  targeta = trx(fly).theta;
  
  a = mod(targeta-closest_walla+pi/2+pi,2*pi)-pi;
  
  data{i} = a;
end
units = parseunits('rad');


function [d,a] = getDist(p1_x,p1_y,p2_x,p2_y,p_x,p_y)

dp1p2 = (p1_x-p2_x).^2 + (p1_y-p2_y).^2;


dotpr_x = (p1_x - p_x).*(p1_x-p2_x); 
dotpr_y = (p1_y - p_y).*(p1_y-p2_y);

t = (dotpr_x+dotpr_y)/dp1p2;
proj_x = p1_x + t.*(p2_x-p1_x);
proj_y = p1_y + t.*(p2_y-p1_y);
p2proj_x = p_x-proj_x;
p2proj_y = p_y-proj_y;

d = sqrt(  (p2proj_x).^2 + (p2proj_y).^2);
