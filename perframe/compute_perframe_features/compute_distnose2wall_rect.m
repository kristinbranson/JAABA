% distance to wall
function [data,units] = compute_distnose2wall_rect(trx,n)

flies = trx.exp2flies{n};
nflies = numel(flies);
data = cell(1,nflies);
for i = 1:nflies,
  fly = flies(i);  
  
  xnose = trx(fly).x_mm + 2*trx(fly).a_mm.*cos(trx(fly).theta_mm);
  ynose = trx(fly).y_mm + 2*trx(fly).a_mm.*sin(trx(fly).theta_mm);

  
  % The ordering of the arguments is to account for weird y directions of images. 
  dtop = getDist(trx.landmark_params{n}.tl_x(fly),trx.landmark_params{n}.tl_y(fly), trx.landmark_params{n}.tr_x(fly),trx.landmark_params{n}.tr_y(fly),...
    xnose,ynose);
  dleft = getDist(trx.landmark_params{n}.bl_x(fly),trx.landmark_params{n}.bl_y(fly), trx.landmark_params{n}.tl_x(fly),trx.landmark_params{n}.tl_y(fly),...
    xnose,ynose);
  dright = getDist(trx.landmark_params{n}.tr_x(fly),trx.landmark_params{n}.tr_y(fly), trx.landmark_params{n}.br_x(fly),trx.landmark_params{n}.br_y(fly),...
    xnose,ynose);
  dbottom = getDist(trx.landmark_params{n}.br_x(fly),trx.landmark_params{n}.br_y(fly), trx.landmark_params{n}.bl_x(fly),trx.landmark_params{n}.bl_y(fly),...
    xnose,ynose);
  
  data{i} = min([dtop;dleft; dright; dbottom],[],1);

end
units = parseunits('mm');


function d = getDist(p1_x,p1_y,p2_x,p2_y,p_x,p_y)

dp1p2 = (p1_x-p2_x).^2 + (p1_y-p2_y).^2;


dotpr_x = (p1_x - p_x).*(p1_x-p2_x); 
dotpr_y = (p1_y - p_y).*(p1_y-p2_y);

t = (dotpr_x+dotpr_y)/dp1p2;
proj_x = p1_x + t.*(p2_x-p1_x);
proj_y = p1_y + t.*(p2_y-p1_y);
d = sqrt(  (p_x-proj_x).^2 + (p_y-proj_y).^2);

% Find whether the mice is in the interior of the square.
side = sign((p2_x-p1_x)*(p_y-p1_y)-(p2_y-p1_y)*(p_x-p1_x));
d = side.*d;