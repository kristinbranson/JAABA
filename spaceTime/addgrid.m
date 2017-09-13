function I = addgrid(I,gsz)

gridpts = 1:gsz(2):size(I,2);
grid_x(:,1) = gridpts;
grid_x(:,3) = gridpts;
grid_x(:,2) = 1;
grid_x(:,4) = size(I,1);
gridpts = 1:gsz(1):size(I,1);
grid_y(:,2) = gridpts;
grid_y(:,4) = gridpts;
grid_y(:,1) = 1;
grid_y(:,3) = size(I,2);

I = insertShape(I,'Line',[grid_x;grid_y]);