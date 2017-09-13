function err = RotateMatrixCriterion(R,U)

dx = R*U(:,1)-[1;0;0;0];
dy = R*U(:,2)-[0;1;0;0];
err = dx'*dx + dy'*dy;