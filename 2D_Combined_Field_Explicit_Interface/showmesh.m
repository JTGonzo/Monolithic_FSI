function showmesh(dof,elem2dof)
trisurf(elem2dof,dof(:,1),dof(:,2),zeros(size(dof,1),1))
view(2); axis equal; axis off;