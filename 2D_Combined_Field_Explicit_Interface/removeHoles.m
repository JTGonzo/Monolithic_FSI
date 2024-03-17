function tri = removeHoles(tri,dof,rect,ellp);
if ~isempty(ellp)
    cent=1/3*(dof(tri(:,1),:)+dof(tri(:,2),:)+dof(tri(:,3),:));
    tri=tri(fd(cent,rect,ellp)<0,:);
end

function d=fd(p,rect,ellp)
d=max(drectangle(p,rect(1),rect(2),rect(3),rect(4)),...
      -dist_ellipse(p,ellp(1),ellp(2),ellp(3),ellp(4)));