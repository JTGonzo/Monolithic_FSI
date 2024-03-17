function figPlot(k,p,pex,dof,fn,str1,str2,gamI)
DelTri = delaunay(dof(:,1),dof(:,2));
if ~isempty(gamI.ellp)
    DelTri = removeHoles(DelTri,dof,gamI.rect,gamI.ellp{1});
end
if ~isempty(pex)
    % figure(2); plot(p-fn*pex(dof(:,1),dof(:,2)));
%     figure(3); showmesh(dof,DelTri);
%     figure(4);showsolution(dof,DelTri,p); title(str1);
%     figure(5);showsolution(dof,DelTri,fn*pex(dof(:,1),dof(:,2))); title(str2);
    figure(k);showsolution(dof,DelTri,p-fn*pex(dof(:,1),dof(:,2))); title([str1,'-',str2]);
    % view(2); axis equal; axis off;
else
    tricontour(dof,DelTri,p,20);
    axis equal; axis off;
end

