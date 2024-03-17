function DeE=deformedMeshData(dof,Elem)
DeE.dof=dof;
% [DeE.bdNEInvDF,DeE.bdNEAbsDF]=generateDF...
%     (Elem.elem2dof,Elem.bdNeuElem,dof,Elem.Quad.eQuad,Elem.Quad.bdQuad);
% [DeE.bdDEInvDF,DeE.bdDEAbsDF]=generateDF...
%     (Elem.elem2dof,Elem.bdDirElem,dof,Elem.Quad.eQuad,Elem.Quad.bdQuad);
% [DeE.itfInvDF,DeE.itfAbsDF]=generateDF...
%     (Elem.elem2dof,Elem.itfElem,dof,Elem.Quad.eQuad,Elem.Quad.bdQuad);

[DeE.bdInvDF,DeE.bdAbsDF]=generateDF...
    (Elem.elem2dof,Elem.bdElem,dof,Elem.Quad.eQuad,Elem.Quad.bdQuad);
% DeE.bdAbsDF=[DeE.bdNEAbsDF,DeE.bdDEAbsDF];

% p1p2 and p3p1 are calculated based on dof and elem2dof
DeE.p1p2 = dof(Elem.elem2dof(Elem.inElem,2),:)...
    -dof(Elem.elem2dof(Elem.inElem,1),:); % edge 12
DeE.p3p1 = dof(Elem.elem2dof(Elem.inElem,1),:)...
    -dof(Elem.elem2dof(Elem.inElem,3),:); % edge 31
% det only. No abs
% after reordering the nodes on each triangle, it is positive
DeE.detT = DeE.p1p2(:,2).*DeE.p3p1(:,1)-DeE.p1p2(:,1).*DeE.p3p1(:,2);

