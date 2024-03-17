function [DeE,Sol]=memvDeE(mi,ALE,Elem,Sol)
% move the mesh.
me=zeros(Elem.nDof,2); mv=me;

% displacement at t_n+1 (next instant)
me(Elem.itfNode,:)=BDdv(Elem.dof(Elem.itfNode,:),Sol.t+Sol.delt,'d');

% velocity at t_n (right now)
mv(Elem.itfNode,:)=BDdv(Elem.dof(Elem.itfNode,:),Sol.t,'v');

tempInd=[ALE.fFreeN;ALE.fFreeN+ALE.nVertex];
be=-ALE.fixedS(tempInd,:)*[me(1:ALE.nVertex,1);me(1:ALE.nVertex,2)];
bv=-ALE.fixedS(tempInd,:)*[mv(1:ALE.nVertex,1);mv(1:ALE.nVertex,2)];
temp=reshape(ALE.fixedS(tempInd,tempInd)\[be,bv],[],4);
me(ALE.fFreeN,:)=temp(:,1:2);
mv(ALE.fFreeN,:)=temp(:,3:4);
DeE.dof=completeDof(Elem.dof+me,Elem); % input: dofs on bdry are right
Sol.mv_qt(:,:,:,mi)=uGraduAtQuadPts(completeDof(mv,Elem),DeE,Elem);
DeE=deformedMeshData(DeE.dof,Elem);
