function [ElemF,ElemS,Sol,DeES]=fixedData(meshString,nGlobalRefine,Phy,Sol)
% u-v formulation

global ALE GM %contains sM sL sE sEbdry fFixS

% ALE is for fluid phase only
ALE=MeshS(meshString(1,:),nGlobalRefine);

deg=2; % P2/P1 and P2
ElemF=fixedMeshData(deg,meshString(1,:),nGlobalRefine);
ElemS=fixedMeshData(deg,meshString(2,:),nGlobalRefine);

% interface dof
% itfNode = A(I) and A = itfNode(J)
[ElemF.itfNode,II,JJ]=unique(ElemF.elem2dof(ElemF.itfElem,ElemF.indLBN)); 

A=ElemS.elem2dof(ElemS.itfElem,ElemS.indLBN);
A=A(:);A=A(end:-1:1);
ElemS.itfNode=A(II);

if max(max(abs(ElemF.dof(ElemF.itfNode,:)-ElemS.dof(ElemS.itfNode,:))))>5e-11
    disp('mistake in fixedData');
    pause; return;
else
    ElemF.dof(ElemF.itfNode,:)=ElemS.dof(ElemS.itfNode,:);
end

DeES=deformedMeshData(ElemS.dof,ElemS);

% p=ElemF.dof(ElemF.freeN,:);
% plot(p(:,1),p(:,2),'r*');
% hold on;
% showmesh(ElemF.dof,ElemF.elem2dof(:,1:3))
% p=ElemS.dof(ElemS.freeN,:);
% hold on;
% plot(p(:,1),p(:,2),'bo');
% hold on; 
% showmesh(ElemS.dof,ElemS.elem2dof(:,1:3))

% GM.uveFreeN=[ElemF.freeN;ElemS.freeN+ElemF.nDof];

[sM,sC,sD] = initMatrix2(DeES,ElemS);
sM=Phy.srho*sM;
% BDF-1 w. trapezoid quadrature [dt/2,dt/2]
sL=Sol.delt*Phy.smu*sC+Sol.delt*Phy.slam*sD; % contains delt factor
clear sC; clear sD;
zero1=sparse(ElemS.nDof,ElemS.nDof);
sE=[sM,zero1;zero1,sM]+0.5*sL; 

% from solid dof to global dof
% t1(i) is the global index of the ith dof in solid structure
t1=(1:2*ElemS.nDof)';
t2=setdiff(t1,[ElemS.itfNode;ElemS.nDof+ElemS.itfNode]);
t1(ElemS.itfNode)=ElemF.itfNode;
t1(ElemS.nDof+ElemS.itfNode)=ElemF.nDof+ElemF.itfNode;
t1(t2)=(ElemF.nDof*2+ElemF.nVertex)+(1:length(t2));

GM.s2g=t1;
GM.nd=max(t1); % number of uvpe dof
GM.sE=sparse(GM.nd,GM.nd);
GM.sE(t1,t1)=sE;
GM.sM=sM;
GM.sL=sL;

% t3=setdiff(ElemS.bdDirElem,ElemS.itfElem);
% t4=unique(ElemS.elem2dof(t3,ElemS.indLBN)); % include the end pts
t4=setdiff(ElemS.bdDirN,ElemS.itfNode); % exclude the end pts
% this one also include the fixed solid boundry.
GM.bdDirN=[ElemF.bdDirN;ElemF.nDof+ElemF.bdDirN;t1(t4);t1(ElemS.nDof+t4)];
GM.bdDirNf=[ElemF.bdDirN;ElemF.nDof+ElemF.bdDirN];

t5=setdiff((1:ElemF.nDof)',ElemF.bdDirN);
t6=setdiff((1:2*ElemS.nDof)',[ElemS.bdDirN;ElemS.nDof+ElemS.bdDirN]);
GM.uvpeFreeN=[t5;ElemF.nDof+t5;2*ElemF.nDof+(1:ElemF.nVertex)';t1(t6)];

% n1=length(GM.uvpeFreeN);
% n2=2*ElemF.nDof+ElemF.nVertex+2*ElemS.nDof;
% n3=2*length(ElemF.itfNode);
% n4=2*length(ElemF.bdDirN)+2*length(t4);
% % n2-n3-n4-n1==0
 
t7=[ElemS.itfNode;ElemS.nDof+ElemS.itfNode];
GM.sEItf=sE(t7,t7);

GM.overlapIx=[ElemF.itfNode;ElemF.nDof+ElemF.itfNode];
GM.nuv=2*ElemF.nDof;
GM.nuvp=2*ElemF.nDof+ElemF.nVertex;