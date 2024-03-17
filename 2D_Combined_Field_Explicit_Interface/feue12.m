function [DeEF,Sol]=feue2(ALE,DeEF,ElemF,ElemS,Phy,Sol)
% u-eta formulation
% conservative formulation
global GM
% First, determine dof which is the dof at t_n+1. 
% me is the displacement
% Sol.p is the p at t_{n+1}. 
% Sol.U(:,:,[1,2]) are the U at t_n+1 and t_n
% Sol.me(:,:,[1,2]) are the mesh displacement at t_n+1 and t_n
% By taking the difference, we obtain Sol.mv which is defined on t_n+1 mesh
% DeE is the current (t_n+1) mesh information (deformed element)
% MiE is the mesh information at t_{n+1/2} (middle element)

% eta is displacement
Sol.eta=Sol.eta+Sol.delt*Sol.etat(:,:,1);

me=zeros(ElemF.nDof,2);
me(ElemF.itfNode,:)=Sol.eta(ElemS.itfNode,:); % me is displacement
tempInd=[ALE.fFreeN;ALE.fFreeN+ALE.nVertex];
b= -ALE.fixedS(tempInd,:)*[me(1:ALE.nVertex,1);me(1:ALE.nVertex,2)]; % x by 2N and 2N by 1 where x is the number of rows in tempInd
me(ALE.fFreeN,:)=reshape(ALE.fixedS(tempInd,tempInd)\b,size(ALE.fFreeN,1),2);

DeEF.dof=completeDof(me+ElemF.dof,ElemF); % input: dofs on bdy are right
% now me is current position

DeEF=deformedMeshData(DeEF.dof,ElemF); % dof at t_n+1
MiEF=deformedMeshData(0.5*(Sol.me(:,:,1)+Sol.me(:,:,2))+ElemF.dof,ElemF); % mid mesh
if min(DeEF.detT)<0, disp('detT<0 in feue2'); pause; end

% fluid matrices
[GM.fS,GM.fM,fC] = initMatrix2old(DeEF,ElemF); % t_n+1 dof
fM=Phy.frho*GM.fM; % fM contains rho

fR=(-Phy.frho*Sol.delt)*initMatrixdivuwv(MiEF,ElemF,Sol); % on mid mesh
[UxP,UyP] = initMatrix3(DeEF,ElemF); % t_n+1 dof
UxPUyP=(-Sol.delt)*[UxP;UyP]; % times -delt

zero1=sparse(ElemF.nDof,ElemF.nDof);
zero2=sparse(ElemF.pnDof,ElemF.pnDof);
coef=Phy.frho*Phy.fmu*Sol.delt;
fE=[[fM,zero1;zero1,fM]+coef*fC+fR,UxPUyP;UxPUyP',zero2];

% fluid forcing term
gg=Phy.g*[sin(Phy.theta),-cos(Phy.theta)];
Mg=fM*repmat(Sol.delt*gg,ElemF.nDof,1); % fM contains rho, delt comes from the time derivative...this term represents <rho*g,phi>, so we can't have a fM over here
bf=Sol.MU+Mg;    % the MU term comes from the time derivative...
if Phy.g~=0 && ~isempty(ElemF.openbdryElem) % strcmp(ElemF.meshString(1:4),'FPCf')==1
    bfOut=uvTractionBC(Sol.delt,DeEF,ElemF,Phy); 
    bf=bf(:)+bfOut(:);
end

% solid forcing term
bs=GM.sM*(Sol.etat(:,:,2)+repmat(Sol.delt*gg,ElemS.nDof,1)); % sM contains rho

switch Sol.SVK
    case 0
        bs=bs(:)-GM.sL*Sol.eta(:); % sL contains delt
    case 1
        [S,bsvk]=initMatrixSVK(ElemS,Phy,Sol.eta);
        zero1=sparse(ElemS.nDof,ElemS.nDof);
        sE=[GM.sM,zero1;zero1,GM.sM]+Sol.delt^2*S; % delt^2
        GM.sE=sparse(GM.nd,GM.nd);
        GM.sE(GM.s2g,GM.s2g)=sE;
        GM.sEItf=sE(GM.sEItfIx,GM.sEItfIx);
        
        bs=bs(:)-Sol.delt*bsvk(:);
end

% E will be a global matrix for (u,v,p,eta1,eta2)
% E=GM.sE; 
% E(1:GM.nuvp,1:GM.nuvp)=fE; % this line still takes the most of the time
E=[fE,GM.sE(1:GM.nuvp,GM.nuvp+1:GM.nd);GM.sE(GM.nuvp+1:GM.nd,:)]; 
E(GM.overlapIx,GM.overlapIx)=E(GM.overlapIx,GM.overlapIx)+GM.sEItf;


bb=zeros(GM.nd,1);
bb(1:2*ElemF.nDof)=bf(:);
bb(GM.s2g)=bb(GM.s2g)+bs;
Ubd=Ubdry(DeEF.dof(ElemF.bdDirN,:),ElemF.meshString,Sol.t);
bb=bb-E(:,GM.bdDirNf)*Ubd(:);

uvpe=zeros(GM.nd,1);
uvpe(GM.uvpeFreeN)=E(GM.uvpeFreeN,GM.uvpeFreeN)\bb(GM.uvpeFreeN);

% save the old data
Sol.etat(:,:,2)=Sol.etat(:,:,1);
Sol.U(:,:,2)=Sol.U(:,:,1);
Sol.me(:,:,2)=Sol.me(:,:,1);

% store the new data
Sol.me(:,:,1)=DeEF.dof-ElemF.dof;
Sol.U(:,:,1)=reshape(uvpe(1:2*ElemF.nDof),[],2);
Sol.U(ElemF.bdDirN,:,1)=Ubd;
Sol.p=uvpe(GM.nuv+1:GM.nuvp);
Sol.etat(:,:,1)=reshape(uvpe(GM.s2g),[],2);
Sol.MU=fM*Sol.U(:,:,1); % contains rho

% error1 = max(abs(Sol.U(:,1,1)-Sol.ActU));
error1 = 0;
Sol.error = [Sol.error;Sol.t error1];