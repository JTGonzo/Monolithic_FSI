function [DeEF,Sol]=feue(ALE,DeEF,ElemF,ElemS,Phy,Sol)
% u-eta formulation
global GM
% First, determine dof which is the dof at t_n+1. 
% Calculate md which is the displacement from dof1 to dof
% me is the node displacement from Elem.dof to DeE{1}.dof

% Sol.p is the p at t_{n+1}. 
% Sol.U(:,:,[1,2]) are the U at t_n+1 and t_n

me=zeros(ElemF.nDof,2);
me(ElemF.itfNode,:)=Sol.eta(ElemS.itfNode,:);
tempInd=[ALE.fFreeN;ALE.fFreeN+ALE.nVertex];
b=-ALE.fixedS(tempInd,:)*[me(1:ALE.nVertex,1,1);me(1:ALE.nVertex,2,1)];
me(ALE.fFreeN,:)=reshape(ALE.fixedS(tempInd,tempInd)\b,[],2);

DeEF{1}.dof=completeDof(ElemF.dof+me,ElemF); % input: dofs on bdy are right
DeEF{1}=deformedMeshData(DeEF{1}.dof,ElemF); % dof at t_n+1
if min(DeEF{1}.detT)<0, disp('detT<0 in fl'); pause; end

md=DeEF{1}.dof-DeEF{2}.dof;
mv=md./Sol.delt;
mv_qt=uGraduAtQuadPts(mv,DeEF{2},ElemF); % the only place that uses DeE{2}


if Sol.NSE~=2 % stokes or u^n grad u^n or L^n+1-L^n
    [fM,fC] = initMatrix2(DeEF{1},ElemF); % t_n+1 dof
    GM.fM=fM;
    fM=Phy.frho*fM;
else % -1/2 div u^n u^n+1 + (u^n - mv) \grad u^n+1
    [fM,fC,fF] = initMatrix4(Sol.U_qt(:,:,:,2)-mv_qt,...
        Sol.Ux_qt(:,:,:,2),Sol.Uy_qt(:,:,:,2),DeEF{1},ElemF);
    GM.fM=fM;
    fM=Phy.frho*fM+Phy.frho*Sol.delt*fF;
end

[UxP,UyP] = initMatrix3(DeEF{1},ElemF); % t_n+1 dof
UxPUyP=(-Sol.delt)*[UxP;UyP]; % times -delt

zero1=sparse(ElemF.nDof,ElemF.nDof);
zero2=sparse(ElemF.nVertex,ElemF.nVertex);
if Sol.NSE==3
    coef=Phy.frho*(Phy.fmu+Phy.fmu1)*Sol.delt;
else
    coef=Phy.frho*Phy.fmu*Sol.delt;
end
fE=[[fM,zero1;zero1,fM]+coef*fC,UxPUyP;UxPUyP',zero2];

gg=repmat(Sol.delt*Phy.g*[sin(Phy.theta),-cos(Phy.theta)],ElemS.nDof,1);
bs=GM.sM*(Sol.etat(:,:,2)+gg); % sM contains rho

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
    case 2
        S=initMatrixSVK2(ElemS,Phy,Sol.eta);
        zero1=sparse(ElemS.nDof,ElemS.nDof);
        sE=[GM.sM,zero1;zero1,GM.sM]+Sol.delt^2*S; % delt^2
        GM.sE=sparse(GM.nd,GM.nd);
        GM.sE(GM.s2g,GM.s2g)=sE;
        GM.sEItf=sE(GM.sEItfIx,GM.sEItfIx);
        
        bs=bs(:)-Sol.delt*(S*(Sol.eta(:)+ElemS.dof(:)));
end

% E will be a global matrix for (u,v,p,eta1,eta2)
% E=GM.sE; 
% E(1:GM.nuvp,1:GM.nuvp)=fE; % this line still takes the most of the time
E=[fE,GM.sE(1:GM.nuvp,GM.nuvp+1:GM.nd);GM.sE(GM.nuvp+1:GM.nd,:)]; 
E(GM.overlapIx,GM.overlapIx)=E(GM.overlapIx,GM.overlapIx)+GM.sEItf;

gg=Sol.delt*Phy.g*[sin(Phy.theta),-cos(Phy.theta)];
temp=Sol.U_qt(:,:,:,2);
temp(:,:,1)=temp(:,:,1)+gg(1);
temp(:,:,2)=temp(:,:,2)+gg(2);

if Sol.NSE==0 % stokes
    % mv_qt(:)=0; % even set mv=0 to check the stability
    fU_ugradu=nonlin_u(-mv_qt,Sol.Ux_qt(:,:,:,2),...
        Sol.Uy_qt(:,:,:,2),DeEF{1},ElemF); % t_n+1 dof
    bf=Phy.frho*(linear_u(temp,DeEF{1},ElemF)-Sol.delt*fU_ugradu);
elseif Sol.NSE==1 % u^n grad u^n
    fU_ugradu=nonlin_u(Sol.U_qt(:,:,:,2)-mv_qt,Sol.Ux_qt(:,:,:,2),...
        Sol.Uy_qt(:,:,:,2),DeEF{1},ElemF); % t_n+1 dof
    bf=Phy.frho*(linear_u(temp,DeEF{1},ElemF)-Sol.delt*fU_ugradu);    
elseif Sol.NSE==2 % -1/2 div u^n u^n+1 + (u^n - mv) \grad u^n+1
    bf=Phy.frho*linear_u(temp,DeEF{1},ElemF);
elseif Sol.NSE==3
    fU_ugradu=nonlin_u(Sol.U_qt(:,:,:,2)-mv_qt,Sol.Ux_qt(:,:,:,2),...
        Sol.Uy_qt(:,:,:,2),DeEF{1},ElemF); % t_n+1 dof
    bf=Phy.frho*(linear_u(temp,DeEF{1},ElemF)-Sol.delt*fU_ugradu);
    bf=bf(:)+(Phy.fmu1*Phy.frho*Sol.delt)*(fC*[Sol.U(:,1,2);Sol.U(:,2,2)]);
else
    disp('not implemented'); pause; 
end

if strcmp(ElemF.meshString(1:4),'FPCf')==1
    bfOut=uvTractionBC(Sol.delt,DeEF{1},ElemF,Phy); % with sig n given
    bf=bf(:)+bfOut(:);
end

bb=zeros(GM.nd,1);
bb(1:2*ElemF.nDof)=bf(:);
bb(GM.s2g)=bb(GM.s2g)+bs;
Ubd=Ubdry(DeEF{1}.dof(ElemF.bdDirN,:),Sol.t);
bb=bb-E(:,GM.bdDirNf)*Ubd(:);

uvpe=zeros(GM.nd,1);
uvpe(GM.uvpeFreeN)=E(GM.uvpeFreeN,GM.uvpeFreeN)\bb(GM.uvpeFreeN);

Sol.U(:,:,1)=reshape(uvpe(1:2*ElemF.nDof),[],2);
Sol.U(ElemF.bdDirN,:,1)=Ubd;
Sol.p=uvpe(GM.nuv+1:GM.nuvp);
Sol.etat(:,:,1)=reshape(uvpe(GM.s2g),[],2);