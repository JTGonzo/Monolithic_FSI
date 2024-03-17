function [DeEF,Sol]=feue2Predictor2(ALE,DeEF,ElemF,ElemS,Phy,Sol,meshString,etat)
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

eta=Sol.eta+Sol.delt*etat;
me=zeros(ElemF.nDof,2);
me(ElemF.itfNode,:)=eta(ElemS.itfNode,:); % me is displacement
tempInd=[ALE.fFreeN;ALE.fFreeN+ALE.nVertex];
b1=-ALE.fixedS(tempInd,:)*[me(1:ALE.nVertex,1);me(1:ALE.nVertex,2)];%-N(tempInd,:)*[me(1:ALE.nVertex,1);me(1:ALE.nVertex,2)]; 
a = ALE.fixedS(tempInd,tempInd)...
    \(b1);%+N(tempInd,:)*[me0(1:ALE.nVertex,1);me0(1:ALE.nVertex,2)]);
me(ALE.fFreeN,:) = reshape(a,length(ALE.fFreeN),2);

dof = me+Sol.dof;

%% being added on 29th Dec 2012 
%% Last Modified on 2nd Jan 2013
% This being added to smoothen the mesh
% this is an implementation simple Lagragian Smoothing
%% Code starts Here

% dof(ALE.fFreeNodes,:) = ElemF.Node2Node(ALE.fFreeNodes,:)*dof(1:ALE.nVertex,:)./(sum(ElemF.Node2Node(ALE.fFreeNodes,:),2)*[1 1]);

%% Code Added on 31st Dec 2012
DeEF.angle = getAngle(dof,ElemF);

if min(min(DeEF.angle))<15 || max(max(DeEF.angle))>150
    min(min(DeEF.angle))
    disp('angle <15');

end
%% Code Ends Here for 30th Dec 2012 Last Modified 31st Dec 2012
%% Code End Here for 29th Dec 2012

DeEF.dof = completeDof(dof,ElemF);
% dof_nxt = completeDof(dof_nxt,ElemF);

DeEF=deformedMeshData(DeEF.dof,ElemF); % dof at t_n+1
MiEF=deformedMeshData(0.5*(Sol.me(:,:,1)+DeEF.dof-ElemF.dof)+ElemF.dof,ElemF); % mid mesh

if min(DeEF.detT)<0
    disp('detT<0 in feue2'); 
    pause;
    exit; 
end

% fluid matrices
[GM.fS,GM.fM,fC] = initMatrix2old(DeEF,ElemF); % t_n+1 dof
fM=Phy.frho*GM.fM; % fM contains rho

fR=(-Phy.frho*Sol.delt)*initMatrixdivuwv(MiEF,ElemF,Sol,DeEF.dof-ElemF.dof); % on mid mesh
[UxP,UyP] = initMatrix3(DeEF,ElemF); % t_n+1 dof
UxPUyP=(-Sol.delt)*[UxP;UyP]; % times -delt

zero1=sparse(ElemF.nDof,ElemF.nDof);
zero2=sparse(ElemF.pnDof,ElemF.pnDof);
coef=Phy.frho*Phy.fmu*Sol.delt;
fE=[[fM,zero1;zero1,fM]+coef*fC+fR,UxPUyP;UxPUyP',zero2];

% fluid forcing term
gg=Phy.g*[sin(Phy.theta),-cos(Phy.theta)];
Mg=fM*repmat(Sol.delt*gg,ElemF.nDof,1); % contains rho delt
bf=Sol.fMU+Mg;    
% if Phy.g~=0 && ~isempty(ElemF.openbdryElem) % strcmp(ElemF.meshString(1:4),'FPCf')==1
%     bfOut=uvTractionBC(Sol.delt,DeEF,ElemF,Phy); 
%     bf=bf(:)+bfOut(:);
% end

% solid forcing term

% this part has been added on 25 12 2012 as part of trail and testing to
% stabilize the structure part

DeESN = deformedMeshData(ElemS.dof+eta,ElemS); % used for linear elasticity part

if Sol.SVK==0  % Linear elasticity of the sol. is being used here...
    [sM,sCN,sDN] = initMatrix2(DeESN,ElemS);
        
    sM=Phy.srho*sM;
    % BDF-1
    sL=Sol.delt*Phy.smu*sCN+Sol.delt*Phy.slam*sDN; % contains delt factor
    zero1=sparse(ElemS.nDof,ElemS.nDof);
    sE=[sM,zero1;zero1,sM]+Sol.delt*sL; % delt^2
    
    GM.sE=sparse(GM.nd,GM.nd);
    GM.sE(GM.s2g,GM.s2g)=sE;
%     GM.sM=sM;
    if Sol.SVK==0, GM.sL=sL; end
    GM.sEItf=sE(GM.sEItfIx,GM.sEItfIx);
else
    GM.sM = Phy.srho*initMatrix2(DeES,ElemS);
end

% the added code end here..

bs=Sol.sMU+sM*(repmat(Sol.delt*gg,ElemS.nDof,1)); % sM contains rho

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
Ubd=Ubdry(DeEF.dof(ElemF.bdDirN,:),ElemF.meshString,Sol.t,Phy);
bb=bb-E(:,GM.bdDirNf)*Ubd(:);

uvpe=zeros(GM.nd,1);
uvpe(GM.uvpeFreeN)=E(GM.uvpeFreeN,GM.uvpeFreeN)\bb(GM.uvpeFreeN);

etat=reshape(uvpe(GM.s2g),[],2);

newEta = Sol.eta+Sol.delt*etat(:,:,1);

error = newEta(10,:)-eta(10,:)

if strcmp(meshString{1},'BenchMarkF')==1
    max(max(Sol.U(:,:,1)-Sol.U(:,:,2)))
    Sol.res(Sol.nit)=max(max(Sol.U(:,:,1)-Sol.U(:,:,2)));
    n1 = find(ElemF.dof(:,1)==1);
    Sol.U(n1,:)
    ElemF.dof(n1,:)
end

if strcmp(meshString{1},'BenchMarkF')==1
   n1 = find(ismember(ElemF.dof,[0.15 0.2],'rows')==1);
   n2 = find(ismember(ElemF.dof,[0.25 0.2],'rows')==1);
   Sol.DelP(Sol.nit,:) = Sol.p([n1;n2],:);
end

function [x,y]= polar2cart (mag, ang_in_deg)
 x = mag .* cos(ang_in_deg*pi/180);
 y = mag .* sin(ang_in_deg*pi/180);
  
 
function [r ad] = cart2polar(x)
 j = sqrt(-1);
 x = x(:,1)+j*x(:,2);
 r = abs(x);
 ar = angle(x);
 ad = ar*180/pi;

function angle = getAngle(dof,Elem)
d12 =sqrt(sum((dof(Elem.elem2dof(:,2),:)-dof(Elem.elem2dof(:,1),:)).^2,2));
d23 =sqrt(sum((dof(Elem.elem2dof(:,3),:)-dof(Elem.elem2dof(:,2),:)).^2,2));
d31 =sqrt(sum((dof(Elem.elem2dof(:,1),:)-dof(Elem.elem2dof(:,3),:)).^2,2));

angle = [acos((d12.^2+d31.^2-d23.^2)./(2*d12.*d31))...
    acos((d12.^2+d23.^2-d31.^2)./(2*d12.*d23))...
    acos((d23.^2+d31.^2-d12.^2)./(2*d23.*d31))]*180/pi;
