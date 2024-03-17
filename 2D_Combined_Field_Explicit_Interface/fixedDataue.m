function [ElemF,ElemS,Sol]=fixedDataue(meshString,Phy,Sol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  the function format has been modified it was actually
% [ElemF,ElemS,Sol,DeES]=fixedDataue(meshString,Phy,Sol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% u-eta formulation

global ALE GM %contains sM sL sE sEbdry fFixS

nGlobalRefine=Sol.nGlobalRefine;

% this step has been commented on 6th Jan 
% in the process of developing non linear fluid mesh

% ALE is for fluid phase only
Elem=fixedMeshData(1,meshString{1},nGlobalRefine);
ALE = MeshS(Elem,Elem.dof);

deg=Sol.deg; 

ElemF=fixedMeshData(deg,meshString{1},nGlobalRefine);

% ALE = SVKFluidMeshUpdate(ElemF);

% mixed formulation
ElemFp=fixedMeshData(deg-1,meshString{1},nGlobalRefine);
if max(abs(ElemF.bdElem-ElemFp.bdElem)>1e-10) | max(abs(ElemF.inElem-ElemFp.inElem))>1e-10
    % need the above equality in initMatrix3
    disp('error: p and u inElem and bdElem does not match');
    pause;
end
ElemF.pelem2dof=ElemFp.elem2dof;
ElemF.pinElem=ElemFp.inElem;
ElemF.pbdElem=ElemFp.bdElem;
ElemF.pitfNode=unique(ElemFp.elem2dof(setdiff(ElemFp.itfElem,ElemFp.citfElem),ElemFp.indLBN)); 
citfNode=unique(ElemFp.elem2dof(ElemFp.citfElem,ElemFp.cindLBN)); 
if length(ElemF.citfElem)==1
    ElemF.pitfNode = unique([ElemF.pitfNode;citfNode']);
else
    ElemF.pitfNode = unique([ElemF.pitfNode;citfNode]);
end
ElemF.pnDof=ElemFp.nDof;
ElemF.pnnt=ElemFp.nnt;
ElemF.pdeg=ElemFp.deg;

ElemS=fixedMeshData(deg,meshString{2},nGlobalRefine);

% % debug zone starts

ElemF = orderingIntfElem(ElemF);
ElemS = orderingIntfElem(ElemS);

% figure;
% showmesh(ElemS.dof,ElemS.elem2dof(:,1:3));
% % for i=1:length(ElemF.itfElem)
% %     p=ElemF.elem2dof(ElemF.itfElem(i),1:3);
% %     p=ElemF.dof(p,:);
% %     hold on;
% %     text(mean(p(:,1)),mean(p(:,2)),num2str(i));
% % end
% figure;
% showmesh(ElemF.dof,ElemF.elem2dof(:,1:3));
% % for i=1:length(ElemS.itfElem)
% %     p=ElemS.elem2dof(ElemS.itfElem(i),1:3);
% %     p=ElemS.dof(p,:);
% %     hold on;
% %     text(mean(p(:,1)),mean(p(:,2)),num2str(i));
% % end
if (strcmp(meshString{1},'FlexCurvedPlateF')==1)
    [r,d] = cart2polar(ElemF.dof(ElemF.itfNode,:)-repmat([0 4.983],size(ElemF.itfNode,1),1));
    a = round([r d]*10^3)*10^(-3);
    [~,itfOrder] = sortrows(a,[1 2]);
    ElemF.itfNode = ElemF.itfNode(itfOrder,:);
    
    [r,d] = cart2polar(ElemS.dof(ElemS.itfNode,:)-repmat([0 4.983],size(ElemS.itfNode,1),1));     
     a = round([r d]*10^3)*10^(-3);
    [~,itfOrder] = sortrows(a,[1 2]);
    ElemS.itfNode = ElemS.itfNode(itfOrder,:);
elseif (strcmp(meshString{1},'FS4f')==1)
    [~,itfOrder]=sortrows(ElemF.dof(ElemF.itfNode,:),1);
    ElemF.itfNode = ElemF.itfNode(itfOrder,:);
    [~,itfOrder]=sortrows(ElemS.dof(ElemS.itfNode,:),1);
    ElemS.itfNode = ElemS.itfNode(itfOrder,:);
else
    [~,itfOrder]=sortrows(ElemF.dof(ElemF.itfNode,:),[2,1]);
    ElemF.itfNode = ElemF.itfNode(itfOrder,:);
    [~,itfOrder]=sortrows(ElemS.dof(ElemS.itfNode,:),[2,1]);
    ElemS.itfNode = ElemS.itfNode(itfOrder,:);
end

if max(max(abs(ElemF.dof(ElemF.itfNode,:)-ElemS.dof(ElemS.itfNode,:))))>6e-4
    disp('mistake in fixedData');
    pause; return;
else
    ElemF.dof(ElemF.itfNode,:)=ElemS.dof(ElemS.itfNode,:);
end

% commented on 25 12 2012 as part of trail and testing to handle the
% instability in the structure part

DeES=deformedMeshData(ElemS.dof,ElemS); % used for linear elasticity part

% % debug zone starts

% p=ElemF.dof(ElemF.freeN,:);
% plot(p(:,1),p(:,2),'r*');
% hold on;
% showmesh(ElemF.dof,ElemF.elem2dof(:,1:3))
% p=ElemS.dof(ElemS.freeN,:);
% hold on;
% plot(p(:,1),p(:,2),'bo');
% hold on; 
% showmesh(ElemS.dof,ElemS.elem2dof(:,1:3))

% % debug zone ends testing
 
% from solid dof to global dof
% t1(i) is the global index of the ith dof in solid structure

if Sol.Conservative==0

    t1=(1:2*ElemS.nDof)';
    GM.s2g1 = (ElemF.nDof*2+ElemF.pnDof)+(1:length(t1));
    GM.nd1 = max(GM.s2g1);
    t2=setdiff(t1,[ElemS.itfNode;ElemS.nDof+ElemS.itfNode]);
    t1(ElemS.itfNode)=ElemF.itfNode;
    t1(ElemS.nDof+ElemS.itfNode)=ElemF.nDof+ElemF.itfNode;
    t1(t2)=(ElemF.nDof*2+ElemF.pnDof)+t2;

    GM.s2g=t1;
    GM.nd=max(t1); % number of uvpe dof
    t4=(ElemS.bdDirN); % exclude the end pts
    GM.bdDirN=[ElemF.bdDirN;ElemF.nDof+ElemF.bdDirN;t1(t4);t1(ElemS.nDof+t4)];
    GM.bdDirNf=[ElemF.bdDirN;ElemF.nDof+ElemF.bdDirN];

    t5=setdiff((1:ElemF.nDof)',ElemF.bdDirN);
    t6=setdiff((1:2*ElemS.nDof)',[ElemS.bdDirN;ElemS.nDof+ElemS.bdDirN;]);% ElemS.itfNode;ElemS.nDof+ElemS.itfNode
    GM.uvpeFreeN=[t5;ElemF.nDof+t5;2*ElemF.nDof+(1:ElemF.pnDof)';GM.s2g1(t6)'];

    t7=[2*ElemF.nDof+ElemF.pnDof+ElemS.itfNode;2*ElemF.nDof+ElemF.pnDof+ElemS.nDof+ElemS.itfNode];
    GM.sEItfIx=t7;
    GM.overlapIx=[ElemF.itfNode;ElemF.nDof+ElemF.itfNode];
    GM.nuv=2*ElemF.nDof;
    GM.nuvp=2*ElemF.nDof+ElemF.pnDof;

elseif Sol.Conservative==1
    
    t1=(1:2*ElemS.nDof)';
    t2=setdiff(t1,[ElemS.itfNode;ElemS.nDof+ElemS.itfNode]);
    t1(ElemS.itfNode)=ElemF.itfNode;
    t1(ElemS.nDof+ElemS.itfNode)=ElemF.nDof+ElemF.itfNode;
    t1(t2)=(ElemF.nDof*2+ElemF.pnDof)+(1:length(t2));

    GM.s2g=t1;
    GM.nd=max(t1); % number of uvpe dof
    t4=setdiff(ElemS.bdDirN,ElemS.itfNode); % exclude the end pts   
    GM.bdDirN=[ElemF.bdDirN;ElemF.nDof+ElemF.bdDirN;t1(t4);t1(ElemS.nDof+t4)];
    GM.bdDirNf=[ElemF.bdDirN;ElemF.nDof+ElemF.bdDirN];

    t5=setdiff((1:ElemF.nDof)',ElemF.bdDirN);
    t6=setdiff((1:2*ElemS.nDof)',[ElemS.bdDirN;ElemS.nDof+ElemS.bdDirN;]);% ElemS.itfNode;ElemS.nDof+ElemS.itfNode
    GM.uvpeFreeN=[t5;ElemF.nDof+t5;2*ElemF.nDof+(1:ElemF.pnDof)';t1(t6)];

    t7=[ElemS.itfNode;ElemS.nDof+ElemS.itfNode];
    GM.sEItfIx=t7;
    GM.overlapIx=[ElemF.itfNode;ElemF.nDof+ElemF.itfNode];
    GM.nuv=2*ElemF.nDof;
    GM.nuvp=2*ElemF.nDof+ElemF.pnDof;
    
end

% add this for filling the elastic bar
% only for FPC
if size(meshString,1)==2 && strcmp(meshString{2}(1:4),'FPCs')==1
    a=ElemF.elem2dof(ElemF.icElem,ElemF.indLBN);
    a=unique(a);
    nd=ElemF.dof(a,:);
    temp=atan2(nd(:,2)-0.2,nd(:,1)-0.25);
    [Y,Ind]=sort(temp,'descend');
%     snd=nd(Ind,:);
%     for i=1:length(snd), text(snd(i,1),snd(i,2),num2str(i)); hold on; end
%     plot(nd(Ind,1),nd(Ind,2),'r*');
%     fill(nd(Ind,1),nd(Ind,2),'m')
%     ElemF.barBDdof=a;
    ElemF.barBDInd{1}=a(Ind); 
end

if size(meshString,1)==2 && strcmp(meshString{2}(1:4),'FS4s')==1
    a=ElemF.elem2dof(ElemF.icElem,ElemF.indLBN);
    a=unique(a);
    nd=ElemF.dof(a,:);
    temp=atan2(nd(:,2)+0.25,nd(:,1)-1.5);
    [Y,Ind]=sort(temp,'descend');
%     snd=nd(Ind,:);
%     for i=1:length(snd), text(snd(i,1),snd(i,2),num2str(i)); hold on; end
%     plot(nd(Ind,1),nd(Ind,2),'r*');
%     fill(nd(Ind,1),nd(Ind,2),'m')
%     ElemF.barBDdof=a;
    ElemF.barBDInd{1}=a(Ind); 
end

if size(meshString,1)==2 && strcmp(meshString{2}(1:4),'HVs1')==1
    my=[-0.45,0.45]; % the first ellp for ElemF is on the top
    for i=1:2
        b=ElemF.itfElemInfo(i,:);
        a=ElemF.elem2dof(ElemF.itfElem(b(1):b(2)),ElemF.indLBN);
        a=unique(a);
        nd=ElemF.dof(a,:);
        temp=atan2(nd(:,2)+my(i),nd(:,1)-3);
        [Y,Ind]=sort(temp,'descend');
        %     snd=nd(Ind,:);
        %     for i=1:length(snd), text(snd(i,1),snd(i,2),num2str(i)); hold on; end
        %     plot(nd(Ind,1),nd(Ind,2),'r*');
        %     fill(nd(Ind,1),nd(Ind,2),'m')
        %     ElemF.barBDdof=a;
        ElemF.barBDInd{i}=a(Ind);
    end
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

function Elem = orderingIntfElem(Elem)
node = Elem.elem2dof(Elem.itfElem,4);
[node,order]=sortrows(Elem.dof(node,:),[1 2]);
Elem.itfElem = Elem.itfElem(order,:);
