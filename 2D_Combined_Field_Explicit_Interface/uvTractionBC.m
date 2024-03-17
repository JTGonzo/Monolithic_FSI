function fpuv=uvTractionBC(delt,DeE,Elem,Phy)
% fpuv(:,1)=-rho*delt*<(gg(1)*x+gg(2)*y),\phi>_\Gam_2, 
% fpuv(:,2)=0
% not for general outward normal n yet. So far, n=(1,0)
% Gam_2 is only the open bdry

gg = delt*Phy.frho*Phy.g*[sin(Phy.theta),-cos(Phy.theta)];

bdQuad=Elem.Quad.bdQuad;
INDLBN=Elem.indLBN;

% % comment added on (June 24, 2011): 
% % In the old code: 
% % bdNeuElem also contain the interface
% % I guess this is the reason why the domain decomp code does not work 
% % when g is added.

% NbdNeuElem=length(Elem.bdNeuElem); 
% bdNeuNode=Elem.elem2dof(Elem.bdNeuElem,INDLBN);
NopenbdryElem=length(Elem.openbdryElem); 
openbdryNode=Elem.elem2dof(Elem.openbdryElem,INDLBN);

ss1=zeros(NopenbdryElem,length(INDLBN));
for i=1:NopenbdryElem 
    pt=DeE.dof(openbdryNode(i,:),:);
    if abs(pt(2,1)-Phy.outbdryx)<1e-8 % x cordi of mid bdry pt
        R_x=(bdQuad.psi_y(:,INDLBN)-bdQuad.psi_x(:,INDLBN))*pt;
        norm_R_x=sqrt(sum(R_x.^2,2)); %|dr/dt|
        ss1(i,:)=(norm_R_x.*bdQuad.delta.*(-gg(1)*pt(:,1)-gg(2)*pt(:,2)))'*bdQuad.psi(:,INDLBN);
    end
end
% assume normal vector is (1,0)
fpuv=[accumarray(openbdryNode(:),ss1(:),[Elem.nDof,1]),zeros(Elem.nDof,1)];