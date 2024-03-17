function [rr]=ldc(Elem,Phy,mix,u,p,DeEF,meshString)
% lift and drag coefficients for flow past a cylinder with flag

NNT=Elem.nnt; % is the number of nodes on each triangle
if (strcmp(meshString{1},'FlexFlatPlateF')==1) || (strcmp(meshString{1},'BenchMarkF')==1)
    elem = Elem.itfElem;
elseif (strcmp(meshString{1}(1:4),'FPCf')==1)
    elem=Elem.icElem;
    elem=[elem;Elem.itfElem];
end
NT = length(elem);
bdQuadN=Elem.Quad.bdQuad.N;
if mix==1
    NNTp=Elem.pnnt;
    xx=Elem.Quad.bdQuad.x;
    ppsi=basefun(Elem.pdeg,xx,1-xx);
else
    NNTp=NNT;
    ppsi=Elem.Quad.bdQuad.psi;
end

rr=zeros(1,2);
rr1 = zeros(1,2);
uhNode=zeros(NT,NNT,2);
for i=1:NNT, for j=1:2, uhNode(:,i,j)=u(Elem.elem2dof(elem,i),j); end, end
phNode=zeros(NT,NNTp);
if mix==1
    for i=1:NNTp, phNode(:,i)=p(Elem.pelem2dof(elem,i)); end
else
    for i=1:NNTp, phNode(:,i)=p(Elem.elem2dof(elem,i)); end
end
% figure; str=['r*';'bo';'yd']
% for i=1:3
%     a=Elem.dof(Elem.elem2dof(Elem.icElem,i),:);
%     hold on; plot(a(:,1),a(:,2),str(i,:));
% end

qtx=zeros(size(elem,1),NNT); qty=qtx;
for i=1:NNT
    qtx(:,i)=DeEF.dof(Elem.elem2dof(elem,i),1);
    qty(:,i)=DeEF.dof(Elem.elem2dof(elem,i),2);
end

for i=1:bdQuadN
    delta=Elem.Quad.bdQuad.delta(i);
    psi_x=Elem.Quad.bdQuad.psi_x(i,:);
    psi_y=Elem.Quad.bdQuad.psi_y(i,:);
    lAt=[qtx*psi_x',qtx*psi_y',qty*psi_x',qty*psi_y'];%(a11,a12;a21,a22)
    detlAt=lAt(:,1).*lAt(:,4)-lAt(:,2).*lAt(:,3);
    % lAt^{-\top}=[a22,-a21,-a12,a11]./detlAt
    gradPhix=( lAt(:,4)*psi_x-lAt(:,3)*psi_y)./repmat(detlAt,1,NNT);
    gradPhiy=(-lAt(:,2)*psi_x+lAt(:,1)*psi_y)./repmat(detlAt,1,NNT);
    F=[sum(uhNode(:,:,1).*gradPhix,2),...
        sum(uhNode(:,:,1).*gradPhiy,2),...
        sum(uhNode(:,:,2).*gradPhix,2),...
        sum(uhNode(:,:,2).*gradPhiy,2)];
    
    % r=sum_i P_i phi_i, dr/dt=dr/dx-dr/dy since t=x and y=1-x  
    tau=[lAt(:,1)-lAt(:,2),lAt(:,3)-lAt(:,4)];
    jacobian=sqrt(sum(tau.^2,2));
    nn=[-tau(:,2),tau(:,1)]./[jacobian,jacobian]; % outward normal   
    % p=[qtx*psi',qty*psi'];
    % quiver(p(:,1),p(:,2),nn(:,1),nn(:,2));

    % \int_\Gam nu (\grad u + \grad u^\top) - p I \nn 
    F=Phy.frho*Phy.fmu*symm2(F);
    
    PP=phNode*ppsi(i,:)';
    F(:,1)=F(:,1)-PP;
    F(:,4)=F(:,4)-PP;
    f(:,[1 4]) = [-PP -PP];
    f(:,[2 3]) = [F(:,2) F(:,3)];
    
    RR = Ab(F,nn);
    
    rr=rr+sum(repmat(delta*jacobian,1,2).*RR,1);
end % end of for i=1:bdQuadN

function C=Ab(A,b)
% C = A b. A is stored as [a11,a12,a21,a22]
C=[A(:,1).*b(:,1)+A(:,2).*b(:,2),-(A(:,3).*b(:,1)+A(:,4).*b(:,2))];

function C=symm2(A)
% C = (A^\top + A)
temp=A(:,2)+A(:,3);
C=[2*A(:,1),temp,temp,2*A(:,4)];