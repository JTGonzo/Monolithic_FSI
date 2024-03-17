function [u_qt,ux_qt,uy_qt] = uGraduAtQuadPts(uh,DeE,Elem)
nnu=size(uh,2); % how many u's. either u, (u,v) or (u,v,p) or (p,q) etc.
NNT=Elem.nnt;   % is the number of nodes on each triangle
NT=Elem.nElem;
uhNode=zeros([length(Elem.inElem),NNT,nnu]); % only used on Triangle elements
eQuad=Elem.Quad.eQuad;
u_qt=zeros(NT,eQuad.N,nnu); 
for k=1:nnu
    for i=1:NNT
        uhNode(:,i,k)=uh(Elem.elem2dof(Elem.inElem,i),k);
    end
    u_qt(Elem.inElem,:,k)=uhNode(:,:,k)*eQuad.psi';
end
if nargout==3
    ux_qt=zeros(NT,eQuad.N,nnu); uy_qt=ux_qt;
    for k=1:nnu
        for i=1:eQuad.N
            temp=uhNode(:,:,k)*[eQuad.psi_x(i,:)',eQuad.psi_y(i,:)'];
            ux_qt(Elem.inElem,i,k)=-(DeE.p3p1(:,2).*temp(:,1)+...
                DeE.p1p2(:,2).*temp(:,2))./DeE.detT;
            uy_qt(Elem.inElem,i,k)= (DeE.p3p1(:,1).*temp(:,1)+...
                DeE.p1p2(:,1).*temp(:,2))./DeE.detT;
        end
    end
end
%--------------------------------------------------------------------------
% handling isoparametric elements
%--------------------------------------------------------------------------
uhNode=zeros(Elem.nBDElem,NNT,nnu);
for k=1:nnu
    for i=1:NNT
        uhNode(:,i,k)=uh(Elem.elem2dof(Elem.bdElem,i),k);
    end
    u_qt(Elem.bdElem,:,k)=uhNode(:,:,k)*eQuad.psi';
end

if nargout==3
    gradu_loc=zeros(2,eQuad.N);    
    for i=1:Elem.nBDElem
        eId=Elem.bdElem(i); % index of the element
        % dofId=Elem.elem2dof(eId,:); % index of the dof on this element
        for n=1:nnu
            for k=1:eQuad.N
                gradu_loc(:,k) = DeE.bdInvDF(:,2*k-1:2*k,i)*...
                    [eQuad.psi_x(k,:);eQuad.psi_y(k,:)]*uhNode(i,:,n)';
            end
            ux_qt(eId,:,n) = gradu_loc(1,:);
            uy_qt(eId,:,n) = gradu_loc(2,:);
        end
    end
end