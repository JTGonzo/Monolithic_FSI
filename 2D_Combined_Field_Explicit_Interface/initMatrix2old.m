function [A,B,C,D]=initMatrix2old(DeE,Elem)

% A is stiff matrix and B is mass matrix
% C is <\grad \phi + \grad \phi^\top), \grad \psi> and 
% D is <\div\phiI , \grad \psi> = <\div\phi, \div\psi>

NNT=Elem.nnt; % is the number of nodes on each triangle
strBaseName=['base',num2str(Elem.deg),'.mat'];
load(strBaseName);
NT = Elem.nElem;   N = Elem.nDof;
A = sparse(N,N);  B = sparse(N,N);
% C1 is <phi_x,psi_x>, C2 is <phi_y, psi_x>, C3 is <phi_y, psi_y> 
C1 = sparse(N,N); C2 = sparse(N,N); C3 = sparse(N,N);
a = [DeE.p3p1(:,1).^2+DeE.p3p1(:,2).^2,...
    DeE.p1p2(:,2).*DeE.p3p1(:,2)+DeE.p1p2(:,1).*DeE.p3p1(:,1),...
    DeE.p1p2(:,1).^2+DeE.p1p2(:,2).^2];
invDT=1./abs(DeE.detT);
for i = 1:NNT
    for j = 1:NNT
        % triangle
        Aij  = UxUx(i,j)*a(:,1)...
               +(UxUy(i,j)+UxUy(j,i))*a(:,2)...
               +UyUy(i,j)*a(:,3);
        C1ij = UxUx(i,j)*(DeE.p3p1(:,2).^2)...
               +(UxUy(i,j)+UxUy(j,i))*(DeE.p3p1(:,2).*DeE.p1p2(:,2))...
               +UyUy(i,j)*(DeE.p1p2(:,2).^2);
        C2ij = -UxUx(i,j)*(DeE.p3p1(:,1).*DeE.p3p1(:,2))...
               -UxUy(i,j)*(DeE.p3p1(:,1).*DeE.p1p2(:,2))...
               -UxUy(j,i)*(DeE.p3p1(:,2).*DeE.p1p2(:,1))...
               -UyUy(i,j)*(DeE.p1p2(:,1).*DeE.p1p2(:,2));
        C3ij = UxUx(i,j)*(DeE.p3p1(:,1).^2)...
               +(UxUy(i,j)+UxUy(j,i))*(DeE.p3p1(:,1).*DeE.p1p2(:,1))...
               +UyUy(i,j)*(DeE.p1p2(:,1).^2);    
        Aij=Aij.*invDT;
        C1ij=C1ij.*invDT;
        C2ij=C2ij.*invDT;
        C3ij=C3ij.*invDT;
        ti=Elem.elem2dof(Elem.inElem,i);
        tj=Elem.elem2dof(Elem.inElem,j);
        A = A + sparse(ti,tj,Aij,N,N);
        Bij = UU(i,j)*abs(DeE.detT);
        B = B + sparse(ti,tj,Bij,N,N);
        C1 = C1 + sparse(ti,tj,C1ij,N,N);
        C2 = C2 + sparse(ti,tj,C2ij,N,N);
        C3 = C3 + sparse(ti,tj,C3ij,N,N);        
    end
end
%--------------------------------------------------------------------------
% handling isoparametric elements
%--------------------------------------------------------------------------
NNT2=NNT^2;
NIso=length(Elem.bdElem); ind=0; 
ii=zeros(NNT2*NIso,1); jj=ii; ssa=ii; ssb=ii;
kappa=Elem.Quad.eQuad.kappa;
psi=Elem.Quad.eQuad.psi;
psi_x=Elem.Quad.eQuad.psi_x;
psi_y=Elem.Quad.eQuad.psi_y;
eQuadN=Elem.Quad.eQuad.N;
for i=1:NIso
    eId=Elem.bdElem(i); % index of the element
    dofId=Elem.elem2dof(eId,:); % index of the dof on this element
    Aiso=zeros(NNT,NNT);
    Biso=zeros(NNT,NNT);
    C1iso=zeros(NNT,NNT);
    C2iso=zeros(NNT,NNT);
    C3iso=zeros(NNT,NNT);
    for k=1:eQuadN
        G=DeE.bdInvDF(:,2*k-1:2*k,i)*([psi_x(k,:);psi_y(k,:)]);
        Aiso=Aiso+(kappa(k)*DeE.bdAbsDF(k,i))*(G'*G); % No *nu
        Biso=Biso+(kappa(k)*DeE.bdAbsDF(k,i))*(psi(k,:)'*psi(k,:));
        C1iso=C1iso+(kappa(k)*DeE.bdAbsDF(k,i))*(G(1,:)'*G(1,:));
        C2iso=C2iso+(kappa(k)*DeE.bdAbsDF(k,i))*(G(2,:)'*G(1,:));
        C3iso=C3iso+(kappa(k)*DeE.bdAbsDF(k,i))*(G(2,:)'*G(2,:));
    end
    ssa(ind+1:ind+NNT2)=reshape(Aiso',NNT2,1);
    ssb(ind+1:ind+NNT2)=reshape(Biso',NNT2,1);
    ssc1(ind+1:ind+NNT2)=reshape(C1iso',NNT2,1);
    ssc2(ind+1:ind+NNT2)=reshape(C2iso',NNT2,1);
    ssc3(ind+1:ind+NNT2)=reshape(C3iso',NNT2,1);
    
    temp=repmat(dofId,NNT,1);
    ii(ind+1:ind+NNT2)=reshape(temp,NNT2,1);
    jj(ind+1:ind+NNT2)=reshape(temp',NNT2,1);
    ind=ind+NNT2;
%         for ti=1:NNT
%             for tj=1:NNT
%                 ind=ind+1;
%                 ii(ind)=elem2dof(eId,ti);jj(ind)=elem2dof(eId,tj);
%                 ss(ind)=Aiso(ti,tj);
%             end
%         end
end
A = A + sparse(ii,jj,ssa,N,N);
B = B + sparse(ii,jj,ssb,N,N);
C1 = C1 + sparse(ii,jj,ssc1,N,N);
C2 = C2 + sparse(ii,jj,ssc2,N,N);
C3 = C3 + sparse(ii,jj,ssc3,N,N);
C=[A+C1,C2;C2',A+C3];
if nargout==4
    D=[C1,C2';C2,C3];
end