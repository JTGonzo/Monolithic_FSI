function [A,B]=initMatrix3(DeE,Elem)
% A is <u_x, p> and B is <u_y, p>
% U deg and P deg-1
NNTu=Elem.nnt; % is the number of u nodes on each triangle
NNTp=Elem.pnnt;
switch Elem.deg
    case 2
        strBaseName=['base',num2str(21),'.mat'];
    case 3
        strBaseName=['base',num2str(32),'.mat'];
    case 4
        strBaseName=['base',num2str(43),'.mat'];
    case 5
        strBaseName=['base',num2str(54),'.mat'];
    otherwise
        disp(['error: in initMatrix3, not implemented yet']);
        pause;
end
load(strBaseName);
NT = Elem.nElem;   Nu = Elem.nDof;   Np=Elem.pnDof;
A = sparse(Nu,Np);  B = sparse(Nu,Np);

% % debug
% p1p2 = DeE.dof(Elem.elem2dof(:,2),:)...
%     -DeE.dof(Elem.elem2dof(:,1),:); % edge 12
% p3p1 = DeE.dof(Elem.elem2dof(:,1),:)...
%     -DeE.dof(Elem.elem2dof(:,3),:); % edge 31
% for i = 1:NNTu
%     for j = 1:NNTp
%         % triangle
%         Aij  =-UxP(i,j)*p3p1(:,2)-UyP(i,j)*p1p2(:,2);
%         Bij  = UxP(i,j)*p3p1(:,1)+UyP(i,j)*p1p2(:,1);
%         ti=Elem.elem2dof(:,i);
%         tj=Elem.elem2dof(:,j); % j=1:3 is the P1 nodes
%         A = A + sparse(ti,tj,Aij,Nu,Np);
%         B = B + sparse(ti,tj,Bij,Nu,Np);
%     end
% end
% return
% end of debug

for i = 1:NNTu
    for j = 1:NNTp
        % triangle
        Aij  =-UxP(i,j)*DeE.p3p1(:,2)-UyP(i,j)*DeE.p1p2(:,2);
        Bij  = UxP(i,j)*DeE.p3p1(:,1)+UyP(i,j)*DeE.p1p2(:,1);
        ti=Elem.elem2dof(Elem.inElem,i);
        tj=Elem.pelem2dof(Elem.pinElem,j); % p index
        A = A + sparse(ti,tj,Aij,Nu,Np);
        B = B + sparse(ti,tj,Bij,Nu,Np);
    end
end
%--------------------------------------------------------------------------
% handling isoparametric elements
%--------------------------------------------------------------------------
NIso=length(Elem.bdElem); ind=0;
NNTup=NNTu*NNTp;
ii=zeros(NNTup*NIso,1); jj=ii; ssa=ii; ssb=ii;
kappa=Elem.Quad.eQuad.kappa;
xy=Elem.Quad.eQuad.xy;
ppsi=basefun(Elem.pdeg,xy(1,:),xy(2,:));% P1
upsi_x=Elem.Quad.eQuad.psi_x;
upsi_y=Elem.Quad.eQuad.psi_y;
eQuadN=Elem.Quad.eQuad.N;
for i=1:NIso
    dofId=Elem.elem2dof(Elem.bdElem(i),:); % index of the dof on this element
    pdofId=Elem.pelem2dof(Elem.pbdElem(i),:);
    Aiso=zeros(NNTu,NNTp);
    Biso=zeros(NNTu,NNTp);
    for k=1:eQuadN
        G=DeE.bdInvDF(:,2*k-1:2*k,i)*([upsi_x(k,:);upsi_y(k,:)]);
        Aiso=Aiso+(kappa(k)*DeE.bdAbsDF(k,i))*(G(1,:)'*ppsi(k,:)); % No *nu
        Biso=Biso+(kappa(k)*DeE.bdAbsDF(k,i))*(G(2,:)'*ppsi(k,:));
    end
    ssa(ind+1:ind+NNTup)=reshape(Aiso',NNTup,1);
    ssb(ind+1:ind+NNTup)=reshape(Biso',NNTup,1);
    
    temp=repmat(dofId,NNTp,1);
    ii(ind+1:ind+NNTup)=reshape(temp,NNTup,1);
    temp=repmat(pdofId,NNTu,1); % p index
    jj(ind+1:ind+NNTup)=reshape(temp',NNTup,1);
    ind=ind+NNTup;
%         for ti=1:NNTu
%             for tj=1:NNTp
%                 ind=ind+1;
%                 ii(ind)=elem2dof(eId,ti);jj(ind)=elem2dof(eId,tj);
%                 ss(ind)=Aiso(ti,tj);
%             end
%         end
end
A = round((A + sparse(ii,jj,ssa,Nu,Np))*10^13)*10^(-13);
B = round((B + sparse(ii,jj,ssb,Nu,Np))*10^13)*10^(-13);