function S=initMatrixSVK2(Elem,Phy,u)
% S is lam <trE^n grad phi_i, grad phi_j> 
% + 2*mu < grad  phi_i, grad phi_j>_{E^n}     
% u is the displacement vector

lam=Phy.slam;
mu=Phy.smu;

NNT=Elem.nnt; % is the number of nodes on each triangle
NT=Elem.nElem;   N=Elem.nDof;
eQuadN=Elem.Quad.eQuad.N;

S=sparse(2*N,2*N);
uhNode=zeros(NT,NNT,2);
for i=1:NNT, for j=1:2, uhNode(:,i,j)=u(Elem.elem2dof(:,i),j); end, end

% qtx=[]; qty=qtx;
% for i=1:NNT
%     qtx=[qtx,Elem.dof(Elem.elem2dof(:,i),1)];
%     qty=[qty,Elem.dof(Elem.elem2dof(:,i),2)];
% end
 
qtx=zeros(Elem.nElem,NNT); qty=qtx;
for i=1:NNT
    qtx(:,i)=Elem.dof(Elem.elem2dof(:,i),1);
    qty(:,i)=Elem.dof(Elem.elem2dof(:,i),2);
end


for i=1:eQuadN
    kappa=Elem.Quad.eQuad.kappa(i);
    % psi=Elem.Quad.eQuad.psi(i,:);
    psi_x=Elem.Quad.eQuad.psi_x(i,:);
    psi_y=Elem.Quad.eQuad.psi_y(i,:);
    lAt=[qtx*psi_x',qtx*psi_y',qty*psi_x',qty*psi_y'];%(a11,a12;a21,a22)
    detlAt=lAt(:,1).*lAt(:,4)-lAt(:,2).*lAt(:,3);
    % lAt^{-\top}=[a22,-a21,-a12,a11]./detlAt
    gradPhix=( lAt(:,4)*psi_x-lAt(:,3)*psi_y)./repmat(detlAt,1,NNT);
    gradPhiy=(-lAt(:,2)*psi_x+lAt(:,1)*psi_y)./repmat(detlAt,1,NNT);

    F=[sum(uhNode(:,:,1).*gradPhix,2)+1,...
       sum(uhNode(:,:,1).*gradPhiy,2),...
       sum(uhNode(:,:,2).*gradPhix,2),...
       sum(uhNode(:,:,2).*gradPhiy,2)+1];
    E=0.5*(AtB(F,F)-repmat([1,0,0,1],NT,1));

    for j=1:NNT
        for k=1:NNT
            % jj=Elem.elem2dof(:,j); kk=Elem.elem2dof(:,k);

            for m=1:4 % (1,0):(1,0),(1,0):(0,1),(0,1):(1,0),(0,1):(0,1)
                jj=Elem.elem2dof(:,j);
                kk=Elem.elem2dof(:,k);

                if (m==1 || m==2)
                    G=[gradPhix(:,j),gradPhiy(:,j),zeros(NT,2)];
                else
                    G=[zeros(NT,2),gradPhix(:,j),gradPhiy(:,j)];
                    jj=jj+N;
                end

                if (m==1 || m==3)
                    H=[gradPhix(:,k),gradPhiy(:,k),zeros(NT,2)];
                else
                    H=[zeros(NT,2),gradPhix(:,k),gradPhiy(:,k)];
                    kk=kk+N;
                end

                localS=lam*(trA(E).*sum(G.*H,2))...
                    +(2*mu)*(sum(E.*symm(G,H),2));

                S=S+sparse(jj,kk,(kappa*detlAt).*localS,2*N,2*N);
            end 
        end % end of for k=1:NNT
    end % end of for j=1:NNT
end % end of for i=1:eQuadN


% function C=AB(A,B)
% % C = A B. But A and B are stored as [a11,a12,a21,a22]
% C=[A(:,1).*B(:,1)+A(:,2).*B(:,3), A(:,1).*B(:,2)+A(:,2).*B(:,4),...
%    A(:,3).*B(:,1)+A(:,4).*B(:,3), A(:,3).*B(:,2)+A(:,4).*B(:,4)];

function C=AtB(A,B)
% C = A^\top B. But A and B are stored as [a11,a12,a21,a22]
C=[A(:,1).*B(:,1)+A(:,3).*B(:,3), A(:,1).*B(:,2)+A(:,3).*B(:,4),...
   A(:,2).*B(:,1)+A(:,4).*B(:,3), A(:,2).*B(:,2)+A(:,4).*B(:,4)];

function t=trA(A)
% A is stored as [a11,a12,a21,a22]    
t=A(:,1)+A(:,4);

function C=symm(A,B)
if nargin==2
    % C = (A^\top B + B^\top A)/2 and A and B are stored as [a11,a12,a21,a22]
    temp=0.5*(A(:,1).*B(:,2)+A(:,3).*B(:,4)+A(:,2).*B(:,1)+A(:,4).*B(:,3));
    C=[A(:,1).*B(:,1)+A(:,3).*B(:,3),temp,...
        temp,A(:,2).*B(:,2)+A(:,4).*B(:,4)];
end
if nargin==1
    % C = (A^\top + A)/2
    temp=0.5*(A(:,2)+A(:,3));
    C=[A(:,1),temp,temp,A(:,4)];
end