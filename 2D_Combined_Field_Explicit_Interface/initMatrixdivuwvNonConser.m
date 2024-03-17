function R=initMatrixdivuwvNonConser(DeE,Elem,Sol,me )
% terms on t_n+1/2 and t_n-1/2mesh
% the input MiE is the middle mesh information. 

% The matrix is
% div (u o\times w-u^n) \cdot v + 1/2 div u^n w v 
% 
% where:
% w is mesh velocity and u^n is old velocity, v is test function, u is test function


mv=(me-1/2*Sol.me(:,:,1)-Sol.me(:,:,2)+1/2*Sol.me(:,:,3))./Sol.delt; % w^n = (a_i^n-1/2*a_i^(n-1)-a_i^(n-2)+1/2*a_i^(n-3))/(Del T)
[mv_qt,mvx_qt,mvy_qt]=uGraduAtQuadPts(mv,DeE,Elem);
[u_qt,ux_qt,uy_qt]=uGraduAtQuadPts(2*Sol.U(:,:,1)-Sol.U(:,:,2),DeE,Elem); % old U

NNT=Elem.nnt; % is the number of nodes on each triangle
NT=Elem.nElem;   N=Elem.nDof;
eQuadN=Elem.Quad.eQuad.N;

qtx=zeros(NT,NNT); qty=qtx;
for i=1:NNT
    qtx(:,i)=DeE.dof(Elem.elem2dof(:,i),1);
    qty(:,i)=DeE.dof(Elem.elem2dof(:,i),2);
end
locR=zeros(NT,eQuadN,NNT*NNT);
for i=1:eQuadN
    kappa=Elem.Quad.eQuad.kappa(i);
    psi=Elem.Quad.eQuad.psi(i,:);
    psi_x=Elem.Quad.eQuad.psi_x(i,:);
    psi_y=Elem.Quad.eQuad.psi_y(i,:);
    
    lAt=[qtx*psi_x',qtx*psi_y',qty*psi_x',qty*psi_y'];%(a11,a12;a21,a22)
    detAt=lAt(:,1).*lAt(:,4)-lAt(:,2).*lAt(:,3);
    
    % lAt^{-\top}=[a22,-a21,-a12,a11]./detlAt
    gradPhix=( lAt(:,4)*psi_x-lAt(:,3)*psi_y)./repmat(detAt,1,NNT);
    gradPhiy=(-lAt(:,2)*psi_x+lAt(:,1)*psi_y)./repmat(detAt,1,NNT);
    
    for j=1:NNT
        for k=1:NNT
            if Sol.Temam == 0
                locR(:,i,(j-1)*NNT+k)=locR(:,i,(j-1)*NNT+k)+...
                    (gradPhix(:,k).*(mv_qt(:,i,1)-u_qt(:,i,1))+...
                     gradPhiy(:,k).*(mv_qt(:,i,2)-u_qt(:,i,2))).*...
                    (psi(j)*kappa*(detAt));
            elseif Sol.Temam ==1
                locR(:,i,(j-1)*NNT+k)=locR(:,i,(j-1)*NNT+k)+...
                    (gradPhix(:,k).*(mv_qt(:,i,1)-u_qt(:,i,1))+...
                     gradPhiy(:,k).*(mv_qt(:,i,2)-u_qt(:,i,2))- ...
                     1/2*psi(k)*(ux_qt(:,i,1)+uy_qt(:,i,2))).*...
                    (psi(j)*kappa*(detAt));
            end
        end % end of for k=1:NNT
    end % end of for j=1:NNT
end % end of for i=1:eQuadN

locR=sum(locR,2);

jj=zeros(NT,NNT*NNT);
kk=zeros(NT,NNT*NNT);
R=sparse(N,N);
for j=1:NNT
    for k=1:NNT;
        jj(:,(j-1)*NNT+k)=Elem.elem2dof(:,j);
        kk(:,(j-1)*NNT+k)=Elem.elem2dof(:,k);
%         jj=Elem.elem2dof(:,j);
%         kk=Elem.elem2dof(:,k);        
%         R=R+sparse(jj,kk,locR(:,:,(j-1)*NNT+k),N,N);
    end
end

R=R+sparse(jj(:),kk(:),locR(:),N,N);
R=[R,sparse(N,N);sparse(N,N),R];