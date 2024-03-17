function ALE=SVKFluidMeshUpdateNonLinearElastic(ALE,dof,u)
% L is linar component in the Div. of stress
% N is the non linear component in the Div. of stress

ALE.Elem.dof = dof;
Elem = ALE.Elem;

% values of rubber
E = 0.01;
nu = 0.45; 
mu = E/(2+2*nu);
lam = E*nu/((1+nu)*(1-2*nu));

eQuadN=Elem.Quad.eQuad.N;
NNT=Elem.nnt; % is the number of nodes on each triangle
N = Elem.nDof;

L1_x = sparse(N,N);
L2_x = sparse(N,N);
L1_y = sparse(N,N);
L2_y = sparse(N,N); 

N1_x = sparse(N,N);
N2_x = sparse(N,N);
N1_y = sparse(N,N);
N2_y = sparse(N,N);   

qtx=zeros(Elem.nElem,NNT); qty=qtx;
for i=1:NNT
    qtx(:,i)=Elem.dof(Elem.elem2dof(:,i),1);
    qty(:,i)=Elem.dof(Elem.elem2dof(:,i),2);
end

ux=zeros(Elem.nElem,NNT); uy=qtx;
for i=1:NNT
    ux(:,i)=u(Elem.elem2dof(:,i),1);
    uy(:,i)=u(Elem.elem2dof(:,i),2);
end

for i=1:eQuadN
    kappa=Elem.Quad.eQuad.kappa(i);
    % psi=Elem.Quad.eQuad.psi(i,:);
    psi_x=Elem.Quad.eQuad.psi_x(i,:);
    psi_y=Elem.Quad.eQuad.psi_y(i,:);
    lAt=[qtx*psi_x',qtx*psi_y',qty*psi_x',qty*psi_y'];%(a11,a12;a21,a22)
    U = [ux*psi_x',ux*psi_y',uy*psi_x',uy*psi_y']; % par. Deri Ux with x, par. Deri Ux with y, par. Deri Uy with x, par. Deri Uy with y
    
    detlAt=lAt(:,1).*lAt(:,4)-lAt(:,2).*lAt(:,3);
    % lAt^{-\top}=[a22,-a21,-a12,a11]./detlAt
    gradPhix=( lAt(:,4)*psi_x-lAt(:,3)*psi_y)./repmat(detlAt,1,NNT);
    gradPhiy=(-lAt(:,2)*psi_x+lAt(:,1)*psi_y)./repmat(detlAt,1,NNT);
    tau=1;1+(max(detlAt)-min(detlAt))./(detlAt).^1.25;

    for j=1:NNT
        for k=1:NNT
            l1_x = (mu.*(2*gradPhix(:,j).*gradPhix(:,k)+gradPhiy(:,j).*gradPhiy(:,k))+lam.*gradPhix(:,j).*gradPhix(:,k))*kappa.*tau.*detlAt;
            l1_y = (lam.*gradPhiy(:,j).*gradPhix(:,k)+mu.*gradPhix(:,j).*gradPhiy(:,k))*kappa.*tau.*detlAt;
            l2_x = (mu.*gradPhiy(:,j).*gradPhix(:,k)+lam.*gradPhix(:,j).*gradPhiy(:,k))*kappa.*tau.*detlAt;
            l2_y = (mu.*(2*gradPhiy(:,j).*gradPhiy(:,k)+gradPhix(:,j).*gradPhix(:,k))+lam.*gradPhiy(:,j).*gradPhiy(:,k))*kappa.*tau.*detlAt;
                      
            n1_x = ((gradPhix(:,k).*U(:,1)+gradPhiy(:,k).*U(:,2)).*gradPhix(:,j).*lam/2+...
                (gradPhix(:,k).*U(:,1)).*gradPhix(:,j).*mu+(gradPhix(:,k).*U(:,2)).*gradPhiy(:,j).*mu)*kappa.*tau.*detlAt;
            n1_y = ((gradPhix(:,k).*U(:,3)+gradPhiy(:,k).*U(:,4)).*gradPhix(:,j).*lam/2+...
                (gradPhix(:,k).*U(:,3)).*gradPhix(:,j).*mu+(gradPhix(:,k).*U(:,4)).*gradPhiy(:,j).*mu)*kappa.*tau.*detlAt;
            n2_x = ((gradPhix(:,k).*U(:,1)+gradPhiy(:,k).*U(:,2)).*gradPhiy(:,j).*lam/2+...
                (gradPhiy(:,k).*U(:,1)).*gradPhix(:,j).*mu+(gradPhiy(:,k).*U(:,2)).*gradPhiy(:,j).*mu)*kappa.*tau.*detlAt;
            n2_y = ((gradPhix(:,k).*U(:,3)+gradPhiy(:,k).*U(:,4)).*gradPhiy(:,j).*lam/2+...
                (gradPhix(:,k).*U(:,3)).*gradPhix(:,j).*mu+(gradPhiy(:,k).*U(:,4)).*gradPhiy(:,j).*mu)*kappa.*tau.*detlAt;
         
            c = Elem.elem2dof(:,k);
            r = Elem.elem2dof(:,j);
            
            L1_x = L1_x + sparse(r,c,l1_x,N,N);
            L2_x = L2_x + sparse(r,c,l2_x,N,N);
            L1_y = L1_y + sparse(r,c,l1_y,N,N);
            L2_y = L2_y + sparse(r,c,l2_y,N,N); 
            
            N1_x = N1_x + sparse(r,c,n1_x,N,N);
            N2_x = N2_x + sparse(r,c,n2_x,N,N);
            N1_y = N1_y + sparse(r,c,n1_y,N,N);
            N2_y = N2_y + sparse(r,c,n2_y,N,N);            
        end % end of for k=1:NNT
    end % end of for j=1:NNT
end % end of for i=1:eQuadN
ALE.L = [L1_x L1_y;L2_x L2_y];
ALE.N = [N1_x N1_y;N2_x N2_y];

ALE.nVertex=Elem.nDof;
ALE.fFreeN = setdiff((1:Elem.nDof)',[Elem.bdNeuN;Elem.bdDirN]);
ALE.fFreeNodes = setdiff((1:Elem.nVertex)',[Elem.bdNeuN;Elem.bdDirN]);
