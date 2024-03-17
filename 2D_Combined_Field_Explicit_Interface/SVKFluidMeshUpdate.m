function ALE=SVKFluidMeshUpdate(Elem)
% L is linar component in the Div. of stress
% N is the non linear component in the Div. of stress

NNT=Elem.nnt; % is the number of nodes on each triangle
N=Elem.nDof;
eQuadN=Elem.Quad.eQuad.N;

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
    tau=1+(max(detlAt)-min(detlAt))./detlAt;

    for j=1:NNT
        for k=1:NNT
            l1_x = (3*gradPhix(:,j).*gradPhix(:,k))*kappa.*tau.*detlAt;
            l2_x = (gradPhiy(:,j).*gradPhix(:,k)+gradPhix(:,j).*gradPhiy(:,k))*kappa.*tau.*detlAt;
            l1_y = (gradPhiy(:,j).*gradPhix(:,k)+gradPhix(:,j).*gradPhiy(:,k))*kappa.*tau.*detlAt;
            l2_y = (3*gradPhiy(:,j).*gradPhiy(:,k))*kappa.*tau.*detlAt;
            n1_x = ((gradPhix(:,j).*gradPhix(:,j)+gradPhiy(:,j).*gradPhiy(:,j)).*gradPhix(:,k)+...
                (gradPhix(:,j).*gradPhix(:,j)).*gradPhix(:,k)+(gradPhix(:,j).*gradPhiy(:,j)).*gradPhiy(:,k))*kappa.*tau.*detlAt;
            n2_x = ((gradPhix(:,j).*gradPhix(:,j)+gradPhiy(:,j).*gradPhiy(:,j)).*gradPhix(:,k)+...
                (gradPhix(:,j).*gradPhix(:,j)).*gradPhix(:,k)+(gradPhix(:,j).*gradPhiy(:,j)).*gradPhiy(:,k))*kappa.*tau.*detlAt;
            n1_y = ((gradPhix(:,j).*gradPhix(:,j)+gradPhiy(:,j).*gradPhiy(:,j)).*gradPhiy(:,k)+...
                (gradPhix(:,j).*gradPhiy(:,j)).*gradPhix(:,k)+(gradPhiy(:,j).*gradPhiy(:,j)).*gradPhiy(:,k))*kappa.*tau.*detlAt;
            n2_y = ((gradPhix(:,j).*gradPhix(:,j)+gradPhiy(:,j).*gradPhiy(:,j)).*gradPhiy(:,k)+...
                (gradPhix(:,j).*gradPhiy(:,j)).*gradPhix(:,k)+(gradPhiy(:,j).*gradPhiy(:,j)).*gradPhiy(:,k))*kappa.*tau.*detlAt;
            r = Elem.elem2dof(:,j);
            c = Elem.elem2dof(:,k);
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
ALE.L = [L1_x L2_x;L1_y L2_y];
ALE.N = [N1_x N2_x;N1_y N2_y];
% the freeN of E is not the same as that of Elem
ALE.nVertex=Elem.nVertex;
ALE.fFreeVertex = setdiff((1:ALE.nVertex)',[Elem.bdNeuN;Elem.bdDirN]);