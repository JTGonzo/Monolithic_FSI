function R=initMatrixdivuwv2(MiE1,MiE2,DeE,Elem,Sol,me )
% terms on t_n+1/2 and t_n-1/2mesh
% the input MiE is the middle mesh information. 

% The matrix is
% div (u o\times w-u^n) \cdot v + 1/2 div u^n w v 
% 
% where:
% w is mesh velocity and u^n is old velocity, v is test function, u is test function


mv1=(me-Sol.me(:,:,1))./Sol.delt;
mv2 = (Sol.me(:,:,1)-Sol.me(:,:,2))./Sol.delt;
[mv1_qt,mvx1_qt,mvy1_qt]=uGraduAtQuadPts(mv1,MiE1,Elem);
[mv2_qt,mvx2_qt,mvy2_qt]=uGraduAtQuadPts(mv2,MiE2,Elem);
[u_qt,ux_qt,uy_qt]=uGraduAtQuadPts(2*Sol.U(:,:,1)-Sol.U(:,:,2),DeE,Elem); % old U

NNT=Elem.nnt; % is the number of nodes on each triangle
NT=Elem.nElem;   N=Elem.nDof;
eQuadN=Elem.Quad.eQuad.N;

qtx1=zeros(NT,NNT); qty1=qtx1;qtx2=zeros(NT,NNT); qty2=qtx2;qtx=zeros(NT,NNT); qty=qtx;
for i=1:NNT
    qtx(:,i)=DeE.dof(Elem.elem2dof(:,i),1);
    qty(:,i)=DeE.dof(Elem.elem2dof(:,i),2);
    
    qtx1(:,i)=MiE1.dof(Elem.elem2dof(:,i),1);
    qty1(:,i)=MiE1.dof(Elem.elem2dof(:,i),2);
    
    qtx2(:,i)=MiE2.dof(Elem.elem2dof(:,i),1);
    qty2(:,i)=MiE2.dof(Elem.elem2dof(:,i),2);
end
locR=zeros(NT,eQuadN,NNT*NNT);
for i=1:eQuadN
    kappa=Elem.Quad.eQuad.kappa(i);
    psi=Elem.Quad.eQuad.psi(i,:);
    psi_x=Elem.Quad.eQuad.psi_x(i,:);
    psi_y=Elem.Quad.eQuad.psi_y(i,:);
    
    lAt=[qtx*psi_x',qtx*psi_y',qty*psi_x',qty*psi_y'];%(a11,a12;a21,a22)
    detAt=lAt(:,1).*lAt(:,4)-lAt(:,2).*lAt(:,3);
    
    lAt1=[qtx1*psi_x',qtx1*psi_y',qty1*psi_x',qty1*psi_y'];%(a11,a12;a21,a22)
    detAt1=lAt1(:,1).*lAt1(:,4)-lAt1(:,2).*lAt1(:,3);

    lAt2=[qtx2*psi_x',qtx2*psi_y',qty2*psi_x',qty2*psi_y'];%(a11,a12;a21,a22)
    detAt2=lAt2(:,1).*lAt2(:,4)-lAt2(:,2).*lAt2(:,3);

    % lAt^{-\top}=[a22,-a21,-a12,a11]./detlAt
    gradPhix=( lAt(:,4)*psi_x-lAt(:,3)*psi_y)./repmat(detAt,1,NNT);
    gradPhiy=(-lAt(:,2)*psi_x+lAt(:,1)*psi_y)./repmat(detAt,1,NNT);
    
    gradPhix1=( lAt1(:,4)*psi_x-lAt1(:,3)*psi_y)./repmat(detAt1,1,NNT);
    gradPhiy1=(-lAt1(:,2)*psi_x+lAt1(:,1)*psi_y)./repmat(detAt1,1,NNT);
    
    gradPhix2=( lAt2(:,4)*psi_x-lAt2(:,3)*psi_y)./repmat(detAt2,1,NNT);
    gradPhiy2=(-lAt2(:,2)*psi_x+lAt2(:,1)*psi_y)./repmat(detAt2,1,NNT);
    
    for j=1:NNT
        for k=1:NNT
            
            % [ \pa_x (u_1 (w_1-u^n))+ \pa_y (u_1 (w_2-v^n)) + 1/2 (u^n_x + v^n_y) u_1 ]v_1 + ...
            % [ \pa_x (u_2 (w_1-u^n))+ \pa_y (u_2 (w_2-v^n)) + 1/2 (u^n_x + v^n_y) u_2 ]v_2
            
            % w is given (mv), v is test fcn, u is test fcn
            % V (j),
            % U (k)
            % It will be zero if it's (1,0):(0,1) or (0,1):(1,0)
            
            if Sol.Temam==1            
                locR(:,i,(j-1)*NNT+k)=locR(:,i,(j-1)*NNT+k)+...
                    (1.5*gradPhix1(:,k).*mv1_qt(:,i,1).*(detAt1)-0.5*gradPhix2(:,k).*mv2_qt(:,i,1).*(detAt2)...
                    -gradPhix(:,k).*u_qt(:,i,1).*(detAt)+...
                     1.5*gradPhiy1(:,k).*mv1_qt(:,i,2).*(detAt1)-0.5*gradPhiy2(:,k).*mv2_qt(:,i,2).*(detAt2)...
                     -gradPhiy(:,k).*u_qt(:,i,2).*(detAt)+...
                     1.5*psi(k)*mvx1_qt(:,i,1).*(detAt1)-0.5*psi(k)*mvx2_qt(:,i,1).*(detAt2)...
                     -psi(k)*ux_qt(:,i,1).*(detAt)...
                            +1.5*psi(k)*mvy1_qt(:,i,2).*(detAt1)-0.5*psi(k)*mvy2_qt(:,i,2).*(detAt2)...
                            -psi(k)*uy_qt(:,i,2).*(detAt)...
                            +1/2*psi(k)*(ux_qt(:,i,1)+uy_qt(:,i,2)).*(detAt)).*... % This line is the Temam Operator
                    (psi(j)*kappa);
            elseif Sol.Temam==0
                locR(:,i,(j-1)*NNT+k)=locR(:,i,(j-1)*NNT+k)+...
                    (1.5*gradPhix1(:,k).*mv1_qt(:,i,1).*(detAt1)-0.5*gradPhix2(:,k).*mv2_qt(:,i,1).*(detAt2)...
                    -gradPhix(:,k).*u_qt(:,i,1).*(detAt)+...
                     1.5*gradPhiy1(:,k).*mv1_qt(:,i,2).*(detAt1)-0.5*gradPhiy2(:,k).*mv2_qt(:,i,2).*(detAt2)...
                     -gradPhiy(:,k).*u_qt(:,i,2).*(detAt)+...
                     1.5*psi(k)*mvx1_qt(:,i,1).*(detAt1)-0.5*psi(k)*mvx2_qt(:,i,1).*(detAt2)...
                     -psi(k)*ux_qt(:,i,1).*(detAt)...
                            +1.5*psi(k)*mvy1_qt(:,i,2).*(detAt1)-0.5*psi(k)*mvy2_qt(:,i,2).*(detAt2)...
                            -psi(k)*uy_qt(:,i,2).*(detAt)).*...
                    (psi(j)*kappa);
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