function b=fU_gradp(DeE,Elem,ph)
% Assembing <grad p, phi>
NNT=Elem.nnt; % is the number of nodes on each triangle
NT=Elem.nElem;
eQuad=Elem.Quad.eQuad;
NinElem=length(Elem.inElem);
NbdElem=Elem.nBDElem;

bt1=zeros(NT,NNT); bt2=bt1;
phNode=zeros(NinElem,NNT); % only used on Triangle elements
phx_quadpt=zeros(NinElem,eQuad.N); phy_quadpt=phx_quadpt;
for i=1:NNT; phNode(:,i)=ph(Elem.elem2dof(Elem.inElem,i)); end
for i=1:eQuad.N
    temp=phNode*[eQuad.psi_x(i,:)',eQuad.psi_y(i,:)'];
    phx_quadpt(:,i) = -(DeE.p3p1(:,2).*temp(:,1)+DeE.p1p2(:,2).*temp(:,2))./DeE.detT;
    phy_quadpt(:,i) =  (DeE.p3p1(:,1).*temp(:,1)+DeE.p1p2(:,1).*temp(:,2))./DeE.detT;
end
kappaPsi=eQuad.psi.*repmat(eQuad.kappa,1,NNT);
for i=1:NNT
    bt1(Elem.inElem,i)=phx_quadpt*kappaPsi(:,i).*abs(DeE.detT);
    bt2(Elem.inElem,i)=phy_quadpt*kappaPsi(:,i).*abs(DeE.detT);    
end

phNode=zeros(NbdElem,NNT); % only used on isoparam. elements
for i=1:NNT; phNode(:,i)=ph(Elem.elem2dof(Elem.bdElem,i)); end
gradph_quadpt=zeros(eQuad.N,2);
for i=1:NbdElem
    for k=1:eQuad.N
        % DF^{-\tau} \grad (\sum_j \psi_j p_j) 
        gradph_quadpt(k,:)=DeE.bdInvDF(:,2*k-1:2*k,i)*...
            ([eQuad.psi_x(k,:);eQuad.psi_y(k,:)]*phNode(i,:)');
    end
    bt1(Elem.bdElem(i),:)=kappaPsi'*(gradph_quadpt(:,1).*DeE.bdAbsDF(1:eQuad.N,i));
    bt2(Elem.bdElem(i),:)=kappaPsi'*(gradph_quadpt(:,2).*DeE.bdAbsDF(1:eQuad.N,i));
end
b = [accumarray(Elem.elem2dof(:),bt1(:),[Elem.nDof 1]),...
     accumarray(Elem.elem2dof(:),bt2(:),[Elem.nDof 1])];