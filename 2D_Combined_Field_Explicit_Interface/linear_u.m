function fu=linear_u(U_qt,DeE,Elem)
% calculate <u, psi>

NNT=Elem.nnt; % is the number of nodes on each triangle
NT=Elem.nElem; 
eQuad=Elem.Quad.eQuad;

fut=zeros(NT,NNT,2); % vel forcing term due to <Ugradu,psi>, <Ugradv,psi>
kappaPsi=eQuad.psi.*repmat(eQuad.kappa,1,NNT);
for i=1:NNT
    fut(Elem.inElem,i,1)=U_qt(Elem.inElem,:,1)*kappaPsi(:,i).*abs(DeE.detT);
    fut(Elem.inElem,i,2)=U_qt(Elem.inElem,:,2)*kappaPsi(:,i).*abs(DeE.detT);
end
%--------------------------------------------------------------------------
% handling isoparametric elements
%--------------------------------------------------------------------------
for i=1:length(Elem.bdElem)
    eId=Elem.bdElem(i); % index of the element
    fut(eId,:,1)=eQuad.psi'*(U_qt(eId,:,1)'.*eQuad.kappa.*...
        DeE.bdAbsDF(1:eQuad.N,i));
    fut(eId,:,2)=eQuad.psi'*(U_qt(eId,:,2)'.*eQuad.kappa.*...
        DeE.bdAbsDF(1:eQuad.N,i));
end
fu=zeros(Elem.nDof,2);
for i=1:2
    fu(:,i) = accumarray(Elem.elem2dof(:),reshape(fut(:,:,i),[],1),[Elem.nDof 1]);
end