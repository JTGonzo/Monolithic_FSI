function fq=fq_divu(u_qt,v_qt,DeE,Elem)
% calculate <U, \grad phi>
NNT=Elem.nnt; % is the number of nodes on each triangle
NT=Elem.nElem;
eQuad=Elem.Quad.eQuad;

fqt=zeros(NT,NNT);   % q forcing term due to <U, grad psi>
kappaPsix=eQuad.psi_x.*repmat(eQuad.kappa,1,NNT);
kappaPsiy=eQuad.psi_y.*repmat(eQuad.kappa,1,NNT);
for i=1:NNT
    fqt(Elem.inElem,i)= ...
         v_qt(Elem.inElem,:)*kappaPsix(:,i).*DeE.p3p1(:,1)...
        +v_qt(Elem.inElem,:)*kappaPsiy(:,i).*DeE.p1p2(:,1)...
        -u_qt(Elem.inElem,:)*kappaPsix(:,i).*DeE.p3p1(:,2)...
        -u_qt(Elem.inElem,:)*kappaPsiy(:,i).*DeE.p1p2(:,2);
end
%--------------------------------------------------------------------------
% handling isoparametric elements
%--------------------------------------------------------------------------
for i=1:Elem.nBDElem
    eId=Elem.bdElem(i); % index of the element
    for k=1:eQuad.N
        fqt(eId,:) = fqt(eId,:) + eQuad.kappa(k)*DeE.bdAbsDF(k,i)*...
            ([u_qt(eId,k),v_qt(eId,k)]*DeE.bdInvDF(:,2*k-1:2*k,i)*...
            [eQuad.psi_x(k,:);eQuad.psi_y(k,:)]);
    end
end

fq = accumarray(Elem.elem2dof(:),fqt(:),[Elem.nDof 1]);