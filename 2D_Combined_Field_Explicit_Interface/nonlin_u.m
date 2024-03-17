function fu=nonlin_u(U_qt,Ux_qt,Uy_qt,DeE,Elem)
% calculate <ugradu, psi>

NNT=Elem.nnt; % is the number of nodes on each triangle
NT = Elem.nElem;
eQuad=Elem.Quad.eQuad;

% use ugradu
ugradu_qt=U_qt(:,:,1).*Ux_qt(:,:,1)+U_qt(:,:,2).*Uy_qt(:,:,1);
ugradv_qt=U_qt(:,:,1).*Ux_qt(:,:,2)+U_qt(:,:,2).*Uy_qt(:,:,2);
% % use (\curl\uu)\times(\uu-\grad q)
% ugradu_qt=-(v_qt-qy_qt).*(vx_qt-uy_qt);
% ugradv_qt=(u_qt-qx_qt).*(vx_qt-uy_qt);

fut=zeros(NT,NNT,2); % velocity forcing term due to <Ugradu,psi>, <Ugradv,psi>
kappaPsi=eQuad.psi.*repmat(eQuad.kappa,1,NNT);
for i=1:NNT
    fut(Elem.inElem,i,1)=ugradu_qt(Elem.inElem,:)*kappaPsi(:,i).*abs(DeE.detT);
    fut(Elem.inElem,i,2)=ugradv_qt(Elem.inElem,:)*kappaPsi(:,i).*abs(DeE.detT);
end
%--------------------------------------------------------------------------
% handling isoparametric elements
%--------------------------------------------------------------------------
for i=1:Elem.nBDElem
    eId=Elem.bdElem(i); % index of the element
    fut(eId,:,1)=eQuad.psi'*(ugradu_qt(eId,:)'.*...
        eQuad.kappa.*DeE.bdAbsDF(1:eQuad.N,i));
    fut(eId,:,2)=eQuad.psi'*(ugradv_qt(eId,:)'.*...
        eQuad.kappa.*DeE.bdAbsDF(1:eQuad.N,i));
end
fu=zeros(Elem.nDof,2);
for i=1:2
    fu(:,i) = accumarray(Elem.elem2dof(:),reshape(fut(:,:,i),[],1),[Elem.nDof 1]);
end