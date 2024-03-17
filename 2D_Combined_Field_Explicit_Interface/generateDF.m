function [bdInvDF,bdAbsDF]=generateDF(elem2dof,bdElem,dof,eQuad,bdQuad)
%--------------------------------------------------------------------------
% Pre-compute the DF for those elem where F is the isoparametric transf.
% DF should be used together with bdElem, gamI
% Assuming no element has two boundary edges     
% bdInvDF is the DF^{-\top} at all the inter-quad and bd-quad pts 
% for the iso-elem touching the boundary
NbdE=length(bdElem); 
bdInvDF=zeros(2,2*(eQuad.N+bdQuad.N),NbdE);
bdAbsDF=zeros(eQuad.N+bdQuad.N,NbdE);
for i=1:NbdE
    dofId=elem2dof(bdElem(i),:); % index of the dof on this element
    % quadptT=psi*dof(dofId,:); % coordinates of the inter quadrature pts
    % qdadptE=bdpsi*dof(dofId,:); % coordinates of the bd quadrature pts
    for k=1:eQuad.N
        DF = [eQuad.psi_x(k,:);eQuad.psi_y(k,:)]*dof(dofId,:); % it's transp. of DF
        bdAbsDF(k,i) = abs(det(DF)); 
        bdInvDF(:,2*k-1:2*k,i)=inv(DF);
    end
    for k=1:bdQuad.N
        DF = [bdQuad.psi_x(k,:);bdQuad.psi_y(k,:)]*dof(dofId,:); % it's transp. of DF
        bdAbsDF(eQuad.N+k,i) = abs(det(DF));
        bdInvDF(:,2*eQuad.N+(2*k-1:2*k),i)=inv(DF);        
    end
end