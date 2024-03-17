function [rr,rr1]=SurfaceTraction(Elem,Neu,Phy,mix,u,p,DeEF,meshString)
% lift and drag coefficients for flow past a cylinder with flag

    delta=Elem.Quad.bdQuad.delta(i);
    psi_x=Elem.Quad.bdQuad.psi_x(i,:);
    psi_y=Elem.Quad.bdQuad.psi_y(i,:);
    lAt=[qtx*psi_x',qtx*psi_y',qty*psi_x',qty*psi_y'];%(a11,a12;a21,a22)
    detlAt=lAt(:,1).*lAt(:,4)-lAt(:,2).*lAt(:,3);
    % lAt^{-\top}=[a22,-a21,-a12,a11]./detlAt
    gradPhix=( lAt(:,4)*psi_x-lAt(:,3)*psi_y)./repmat(detlAt,1,NNT);
    gradPhiy=(-lAt(:,2)*psi_x+lAt(:,1)*psi_y)./repmat(detlAt,1,NNT);
    F=[sigma.*psi];
    
    % r=sum_i P_i phi_i, dr/dt=dr/dx-dr/dy since t=x and y=1-x  
    tau=[lAt(:,1)-lAt(:,2),lAt(:,3)-lAt(:,4)];
    jacobian=sqrt(sum(tau.^2,2));
    nn=[-tau(:,2),tau(:,1)]./[jacobian,jacobian]; % outward normal   
    
    bb = zeros(length(ElemF.bdNeuElem),1);
    
    x(i) = ElemF.dof(ElemF.elem2dof(ElemF.bdElem,
    [psi,~,~]=basefun(ElemF.deg,x,y);
    bb = bb+psi;