function [Linfe,L2e]=Error(u,uexact,node,elem2node,det)
[w,q]= legendre();
Linfe=max(abs(u-uexact(node(:,1),node(:,2))));
L2e=0;
% H1e=0;
uh = [u(elem2node(:,1),:), u(elem2node(:,2),:), u(elem2node(:,3),:)];

for i=1:size(q,1)   
    phi=evalPhi(q(i,1),q(i,2)); 
    qt = node(elem2node(:,1),:)*phi(1,1)+node(elem2node(:,2),:)*phi(1,2)+ node(elem2node(:,3),:)*phi(1,3); 
    uex_qt=uexact(qt(:,1),qt(:,2));
    uh_qt = uh*phi';
%       uhx_qt = sum (gradP.*uh,2);
    L2e = L2e+0.5*sum(w(i)*(uex_qt-uh_qt).^2.*abs(det));
%       H1e = H1e+sum(w(i)*(uexx_qt-uhx_qt).^2.*det);
    Linfe=max(Linfe, max(abs(uex_qt-uh_qt)));
end
L2e=sqrt(L2e);

function [w,pt] = legendre()
pt = [1/6,1/6;1/6,2/3;2/3,1/6];
w=[1/3,1/3,1/3];