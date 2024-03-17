figure(4); clf;
strea=streamFunction(DeEF{1},ElemF,Sol);
%lab=[linspace(-2,-0.2,20),linspace(-0.195,-0.05,8),linspace(-0.0495,0,5)];
a=5;
lab=-a.^(linspace(-6,log(1.3)/log(a),20));
tricontour(DeEF{1}.dof,ElemF.elem2dof(:,1:3),strea,lab);
% clabel(c,h);
showmesh(ElemS.dof(1:ElemS.nVertex,:)+Sol.eta(1:ElemS.nVertex,:,1),...
    ElemS.elem2dof(:,1:3));
axis equal;  axis([0.5,3,-0.5,0.3]);
drawnow;
