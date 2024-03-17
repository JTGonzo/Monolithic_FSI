function vort=vorticityFunction(U,deg,elem2dof,inElem,bdElem,...
    bdiDF,bdAbsDF,eQuad,p1p2,p3p1,detT,Ndof)
global M    
% vorticity function
%inputs:
%outputs:
b1=fU_gradp(deg,elem2dof,inElem,bdElem,bdiDF,bdAbsDF,...
    eQuad,p1p2,p3p1,detT,U(:,1),Ndof);
b2=fU_gradp(deg,elem2dof,inElem,bdElem,bdiDF,bdAbsDF,...
    eQuad,p1p2,p3p1,detT,U(:,2),Ndof);
vort=M\(b2(:,1)-b1(:,2));