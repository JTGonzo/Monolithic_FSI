function strea=streamFunction(DeE,Elem,Sol)
% dof,U_qt,deg,elem2dof,inElem,bdElem,bdiDF,...
%         bdAbsDF,eQuad,p1p2,p3p1,Ndof)
global GM    
% streamline function

set0index=1;

% fq_divu calculate <U,grad phi>. We now let u->v, v->-u;
fq=fq_divu(Sol.U_qt(:,:,2,1),-Sol.U_qt(:,:,1,1),DeE,Elem);
fq(set0index)=0;

s1=GM.fS(set0index,:); s2=GM.fS(:,set0index);
GM.fS(set0index,:)=0; GM.fS(:,set0index)=0; GM.fS(set0index,set0index)=1;
strea=GM.fS\fq;
GM.fS(set0index,:)=s1; GM.fS(:,set0index)=s2;    