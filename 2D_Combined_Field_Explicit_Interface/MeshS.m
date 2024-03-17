function E=MeshS(Elem,dof)
% return value contains E.nVertex, N.freeN, E.fixedS
% elastic stiffness matrix for the fluid. Use P1.

deg =1;

Elem.dof = dof;
% Elem=fixedMeshData(deg,meshString,nGlobalRefine);
Elem.p1p2 = Elem.dof(Elem.elem2dof(:,2),:)-Elem.dof(Elem.elem2dof(:,1),:); 
Elem.p3p1 = Elem.dof(Elem.elem2dof(:,1),:)-Elem.dof(Elem.elem2dof(:,3),:); 
Elem.detT = Elem.p1p2(:,2).*Elem.p3p1(:,1)-Elem.p1p2(:,1).*Elem.p3p1(:,2);

% Since it is P1, there is no iso element
NNT=Elem.nnt; % is the number of nodes on each triangle
strBaseName=['base',num2str(deg),'.mat'];
load(strBaseName);
N = Elem.nDof;
A = sparse(N,N);  B = sparse(N,N);
% C1 is <phi_x,psi_x>, C2 is <phi_y, psi_x>, C3 is <phi_y, psi_y> 
C1 = sparse(N,N); C2 = sparse(N,N); C3 = sparse(N,N);
a = [Elem.p3p1(:,1).^2+Elem.p3p1(:,2).^2,...
    Elem.p1p2(:,2).*Elem.p3p1(:,2)+Elem.p1p2(:,1).*Elem.p3p1(:,1),...
    Elem.p1p2(:,1).^2+Elem.p1p2(:,2).^2];
tau=(1+(max(Elem.detT)-min(Elem.detT))./Elem.detT);
for i = 1:NNT
    for j = 1:NNT
        % triangle
        Aij = (UxUx(i,j)*a(:,1)+(UxUy(i,j)+UxUy(j,i))*a(:,2)...
              +UyUy(i,j)*a(:,3))./abs(Elem.detT).*tau;
        C1ij = (UxUx(i,j)*(Elem.p3p1(:,2).^2)+(UxUy(i,j)+UxUy(j,i))*(Elem.p3p1(:,2).*Elem.p1p2(:,2))...
               +UyUy(i,j)*(Elem.p1p2(:,2).^2))./abs(Elem.detT).*tau;
        C2ij = (-UxUx(i,j)*(Elem.p3p1(:,1).*Elem.p3p1(:,2))-UxUy(i,j)*(Elem.p3p1(:,1).*Elem.p1p2(:,2))...
                -UxUy(j,i)*(Elem.p3p1(:,2).*Elem.p1p2(:,1))-UyUy(i,j)*(Elem.p1p2(:,1).*Elem.p1p2(:,2)))./abs(Elem.detT).*tau;   
        C3ij = (UxUx(i,j)*(Elem.p3p1(:,1).^2)+(UxUy(i,j)+UxUy(j,i))*(Elem.p3p1(:,1).*Elem.p1p2(:,1))...
               +UyUy(i,j)*(Elem.p1p2(:,1).^2))./abs(Elem.detT).*tau;

           
        A = A + sparse(Elem.elem2dof(:,i),Elem.elem2dof(:,j),Aij,N,N);
        C1 = C1 + sparse(Elem.elem2dof(:,i),Elem.elem2dof(:,j),C1ij,N,N);
        C2 = C2 + sparse(Elem.elem2dof(:,i),Elem.elem2dof(:,j),C2ij,N,N);
        C3 = C3 + sparse(Elem.elem2dof(:,i),Elem.elem2dof(:,j),C3ij,N,N);        
    end
end
E.Elem = Elem;
E.fixedS=[A+C1,C2;C2',A+C3]+[C1,C2';C2,C3];
E.nVertex=Elem.nVertex;
% the freeN of E is not the same as that of Elem
E.fFreeN = setdiff((1:E.nVertex)',[Elem.bdNeuN;Elem.bdDirN]);