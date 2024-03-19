%% Coding Project 2 
%%
% declare the global parameters
clear all; close all; clc
disp('test')

tic;
Sol.delt = 0.001;

%% Defining Fluid Properties
U_0 = 1.0;                             % Flow Velocity
Phy.rhof = 1000;                       % Fluid Density
%Phy.fmu = U_0*Phy.rhof*0.1/Re;        % Dynamic Viscosity 
%Re = 20;                              % Reynold's Number
Phy.fmu = 1.0;
Re = U_0*Phy.rhof*1/Phy.fmu;

%% Defining Material Properties
Phy.rhos = 1000;                                % Solid Density
E = 5.6e6;                                      % Young's Modulus
snu = 0.4;                                      % Poisson's Ratio
Phy.lambda = snu * E / ((1-2*snu)*(1+snu));     % %Lame's Constants
Phy.smu = E/(2*(1+snu));                        % Shear Modulus

%% Readign Input Files
wrkDir = './' ;
problemString = 'FlexFil' ;
elemType = '4Quad' ;
problemType = '2D' ;

% Coordinate, connectivity and boundary data
dataStr = strcat(wrkDir,'Data_flexible_filament.mat') ;
load(dataStr);
ndof = size(crd,1) ;

coordinates = crd(:,2:3);
%nodes = conn;
nodes = [BCFluid; BCStructure];

%% Converting to a 3D triangular mesh
% converting from quad to tri elements
tri = Qua4toTri3(nodes);
triF = Qua4toTri3(BCFluid);
triS = Qua4toTri3(BCStructure);
plotMesh(coordinates,tri)

% converting the element and corodinates to 3D
tri3D = [tri, tri+length(coordinates)];
temp1 = crd(:,1)+length(crd);
temp2 = crd(:,4)+0.001;
crd3D = [crd;[temp1,crd(:,2:3),temp2]];

%% Converting BC to nodal format
fluidNodes = [unique(BCFluid);unique(BCFluid)+length(crd)];
solidNodes = [unique(BCStructure);unique(BCStructure)+length(crd)];
int = intersect(fluidNodes,solidNodes);
side1Nodes = crd(:,1);
side2Nodes = crd(:,1)+ length(crd);
itfNodes = intersect(fluidNodes,solidNodes);
inletNodes = [unique(BCLeft);unique(BCLeft)+length(crd)];
outletNodes = [unique(BCRight);unique(BCRight)+length(crd)];
topNodes = [unique(BCTop);unique(BCTop)+length(crd)];
bottomNodes = [unique(BCBottom);unique(BCBottom)+length(crd)];

%% Formatting Additional Variables
Elem.dof = crd3D(:,2:end);
Elem.nDof = length(Elem.dof);
ElemF.elem2dof = [triF, triF+length(coordinates)];
ElemF.ndof = length(ElemF.elem2dof);
ElemS.elem2dof = [triS, triS+length(coordinates)];
ElemS.ndof = length(ElemF.elem2dof);

Sol.u = zeros(Elem.nDof,3,4);
Sol.p = zeros(Elem.nDof,1);
disp = zeros(Elem.nDof,3,3);
hinged_1 = find((Elem.dof(:,1)).^2+(Elem.dof(:,2)).^2<=0.50001^2 & (abs(Elem.dof(:,3)-20)<1E-8 | abs(Elem.dof(:,3))<1E-8));
    
nodeId = find(abs(Elem.dof(:,1))<1E-8 & abs(Elem.dof(:,2))<1E-8 & abs(Elem.dof(:,3)-10)<1E-8);

tipDisp = [];

%% Building matrices for flow domain 

%% Compute geometric quantities and gradient of local basis

points = [3.3333333333333331E-01,  3.3333333333333331E-01, -5.7735026918962584E-01
    1.3333333333333333E+00,  3.3333333333333331E-01, -5.7735026918962584E-01
    3.3333333333333331E-01,  1.3333333333333333E+00, -5.7735026918962584E-01
    3.3333333333333331E-01,  3.3333333333333331E-01,  5.7735026918962584E-01
    1.3333333333333333E+00,  3.3333333333333331E-01,  5.7735026918962584E-01
    3.3333333333333331E-01,  1.3333333333333333E+00,  5.7735026918962584E-01];
w = [6.6666666666666667E-01,  6.6666666666666667E-01,  6.6666666666666667E-01,6.6666666666666667E-01,  6.6666666666666667E-01,  6.6666666666666667E-01];

nQuad = length(w);

phi(:,1) = 1/4*(2-points(:,1)-points(:,2)).*(1-points(:,3));
phi(:,2) = 1/4*(points(:,1)).*(1-points(:,3));
phi(:,3) = 1/4*(points(:,2)).*(1-points(:,3));
phi(:,4) = 1/4*(2-points(:,1)-points(:,2)).*(1+points(:,3));
phi(:,5) = 1/4*(points(:,1)).*(1+points(:,3));
phi(:,6) = 1/4*(points(:,2)).*(1+points(:,3));

phix(:,1) =  -1/4.*(1-points(:,3));
phix(:,2) =   1/4.*(1-points(:,3));
phix(:,3) =   zeros(6,1);
phix(:,4) =  -1/4.*(1+points(:,3));
phix(:,5) =   1/4.*(1+points(:,3));
phix(:,6) =   zeros(6,1);

phiy(:,1) =  -1/4.*(1-points(:,3));
phiy(:,2) =   zeros(6,1);
phiy(:,3) =   1/4.*(1-points(:,3));
phiy(:,4) =  -1/4.*(1+points(:,3));
phiy(:,5) =   zeros(6,1);
phiy(:,6) =   1/4.*(1+points(:,3));

phiz(:,1) = -1/4*(2-points(:,1)-points(:,2));
phiz(:,2) = -1/4*(points(:,1));
phiz(:,3) = -1/4*(points(:,2));
phiz(:,4) =  1/4*(2-points(:,1)-points(:,2));
phiz(:,5) =  1/4*(points(:,1));
phiz(:,6) =  1/4*(points(:,2));

phix = phix';
phiy = phiy';
phiz = phiz';

nnt = size(ElemF.elem2dof,2);
ElemF.nnt = nnt;
ElemS.nnt = nnt;
ElemF.nElem = size(ElemF.elem2dof,1);
ElemS.nElem = size(ElemS.elem2dof,1);
Elem.nDof = size(Elem.dof,1);

iif = zeros(ElemF.nnt^2*ElemF.nElem,1); jjf = zeros(ElemF.nnt^2*ElemF.nElem,1);
index = 0;
for i = 1:nnt
    for j = 1:nnt
        iif(index+1:index+ElemF.nElem) = double(ElemF.elem2dof(:,i)); 
        jjf(index+1:index+ElemF.nElem) = double(ElemF.elem2dof(:,j));  
        index = index + ElemF.nElem;
    end
end

iis = zeros(ElemS.nnt^2*ElemS.nElem,1); jjs = zeros(ElemS.nnt^2*ElemS.nElem,1);
index = 0;
for i = 1:nnt
    for j = 1:nnt
        iis(index+1:index+ElemS.nElem) = double(ElemS.elem2dof(:,i)); 
        jjs(index+1:index+ElemS.nElem) = double(ElemS.elem2dof(:,j));  
        index = index + ElemS.nElem;
    end
end

for timeStep = 1:4000
    
    timeStep
    
disp(solidNodes,:,1) = disp(solidNodes,:,2) + Sol.delt * (1.5*Sol.u(solidNodes,:,2) - 0.5 * Sol.u(solidNodes,:,3));

aleDisp = zeros(Elem.nDof,3);

%% ALE MESH

aleDisp(itfNodes,:) = disp(itfNodes,:,1);

xxf = zeros(size(ElemF.elem2dof));
yyf = zeros(size(ElemF.elem2dof));
zzf = zeros(size(ElemF.elem2dof));

for i=1:ElemF.nnt
    xxf(:,i) = Elem.dof(ElemF.elem2dof(:,i),1);
    yyf(:,i) = Elem.dof(ElemF.elem2dof(:,i),2);
    zzf(:,i) = Elem.dof(ElemF.elem2dof(:,i),3);
end


sA1 = zeros(nnt^2*ElemF.nElem,nQuad);
sA2 = zeros(nnt^2*ElemF.nElem,nQuad);
sA3 = zeros(nnt^2*ElemF.nElem,nQuad);
sA4 = zeros(nnt^2*ElemF.nElem,nQuad);
sA5 = zeros(nnt^2*ElemF.nElem,nQuad);
sA6 = zeros(nnt^2*ElemF.nElem,nQuad);

for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
        index = 0;
    for i = 1:ElemF.nnt
        for j = 1:ElemF.nnt
            Aij_1 = w(p)*(DphiDx(:,i).*DphiDx(:,j));
            Aij_2 = w(p)*(DphiDy(:,i).*DphiDy(:,j));
            Aij_3 = w(p)*(DphiDz(:,i).*DphiDz(:,j));
            Aij_4 = w(p)*(DphiDx(:,i).*DphiDy(:,j));
            Aij_5 = w(p)*(DphiDy(:,i).*DphiDz(:,j));
            Aij_6 = w(p)*(DphiDz(:,i).*DphiDx(:,j));
            Aij_1 = Aij_1.*volume;
            Aij_2 = Aij_2.*volume;
            Aij_3 = Aij_3.*volume;
            Aij_4 = Aij_4.*volume;
            Aij_5 = Aij_5.*volume;
            Aij_6 = Aij_6.*volume;
            sA1(index+1:index+ElemF.nElem,p) = Aij_1;
            sA2(index+1:index+ElemF.nElem,p) = Aij_2;
            sA3(index+1:index+ElemF.nElem,p) = Aij_3;
            sA4(index+1:index+ElemF.nElem,p) = Aij_4;
            sA5(index+1:index+ElemF.nElem,p) = Aij_5;
            sA6(index+1:index+ElemF.nElem,p) = Aij_6;
            index = index + ElemF.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA6 = sum(sA6,2);

% assemble the matrix
A1 = sparse(iif,jjf,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iif,jjf,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iif,jjf,sA3,Elem.nDof,Elem.nDof);
A4 = sparse(iif,jjf,sA4,Elem.nDof,Elem.nDof);
A5 = sparse(iif,jjf,sA5,Elem.nDof,Elem.nDof);
A6 = sparse(iif,jjf,sA6,Elem.nDof,Elem.nDof);

KAle = [2*A1+A2+A3 A4' A6; A4 A1+2*A2+A3 A5'; A6' A5 A1+A2+2*A3] + [A1 A4 A6'; A4' A2 A5; A6 A5' A3];

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 Dphip

RHSAle = - KAle * aleDisp(:);

if flag
    freeNodes = unique([inletNodes;outletNodes;topNodes;bottomNodes;itfNodes;solidNodes;hinged_1;side1Nodes;side2Nodes ]);
else
    freeNodes = find(Elem.dof(:,2)==1 | Elem.dof(:,2)== 0 | Elem.dof(:,2)<0 | Elem.dof(:,1)==0 | Elem.dof(:,1)==1 );
end

freeNodes = setdiff(1:size(Elem.dof,1),freeNodes);
freeNodes = [freeNodes';freeNodes' + size(Elem.dof,1)];

disp(freeNodes) = KAle(freeNodes,freeNodes)\RHSAle(freeNodes);

meshVel = (1.5 * disp(:,:,1) - 2 * disp(:,:,2) + 0.5 * disp(:,:,3))/Sol.delt ;

%% FLUID

xxf = zeros(size(ElemF.elem2dof));
yyf = zeros(size(ElemF.elem2dof));
zzf = zeros(size(ElemF.elem2dof));
ux = zeros(size(ElemF.elem2dof));
uy = zeros(size(ElemF.elem2dof));
uz = zeros(size(ElemF.elem2dof));
locDispx = zeros(size(ElemF.elem2dof));
locDispy = zeros(size(ElemF.elem2dof));
locDispz = zeros(size(ElemF.elem2dof));

locDispPrevx = zeros(size(ElemF.elem2dof));
locDispPrevy = zeros(size(ElemF.elem2dof));
locDispPrevz = zeros(size(ElemF.elem2dof));

dof = Elem.dof+disp(:,:,1);

for i=1:ElemF.nnt
    xxf(:,i) = dof(ElemF.elem2dof(:,i),1);
    yyf(:,i) = dof(ElemF.elem2dof(:,i),2);
    zzf(:,i) = dof(ElemF.elem2dof(:,i),3);
    ux(:,i) =  2.25 * Sol.u(ElemF.elem2dof(:,i),1,2) - 1.5 * Sol.u(ElemF.elem2dof(:,i),1,3) + 0.25 * Sol.u(ElemF.elem2dof(:,i),1,4) - meshVel(ElemF.elem2dof(:,i),1);
    uy(:,i) =  2.25 * Sol.u(ElemF.elem2dof(:,i),2,2) - 1.5 * Sol.u(ElemF.elem2dof(:,i),2,3) + 0.25 * Sol.u(ElemF.elem2dof(:,i),2,4) - meshVel(ElemF.elem2dof(:,i),2);
    uz(:,i) =  2.25 * Sol.u(ElemF.elem2dof(:,i),3,2) - 1.5 * Sol.u(ElemF.elem2dof(:,i),3,3) + 0.25 * Sol.u(ElemF.elem2dof(:,i),3,4) - meshVel(ElemF.elem2dof(:,i),3);
end

%% Assemble of the time derivative mass matrix
    % Defining the quadrature

    if flag
        dirichletNodes = unique([inletNodes;topNodes;bottomNodes ]);
    else
        dirichletNodes = find(abs(Elem.dof(:,1))<1E-8 & Elem.dof(:,2)>0 & Elem.dof(:,2)<1);
    end
    
%     if (timeStep*Sol.delt) <=2
%         Sol.u(dirichletNodes,1,1) = (1-cos(pi/2*timeStep*Sol.delt))/2 * 1.5 * U_0 * 4.0/0.1681 * ( Elem.dof(dirichletNodes,2) .* (0.41 - Elem.dof(dirichletNodes,2)) );
%     else
%         Sol.u(dirichletNodes,1,1) = 1.5 * U_0 * 4.0/0.1681 * ( Elem.dof(dirichletNodes,2) .* (0.41 - Elem.dof(dirichletNodes,2)) );
%     end

Sol.u(dirichletNodes,1,1) = 1.5*U_0/100*(Elem.dof(dirichletNodes,3)).*(20-Elem.dof(dirichletNodes,3));
    
% compute non-zeros
sA = [];
for p = 1:nQuad    
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    negJacobian = find(volume<0);
    if ~isempty(negJacobian)
        disp('Mesh deformed, Negative Jacobian');
        exit
    end
    
    index = 0;
    for i = 1:nnt
        for j = 1:nnt
            Mij = w(p)*(phi(p,i)'*phi(p,j));
            Mij = Mij.*volume;
            sA(index+1:index+ElemF.nElem,p) = Mij;
            index = index + ElemF.nElem;
        end
    end
end
sA = sum(sA,2);
% assemble the matrix
Mf = sparse(iif,jjf,sA,Elem.nDof,Elem.nDof);
Mf = Mf';
ZeroF = sparse(Elem.nDof,Elem.nDof);
Mf = [Mf ZeroF ZeroF ZeroF;ZeroF Mf ZeroF ZeroF;ZeroF ZeroF Mf ZeroF; ZeroF ZeroF ZeroF ZeroF];
Mflow = 1.5 * Phy.rhof * Mf/Sol.delt;
MatlabFlow = Mflow;
clear Mij sA

u1 = Sol.u(:,:,1);
u2 = Sol.u(:,:,2);
u3 = Sol.u(:,:,3);
u4 = Sol.u(:,:,4);

u1 = [u1(:);zeros(Elem.nDof,1)];
u2 = [u2(:);zeros(Elem.nDof,1)];
u3 = [u3(:);zeros(Elem.nDof,1)];
u4 = [u4(:);zeros(Elem.nDof,1)];

if (size(Sol.p,1) == 1 & size(Sol.p,2)>1)
    Sol.p = [Sol.p]';
end
    
p1 = [zeros(Elem.nDof*3,1);Sol.p];
RHS = -1.5*Phy.rhof/Sol.delt * Mf * u1(:) + 2 * Phy.rhof/Sol.delt * Mf * u2(:) - 0.5 * Phy.rhof/Sol.delt * Mf * u3(:);

sA1 = zeros(nnt^2*ElemF.nElem,nQuad);
sA2 = zeros(nnt^2*ElemF.nElem,nQuad);
sA3 = zeros(nnt^2*ElemF.nElem,nQuad);

sA4 = zeros(nnt^2*ElemF.nElem,nQuad);
sA5 = zeros(nnt^2*ElemF.nElem,nQuad);
sA6 = zeros(nnt^2*ElemF.nElem,nQuad);


for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
    xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
    xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
    xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
    xsInv(:,9) = 1./xsInv(:,9) ;
    xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
    xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
    xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
    xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
    xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
    xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
    xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
    xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
    xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
    c1 = 1.154700538379252E0 ;
    c2 = 5.773502691896259E-01;
    
    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
    
    locUX  = sum(repmat(phi(p,:),ElemF.nElem,1).*ux,2);
    locUY  = sum(repmat(phi(p,:),ElemF.nElem,1).*uy,2);
    locUZ  = sum(repmat(phi(p,:),ElemF.nElem,1).*uz,2);
    
    VG1		= gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
    VG2		= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
    VG3		= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

    VGV		= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

    GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
        2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
    tauM = Phy.rhof^2 * (2/Sol.delt)^2 + Phy.fmu^2 * GG * 36 + Phy.rhof^2 * abs(VGV);    
    
    tauM = tauM.^(0.5);
    
    tauC = 0.125 * tauM ./(Phy.rhof^2 * (gijDown(:,1) + gijDown(:,2) + gijDown(:,3)));
    
    tauM = 1./tauM;

    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    
        locUX  = sum(repmat(phi(p,:),ElemF.nElem,1).*ux,2);
        locUY  = sum(repmat(phi(p,:),ElemF.nElem,1).*uy,2);
        locUZ  = sum(repmat(phi(p,:),ElemF.nElem,1).*uz,2);
        
        index = 0;
    for i = 1:ElemF.nnt
        for j = 1:ElemF.nnt
            Aij_1 = w(p)*(locUX.*DphiDx(:,i)*phi(p,j));
            Aij_2 = w(p)*(locUY.*DphiDy(:,i)*phi(p,j));
            Aij_3 = w(p)*(locUZ.*DphiDz(:,i)*phi(p,j));
            
            Aij_4 = w(p)*(locUX.*DphiDx(:,j)*phi(p,i));
            Aij_5 = w(p)*(locUY.*DphiDy(:,j)*phi(p,i));
            Aij_6 = w(p)*(locUZ.*DphiDz(:,j)*phi(p,i));
            
            Aij_1 = Aij_1.*volume.*tauM;
            Aij_2 = Aij_2.*volume.*tauM;
            Aij_3 = Aij_3.*volume.*tauM;
            
            Aij_4 = Aij_4.*volume;
            Aij_5 = Aij_5.*volume;
            Aij_6 = Aij_6.*volume;
            
            sA1(index+1:index+ElemF.nElem,p) = Aij_1;
            sA2(index+1:index+ElemF.nElem,p) = Aij_2;
            sA3(index+1:index+ElemF.nElem,p) = Aij_3;
            sA4(index+1:index+ElemF.nElem,p) = Aij_4;
            sA5(index+1:index+ElemF.nElem,p) = Aij_5;
            sA6(index+1:index+ElemF.nElem,p) = Aij_6;
            
            
            index = index + ElemF.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);

sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA6 = sum(sA6,2);

% assemble the matrix
A1 = sparse(iif,jjf,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iif,jjf,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iif,jjf,sA3,Elem.nDof,Elem.nDof);

A4 = sparse(iif,jjf,sA4,Elem.nDof,Elem.nDof);
A5 = sparse(iif,jjf,sA5,Elem.nDof,Elem.nDof);
A6 = sparse(iif,jjf,sA6,Elem.nDof,Elem.nDof);

TauMStabK1 = 1.5*Phy.rhof^2/Sol.delt*[A1+A2+A3 ZeroF ZeroF ZeroF; ZeroF A1+A2+A3 ZeroF ZeroF;ZeroF ZeroF A1+A2+A3 ZeroF; ZeroF ZeroF ZeroF ZeroF];

Conv = Phy.rhof*[A4+A5+A6 ZeroF ZeroF ZeroF; ZeroF A4+A5+A6 ZeroF ZeroF; ZeroF ZeroF A4+A5+A6 ZeroF; ZeroF ZeroF ZeroF ZeroF];

MatlabFlow = MatlabFlow + TauMStabK1 + Conv;

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 Dphip

RHS = RHS - (Conv) * u1(:);

%% Assemble stiffness matrix for viscous stress component
% generate sparse pattern

% Defining the quadrature
% compute non-zeros
sA1 = zeros(nnt^2*ElemF.nElem,nQuad);
sA2 = zeros(nnt^2*ElemF.nElem,nQuad);
sA3 = zeros(nnt^2*ElemF.nElem,nQuad);
sA4 = zeros(nnt^2*ElemF.nElem,nQuad);
sA5 = zeros(nnt^2*ElemF.nElem,nQuad);
sA6 = zeros(nnt^2*ElemF.nElem,nQuad);

for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
       
        index = 0;
    for i = 1:ElemF.nnt
        for j = 1:ElemF.nnt
            Aij_1 = w(p)*(DphiDx(:,i).*DphiDx(:,j));
            Aij_2 = w(p)*(DphiDy(:,i).*DphiDy(:,j));
            Aij_3 = w(p)*(DphiDz(:,i).*DphiDz(:,j));
            Aij_4 = w(p)*(DphiDx(:,i).*DphiDy(:,j));
            Aij_5 = w(p)*(DphiDy(:,i).*DphiDz(:,j));
            Aij_6 = w(p)*(DphiDz(:,i).*DphiDx(:,j));
            Aij_1 = Aij_1.*volume;
            Aij_2 = Aij_2.*volume;
            Aij_3 = Aij_3.*volume;
            Aij_4 = Aij_4.*volume;
            Aij_5 = Aij_5.*volume;
            Aij_6 = Aij_6.*volume;
            sA1(index+1:index+ElemF.nElem,p) = Aij_1;
            sA2(index+1:index+ElemF.nElem,p) = Aij_2;
            sA3(index+1:index+ElemF.nElem,p) = Aij_3;
            sA4(index+1:index+ElemF.nElem,p) = Aij_4;
            sA5(index+1:index+ElemF.nElem,p) = Aij_5;
            sA6(index+1:index+ElemF.nElem,p) = Aij_6;
            index = index + ElemF.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA6 = sum(sA6,2);

% assemble the matrix
A1 = sparse(iif,jjf,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iif,jjf,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iif,jjf,sA3,Elem.nDof,Elem.nDof);
A4 = sparse(iif,jjf,sA4,Elem.nDof,Elem.nDof);
A5 = sparse(iif,jjf,sA5,Elem.nDof,Elem.nDof);
A6 = sparse(iif,jjf,sA6,Elem.nDof,Elem.nDof);

Kf = Phy.fmu*[2*A1+A2+A3 A4' A6 ZeroF; A4 A1+2*A2+A3 A5' ZeroF; A6' A5 A1+A2+2*A3 ZeroF; ZeroF ZeroF ZeroF ZeroF];

MatlabFlow = MatlabFlow + Kf;

RHS = RHS - Kf * u1(:);

% Kf = Muf*[2*A1+A2+A3 A4' A6; A4 A1+2*A2+A3 A5'; A6' A5 A1+A2+2*A3]+lambda*[A1 A4 A6'; A4' A2 A5; A6 A5' A3];

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 Dphip

sA1 = zeros(nnt^2*ElemF.nElem,nQuad);
sA2 = zeros(nnt^2*ElemF.nElem,nQuad);
sA3 = zeros(nnt^2*ElemF.nElem,nQuad);
sA4 = zeros(nnt^2*ElemF.nElem,nQuad);
sA5 = zeros(nnt^2*ElemF.nElem,nQuad);
sA6 = zeros(nnt^2*ElemF.nElem,nQuad);

for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
    xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
    xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
    xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
    xsInv(:,9) = 1./xsInv(:,9) ;
    xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
    xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
    xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
    xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
    xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
    xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
    xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
    xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
    xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
    c1 = 1.154700538379252E0 ;
    c2 = 5.773502691896259E-01;
    
    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
    
    locUX  = sum(repmat(phi(p,:),ElemF.nElem,1).*ux,2);
    locUY  = sum(repmat(phi(p,:),ElemF.nElem,1).*uy,2);
    locUZ  = sum(repmat(phi(p,:),ElemF.nElem,1).*uz,2);
    
    VG1		= gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
    VG2		= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
    VG3		= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

    VGV		= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

    GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
        2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
    tauM = Phy.rhof^2 * (2/Sol.delt)^2 + Phy.fmu^2 * GG * 36 + Phy.rhof^2 * abs(VGV);
    
    tauM = tauM.^(0.5);
    
    tauC = 0.125 * tauM ./(Phy.rhof^2 * (gijDown(:,1) + gijDown(:,2) + gijDown(:,3)));
    
    tauM = 1./tauM;

    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    
        index = 0;
    for i = 1:ElemF.nnt
        for j = 1:ElemF.nnt
            Aij_1 = w(p)*(DphiDx(:,i).*DphiDx(:,j));
            Aij_2 = w(p)*(DphiDy(:,i).*DphiDy(:,j));
            Aij_3 = w(p)*(DphiDz(:,i).*DphiDz(:,j));
            Aij_4 = w(p)*(DphiDx(:,i).*DphiDy(:,j));
            Aij_5 = w(p)*(DphiDy(:,i).*DphiDz(:,j));
            Aij_6 = w(p)*(DphiDz(:,i).*DphiDx(:,j));
            Aij_1 = Aij_1.*volume.*tauC;
            Aij_2 = Aij_2.*volume.*tauC;
            Aij_3 = Aij_3.*volume.*tauC;
            Aij_4 = Aij_4.*volume.*tauC;
            Aij_5 = Aij_5.*volume.*tauC;
            Aij_6 = Aij_6.*volume.*tauC;
            sA1(index+1:index+ElemF.nElem,p) = Aij_1;
            sA2(index+1:index+ElemF.nElem,p) = Aij_2;
            sA3(index+1:index+ElemF.nElem,p) = Aij_3;
            sA4(index+1:index+ElemF.nElem,p) = Aij_4;
            sA5(index+1:index+ElemF.nElem,p) = Aij_5;
            sA6(index+1:index+ElemF.nElem,p) = Aij_6;
            index = index + ElemF.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA6 = sum(sA6,2);

% assemble the matrix
A1 = sparse(iif,jjf,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iif,jjf,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iif,jjf,sA3,Elem.nDof,Elem.nDof);
A4 = sparse(iif,jjf,sA4,Elem.nDof,Elem.nDof);
A5 = sparse(iif,jjf,sA5,Elem.nDof,Elem.nDof);
A6 = sparse(iif,jjf,sA6,Elem.nDof,Elem.nDof);

TauCStab = [A1 A4 A6' ZeroF; A4' A2 A5 ZeroF;A6 A5' A3 ZeroF; ZeroF ZeroF ZeroF ZeroF];

MatlabFlow = MatlabFlow + TauCStab;

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 Dphip

RHS = RHS - TauCStab*u1(:);

RHS = RHS - TauMStabK1 * u1(:);

RHS = RHS + 2/1.5* TauMStabK1*u2(:) - 0.5/1.5 * TauMStabK1*u3(:);


sA1 = zeros(nnt^2*ElemF.nElem,nQuad);
sA2 = zeros(nnt^2*ElemF.nElem,nQuad);
sA3 = zeros(nnt^2*ElemF.nElem,nQuad);
sA4 = zeros(nnt^2*ElemF.nElem,nQuad);
sA5 = zeros(nnt^2*ElemF.nElem,nQuad);
sA6 = zeros(nnt^2*ElemF.nElem,nQuad);
sA7 = zeros(nnt^2*ElemF.nElem,nQuad);
sA8 = zeros(nnt^2*ElemF.nElem,nQuad);
sA9 = zeros(nnt^2*ElemF.nElem,nQuad);

for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
    xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
    xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
    xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
    xsInv(:,9) = 1./xsInv(:,9) ;
    xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
    xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
    xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
    xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
    xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
    xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
    xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
    xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
    xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
    c1 = 1.154700538379252E0 ;
    c2 = 5.773502691896259E-01;
    
    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
    
    locUX  = sum(repmat(phi(p,:),ElemF.nElem,1).*ux,2);
    locUY  = sum(repmat(phi(p,:),ElemF.nElem,1).*uy,2);
    locUZ  = sum(repmat(phi(p,:),ElemF.nElem,1).*uz,2);
    
    VG1		= gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
    VG2		= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
    VG3		= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

    VGV		= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

    GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
        2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
    tauM = Phy.rhof^2 * (2/Sol.delt)^2 + Phy.fmu^2 * GG * 36 + Phy.rhof^2 * abs(VGV);    
    
    tauM = tauM.^(0.5);
    
    tauC = 0.125 * tauM ./(Phy.rhof^2 * (gijDown(:,1) + gijDown(:,2) + gijDown(:,3)));
    
    tauM = 1./tauM;

    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    
    index = 0;
    for i = 1:ElemF.nnt
        for j = 1:ElemF.nnt
            Aij_1 = w(p)*locUX.*locUX.*(DphiDx(:,i).*DphiDx(:,j));
            Aij_2 = w(p)*locUX.*locUY.*(DphiDx(:,i).*DphiDy(:,j));
            Aij_3 = w(p)*locUX.*locUZ.*(DphiDx(:,i).*DphiDz(:,j));
            Aij_4 = w(p)*locUY.*locUX.*(DphiDy(:,i).*DphiDx(:,j));
            Aij_5 = w(p)*locUY.*locUY.*(DphiDy(:,i).*DphiDy(:,j));
            Aij_6 = w(p)*locUY.*locUZ.*(DphiDy(:,i).*DphiDz(:,j));
            Aij_7 = w(p)*locUZ.*locUX.*(DphiDz(:,i).*DphiDx(:,j));
            Aij_8 = w(p)*locUZ.*locUY.*(DphiDz(:,i).*DphiDy(:,j));
            Aij_9 = w(p)*locUZ.*locUZ.*(DphiDz(:,i).*DphiDz(:,j));
            Aij_1 = Aij_1.*volume.*tauM;
            Aij_2 = Aij_2.*volume.*tauM;
            Aij_3 = Aij_3.*volume.*tauM;
            Aij_4 = Aij_4.*volume.*tauM;
            Aij_5 = Aij_5.*volume.*tauM;
            Aij_6 = Aij_6.*volume.*tauM;
            Aij_7 = Aij_7.*volume.*tauM;
            Aij_8 = Aij_8.*volume.*tauM;
            Aij_9 = Aij_9.*volume.*tauM;
            sA1(index+1:index+ElemF.nElem,p) = Aij_1;
            sA2(index+1:index+ElemF.nElem,p) = Aij_2;
            sA3(index+1:index+ElemF.nElem,p) = Aij_3;
            sA4(index+1:index+ElemF.nElem,p) = Aij_4;
            sA5(index+1:index+ElemF.nElem,p) = Aij_5;
            sA6(index+1:index+ElemF.nElem,p) = Aij_6;
            sA7(index+1:index+ElemF.nElem,p) = Aij_7;
            sA8(index+1:index+ElemF.nElem,p) = Aij_8;
            sA9(index+1:index+ElemF.nElem,p) = Aij_9;
            index = index + ElemF.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA6 = sum(sA6,2);
sA7 = sum(sA7,2);
sA8 = sum(sA8,2);
sA9 = sum(sA9,2);

% assemble the matrix
A1 = sparse(iif,jjf,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iif,jjf,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iif,jjf,sA3,Elem.nDof,Elem.nDof);
A4 = sparse(iif,jjf,sA4,Elem.nDof,Elem.nDof);
A5 = sparse(iif,jjf,sA5,Elem.nDof,Elem.nDof);
A6 = sparse(iif,jjf,sA6,Elem.nDof,Elem.nDof);
A7 = sparse(iif,jjf,sA7,Elem.nDof,Elem.nDof);
A8 = sparse(iif,jjf,sA8,Elem.nDof,Elem.nDof);
A9 = sparse(iif,jjf,sA9,Elem.nDof,Elem.nDof);

TauMStabK2 = Phy.rhof^2*[A1+A2+A3+A4+A5+A6+A7+A8+A9 ZeroF ZeroF ZeroF; ZeroF A1+A2+A3+A4+A5+A6+A7+A8+A9 ZeroF ZeroF;ZeroF ZeroF A1+A2+A3+A4+A5+A6+A7+A8+A9 ZeroF; ZeroF ZeroF ZeroF ZeroF];

MatlabFlow = MatlabFlow + TauMStabK2;

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 Dphip

RHS = RHS - TauMStabK2 * u1(:);

sA1 = zeros(nnt^2*ElemF.nElem,nQuad);
sA2 = zeros(nnt^2*ElemF.nElem,nQuad);
sA3 = zeros(nnt^2*ElemF.nElem,nQuad);
sA4 = zeros(nnt^2*ElemF.nElem,nQuad);
sA5 = zeros(nnt^2*ElemF.nElem,nQuad);
sA6 = zeros(nnt^2*ElemF.nElem,nQuad);
sA7 = zeros(nnt^2*ElemF.nElem,nQuad);
sA8 = zeros(nnt^2*ElemF.nElem,nQuad);
sA9 = zeros(nnt^2*ElemF.nElem,nQuad);

for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
    xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
    xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
    xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
    xsInv(:,9) = 1./xsInv(:,9) ;
    xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
    xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
    xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
    xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
    xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
    xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
    xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
    xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
    xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
    c1 = 1.154700538379252E0 ;
    c2 = 5.773502691896259E-01;
    
    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
    
    locUX  = sum(repmat(phi(p,:),ElemF.nElem,1).*ux,2);
    locUY  = sum(repmat(phi(p,:),ElemF.nElem,1).*uy,2);
    locUZ  = sum(repmat(phi(p,:),ElemF.nElem,1).*uz,2);
    
    VG1		= gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
    VG2		= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
    VG3		= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

    VGV		= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

    GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
        2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
    tauM = Phy.rhof^2 * (2/Sol.delt)^2 + Phy.fmu^2 * GG * 36 + Phy.rhof^2 * abs(VGV);
    
    tauM = tauM.^(0.5);
    
    tauC = 0.125 * tauM ./(Phy.rhof^2 * (gijDown(:,1) + gijDown(:,2) + gijDown(:,3)));
    
    tauM = 1./tauM;

    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    
    index = 0;
    for i = 1:ElemF.nnt
        for j = 1:ElemF.nnt
            Aij_1 = w(p)*locUX.*(DphiDx(:,i).*DphiDx(:,j));
            Aij_2 = w(p)*locUX.*(DphiDx(:,i).*DphiDy(:,j));
            Aij_3 = w(p)*locUX.*(DphiDx(:,i).*DphiDz(:,j));
            Aij_4 = w(p)*locUY.*(DphiDy(:,i).*DphiDx(:,j));
            Aij_5 = w(p)*locUY.*(DphiDy(:,i).*DphiDy(:,j));
            Aij_6 = w(p)*locUY.*(DphiDy(:,i).*DphiDz(:,j));
            Aij_7 = w(p)*locUZ.*(DphiDz(:,i).*DphiDx(:,j));
            Aij_8 = w(p)*locUZ.*(DphiDz(:,i).*DphiDy(:,j));
            Aij_9 = w(p)*locUZ.*(DphiDz(:,i).*DphiDz(:,j));
            Aij_1 = Aij_1.*volume.*tauM;
            Aij_2 = Aij_2.*volume.*tauM;
            Aij_3 = Aij_3.*volume.*tauM;
            Aij_4 = Aij_4.*volume.*tauM;
            Aij_5 = Aij_5.*volume.*tauM;
            Aij_6 = Aij_6.*volume.*tauM;
            Aij_7 = Aij_7.*volume.*tauM;
            Aij_8 = Aij_8.*volume.*tauM;
            Aij_9 = Aij_9.*volume.*tauM;
            sA1(index+1:index+ElemF.nElem,p) = Aij_1;
            sA2(index+1:index+ElemF.nElem,p) = Aij_2;
            sA3(index+1:index+ElemF.nElem,p) = Aij_3;
            sA4(index+1:index+ElemF.nElem,p) = Aij_4;
            sA5(index+1:index+ElemF.nElem,p) = Aij_5;
            sA6(index+1:index+ElemF.nElem,p) = Aij_6;
            sA7(index+1:index+ElemF.nElem,p) = Aij_7;
            sA8(index+1:index+ElemF.nElem,p) = Aij_8;
            sA9(index+1:index+ElemF.nElem,p) = Aij_9;
            index = index + ElemF.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA6 = sum(sA6,2);
sA7 = sum(sA7,2);
sA8 = sum(sA8,2);
sA9 = sum(sA9,2);

% assemble the matrix
A1 = sparse(iif,jjf,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iif,jjf,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iif,jjf,sA3,Elem.nDof,Elem.nDof);
A4 = sparse(iif,jjf,sA4,Elem.nDof,Elem.nDof);
A5 = sparse(iif,jjf,sA5,Elem.nDof,Elem.nDof);
A6 = sparse(iif,jjf,sA6,Elem.nDof,Elem.nDof);
A7 = sparse(iif,jjf,sA7,Elem.nDof,Elem.nDof);
A8 = sparse(iif,jjf,sA8,Elem.nDof,Elem.nDof);
A9 = sparse(iif,jjf,sA9,Elem.nDof,Elem.nDof);

FTauMStab = Phy.rhof*[ZeroF ZeroF ZeroF A1+A4+A7 ; ZeroF ZeroF ZeroF A2+A5+A8;ZeroF ZeroF ZeroF A3+A6+A9; ZeroF ZeroF ZeroF ZeroF];

MatlabFlow = MatlabFlow + FTauMStab;

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 Dphip

RHS = RHS - FTauMStab * p1;

sA1 = zeros(nnt^2*ElemF.nElem,nQuad);
sA2 = zeros(nnt^2*ElemF.nElem,nQuad);
sA3 = zeros(nnt^2*ElemF.nElem,nQuad);

for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    
        index = 0;
    for i = 1:ElemF.nnt
        for j = 1:ElemF.nnt
            Aij_1 = w(p)*(DphiDx(:,i).*phi(p,j));
            Aij_2 = w(p)*(DphiDy(:,i).*phi(p,j));
            Aij_3 = w(p)*(DphiDz(:,i).*phi(p,j));

            Aij_1 = Aij_1.*volume;
            Aij_2 = Aij_2.*volume;
            Aij_3 = Aij_3.*volume;
            
            sA1(index+1:index+ElemF.nElem,p) = Aij_1;
            sA2(index+1:index+ElemF.nElem,p) = Aij_2;
            sA3(index+1:index+ElemF.nElem,p) = Aij_3;
            
            index = index + ElemF.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);

% assemble the matrix
A1 = sparse(iif,jjf,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iif,jjf,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iif,jjf,sA3,Elem.nDof,Elem.nDof);

Fp = [ZeroF ZeroF ZeroF A1; ZeroF ZeroF ZeroF A2; ZeroF ZeroF ZeroF A3; ZeroF ZeroF ZeroF ZeroF];

MatlabFlow = MatlabFlow - Fp;

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 Dphip

RHS = RHS+ Fp * p1 ;

GMassCons = Phy.rhof*Fp';

MatlabFlow = MatlabFlow + GMassCons;

RHS = RHS - GMassCons * u1;

sA1 = zeros(nnt^2*ElemF.nElem,nQuad);
sA2 = zeros(nnt^2*ElemF.nElem,nQuad);
sA3 = zeros(nnt^2*ElemF.nElem,nQuad);

for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
    xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
    xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
    xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
    xsInv(:,9) = 1./xsInv(:,9) ;
    xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
    xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
    xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
    xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
    xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
    xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
    xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
    xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
    xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
    c1 = 1.154700538379252E0 ;
    c2 = 5.773502691896259E-01;
    
    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
    
    locUX  = sum(repmat(phi(p,:),ElemF.nElem,1).*ux,2);
    locUY  = sum(repmat(phi(p,:),ElemF.nElem,1).*uy,2);
    locUZ  = sum(repmat(phi(p,:),ElemF.nElem,1).*uz,2);
    
    VG1		= gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
    VG2		= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
    VG3		= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

    VGV		= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

    GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
        2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
    tauM = Phy.rhof^2 * (2/Sol.delt)^2 + Phy.fmu^2 * GG * 36 + Phy.rhof^2 * abs(VGV);
    
    tauM = tauM.^(0.5);
    
    tauC = 0.125 * tauM ./(Phy.rhof^2 * (gijDown(:,1) + gijDown(:,2) + gijDown(:,3)));
    
    tauM = 1./tauM;

    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    
    index = 0;
    for i = 1:ElemF.nnt
        for j = 1:ElemF.nnt
            Aij_1 = w(p)*(DphiDx(:,i).*phi(p,j));
            Aij_2 = w(p)*(DphiDy(:,i).*phi(p,j));
            Aij_3 = w(p)*(DphiDz(:,i).*phi(p,j));
            Aij_1 = Aij_1.*volume.*tauM;
            Aij_2 = Aij_2.*volume.*tauM;
            Aij_3 = Aij_3.*volume.*tauM;
            sA1(index+1:index+ElemF.nElem,p) = Aij_1;
            sA2(index+1:index+ElemF.nElem,p) = Aij_2;
            sA3(index+1:index+ElemF.nElem,p) = Aij_3;
            index = index + ElemF.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);

% assemble the matrix
A1 = sparse(iif,jjf,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iif,jjf,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iif,jjf,sA3,Elem.nDof,Elem.nDof);

GTauMStab1 = 1.5/Sol.delt*Phy.rhof*[ZeroF ZeroF ZeroF ZeroF ; ZeroF ZeroF ZeroF ZeroF;ZeroF ZeroF ZeroF ZeroF; A1 A2 A3 ZeroF];

MatlabFlow = MatlabFlow + GTauMStab1;

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 Dphip

RHS = RHS - GTauMStab1 * u1;
RHS = RHS + 2/1.5 * GTauMStab1 * u2 - 0.5/1.5 * GTauMStab1 * u3;

sA1 = zeros(nnt^2*ElemF.nElem,nQuad);
sA2 = zeros(nnt^2*ElemF.nElem,nQuad);
sA3 = zeros(nnt^2*ElemF.nElem,nQuad);
sA4 = zeros(nnt^2*ElemF.nElem,nQuad);
sA5 = zeros(nnt^2*ElemF.nElem,nQuad);
sA6 = zeros(nnt^2*ElemF.nElem,nQuad);
sA7 = zeros(nnt^2*ElemF.nElem,nQuad);
sA8 = zeros(nnt^2*ElemF.nElem,nQuad);
sA9 = zeros(nnt^2*ElemF.nElem,nQuad);

for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
    xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
    xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
    xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
    xsInv(:,9) = 1./xsInv(:,9) ;
    xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
    xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
    xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
    xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
    xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
    xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
    xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
    xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
    xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
    c1 = 1.154700538379252E0 ;
    c2 = 5.773502691896259E-01;
    
    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
    
    locUX  = sum(repmat(phi(p,:),ElemF.nElem,1).*ux,2);
    locUY  = sum(repmat(phi(p,:),ElemF.nElem,1).*uy,2);
    locUZ  = sum(repmat(phi(p,:),ElemF.nElem,1).*uz,2);
    
    VG1		= gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
    VG2		= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
    VG3		= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

    VGV		= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

    GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
        2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
    tauM = Phy.rhof^2 * (2/Sol.delt)^2 + Phy.fmu^2 * GG * 36 + Phy.rhof^2 * abs(VGV);    
    
    tauM = tauM.^(0.5);
    
    tauC = 0.125 * tauM ./(Phy.rhof^2 * (gijDown(:,1) + gijDown(:,2) + gijDown(:,3)));
    
    tauM = 1./tauM;

    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    
    index = 0;
    for i = 1:ElemF.nnt
        for j = 1:ElemF.nnt
            Aij_1 = w(p)*locUX.*(DphiDx(:,i).*DphiDx(:,j));
            Aij_2 = w(p)*locUY.*(DphiDx(:,i).*DphiDy(:,j));
            Aij_3 = w(p)*locUZ.*(DphiDx(:,i).*DphiDz(:,j));
            Aij_4 = w(p)*locUX.*(DphiDy(:,i).*DphiDx(:,j));
            Aij_5 = w(p)*locUY.*(DphiDy(:,i).*DphiDy(:,j));
            Aij_6 = w(p)*locUZ.*(DphiDy(:,i).*DphiDz(:,j));
            Aij_7 = w(p)*locUX.*(DphiDz(:,i).*DphiDx(:,j));
            Aij_8 = w(p)*locUY.*(DphiDz(:,i).*DphiDy(:,j));
            Aij_9 = w(p)*locUZ.*(DphiDz(:,i).*DphiDz(:,j));
            Aij_1 = Aij_1.*volume.*tauM;
            Aij_2 = Aij_2.*volume.*tauM;
            Aij_3 = Aij_3.*volume.*tauM;
            Aij_4 = Aij_4.*volume.*tauM;
            Aij_5 = Aij_5.*volume.*tauM;
            Aij_6 = Aij_6.*volume.*tauM;
            Aij_7 = Aij_7.*volume.*tauM;
            Aij_8 = Aij_8.*volume.*tauM;
            Aij_9 = Aij_9.*volume.*tauM;
            sA1(index+1:index+ElemF.nElem,p) = Aij_1;
            sA2(index+1:index+ElemF.nElem,p) = Aij_2;
            sA3(index+1:index+ElemF.nElem,p) = Aij_3;
            sA4(index+1:index+ElemF.nElem,p) = Aij_4;
            sA5(index+1:index+ElemF.nElem,p) = Aij_5;
            sA6(index+1:index+ElemF.nElem,p) = Aij_6;
            sA7(index+1:index+ElemF.nElem,p) = Aij_7;
            sA8(index+1:index+ElemF.nElem,p) = Aij_8;
            sA9(index+1:index+ElemF.nElem,p) = Aij_9;
            index = index + ElemF.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA6 = sum(sA6,2);
sA7 = sum(sA7,2);
sA8 = sum(sA8,2);
sA9 = sum(sA9,2);

% assemble the matrix
A1 = sparse(iif,jjf,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iif,jjf,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iif,jjf,sA3,Elem.nDof,Elem.nDof);
A4 = sparse(iif,jjf,sA4,Elem.nDof,Elem.nDof);
A5 = sparse(iif,jjf,sA5,Elem.nDof,Elem.nDof);
A6 = sparse(iif,jjf,sA6,Elem.nDof,Elem.nDof);
A7 = sparse(iif,jjf,sA7,Elem.nDof,Elem.nDof);
A8 = sparse(iif,jjf,sA8,Elem.nDof,Elem.nDof);
A9 = sparse(iif,jjf,sA9,Elem.nDof,Elem.nDof);

GTauMStab1 = Phy.rhof*[ZeroF ZeroF ZeroF ZeroF ; ZeroF ZeroF ZeroF ZeroF;ZeroF ZeroF ZeroF ZeroF; A1+A2+A3 A4+A5+A6 A7+A8+A9 ZeroF];

MatlabFlow = MatlabFlow + GTauMStab1;

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 Dphip

RHS = RHS - GTauMStab1 * u1;

sA1 = zeros(nnt^2*ElemF.nElem,nQuad);
sA2 = zeros(nnt^2*ElemF.nElem,nQuad);
sA3 = zeros(nnt^2*ElemF.nElem,nQuad);

for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxf*[phix(:,p)], xxf*[phiy(:,p)], xxf*[phiz(:,p)], yyf*[phix(:,p)], yyf*[phiy(:,p)], yyf*[phiz(:,p)], zzf*[phix(:,p)], zzf*[phiy(:,p)], zzf*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    xsInv(:,1) = J(:,5).*J(:,9)-J(:,6).*J(:,8) ;
    xsInv(:,2) = J(:,3).*J(:,8)-J(:,2).*J(:,9) ;
    xsInv(:,3) = J(:,2).*J(:,6)-J(:,3).*J(:,5) ;
    xsInv(:,9) = J(:,1).*xsInv(:,1) + J(:,4).*xsInv(:,2) + J(:,7).*xsInv(:,3);
    xsInv(:,9) = 1./xsInv(:,9) ;
    xsInv(:,1) = xsInv(:,9).*xsInv(:,1);
    xsInv(:,2) = xsInv(:,9).*xsInv(:,2);
    xsInv(:,3) = xsInv(:,9).*xsInv(:,3);
    
    xsInv(:,4) = xsInv(:,9).* (-J(:,4).*J(:,9)+J(:,6).*J(:,7));
    xsInv(:,5) = xsInv(:,9).* ( J(:,1).*J(:,9)-J(:,3).*J(:,7));
    xsInv(:,6) = xsInv(:,9).* ( J(:,3).*J(:,4)-J(:,1).*J(:,6));
    
    xsInv(:,7) = xsInv(:,9).* ( J(:,4).*J(:,8)-J(:,7).*J(:,5));
    xsInv(:,8) = xsInv(:,9).* ( J(:,2).*J(:,7)-J(:,1).*J(:,8));
    xsInv(:,9) = xsInv(:,9).* ( J(:,1).*J(:,5)-J(:,2).*J(:,4));
    
    c1 = 1.154700538379252E0 ;
    c2 = 5.773502691896259E-01;
    
    a1 = c1 * xsInv(:,1) + c2 * xsInv(:,4);
    a2 = c1 * xsInv(:,4) + c2 * xsInv(:,1);
    
    gijDown(:,1) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,7) ;
    
    a1 = c1 * xsInv(:,2) + c2 * xsInv(:,5) ;
    a2 = c1 * xsInv(:,5) + c2 * xsInv(:,2) ;
    gijDown(:,2) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,8);
    gijDown(:,4) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,8);
    
    a1 = c1 * xsInv(:,3) + c2 * xsInv(:,6) ;
    a2 = c1 * xsInv(:,6) + c2 * xsInv(:,3) ;
    gijDown(:,3) = xsInv(:,3) .* a1 + xsInv(:,6) .* a2 + xsInv(:,9) .* xsInv(:,9);
    gijDown(:,5) = xsInv(:,2) .* a1 + xsInv(:,5) .* a2 + xsInv(:,8) .* xsInv(:,9);
    gijDown(:,6) = xsInv(:,1) .* a1 + xsInv(:,4) .* a2 + xsInv(:,7) .* xsInv(:,9);
    
    locUX  = sum(repmat(phi(p,:),ElemF.nElem,1).*ux,2);
    locUY  = sum(repmat(phi(p,:),ElemF.nElem,1).*uy,2);
    locUZ  = sum(repmat(phi(p,:),ElemF.nElem,1).*uz,2);
    
    VG1		= gijDown(:,1) .* locUX + gijDown(:,4) .* locUY + gijDown(:,6) .* locUZ;
    VG2		= gijDown(:,4) .* locUX + gijDown(:,2) .* locUY + gijDown(:,5) .* locUZ;
    VG3		= gijDown(:,6) .* locUX	+ gijDown(:,5) .* locUY	+ gijDown(:,3) .* locUZ;

    VGV		= VG1 .* locUX + VG2 .* locUY + VG3 .* locUZ;

    GG = gijDown(:,1).^2 + gijDown(:,2).^2 + gijDown(:,3).^2 + ...
        2 * ( gijDown(:,4).^2 + gijDown(:,5).^2 + gijDown(:,6).^2 );
    
    tauM = Phy.rhof^2 * (2/Sol.delt)^2 + Phy.fmu^2 * GG * 36 + Phy.rhof^2 * abs(VGV);
    
    tauM = tauM.^(0.5);
    
    tauC = 0.125 * tauM ./(Phy.rhof^2 * (gijDown(:,1) + gijDown(:,2) + gijDown(:,3)));
    
    tauM = 1./tauM;

    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemF.nnt);
    
    index = 0;
    for i = 1:ElemF.nnt
        for j = 1:ElemF.nnt
            Aij_1 = w(p)*(DphiDx(:,i).*DphiDx(:,j));
            Aij_2 = w(p)*(DphiDy(:,i).*DphiDy(:,j));
            Aij_3 = w(p)*(DphiDz(:,i).*DphiDz(:,j));
            Aij_1 = Aij_1.*volume.*tauM;
            Aij_2 = Aij_2.*volume.*tauM;
            Aij_3 = Aij_3.*volume.*tauM;
            sA1(index+1:index+ElemF.nElem,p) = Aij_1;
            sA2(index+1:index+ElemF.nElem,p) = Aij_2;
            sA3(index+1:index+ElemF.nElem,p) = Aij_3;
            index = index + ElemF.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);

% assemble the matrix
A1 = sparse(iif,jjf,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iif,jjf,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iif,jjf,sA3,Elem.nDof,Elem.nDof);

CTauMStab = [ZeroF ZeroF ZeroF ZeroF ; ZeroF ZeroF ZeroF ZeroF;ZeroF ZeroF ZeroF ZeroF; ZeroF ZeroF ZeroF A1+A2+A3];

MatlabFlow = MatlabFlow + CTauMStab;

RHS = RHS - CTauMStab * p1;

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 A7 A8 A9 Dphip

%% SOLID

xxs = zeros(size(ElemS.elem2dof));
yys = zeros(size(ElemS.elem2dof));
zzs = zeros(size(ElemS.elem2dof));

for i=1:ElemS.nnt
    xxs(:,i) = Elem.dof(ElemS.elem2dof(:,i),1);
    yys(:,i) = Elem.dof(ElemS.elem2dof(:,i),2);
    zzs(:,i) = Elem.dof(ElemS.elem2dof(:,i),3);
end

ElemS.nElem = size(ElemS.elem2dof,1);
        
% compute non-zeros
sA = [];
for p = 1:nQuad    
    J = [xxs*[phix(:,p)], xxs*[phiy(:,p)], xxs*[phiz(:,p)], yys*[phix(:,p)], yys*[phiy(:,p)], yys*[phiz(:,p)], zzs*[phix(:,p)], zzs*[phiy(:,p)], zzs*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    negJacobian = find(volume<0);
    if ~isempty(negJacobian)
        disp('Mesh deformed, Negative Jacobian');
        exit
    end
    
    index = 0;
    for i = 1:nnt
        for j = 1:nnt
            Mij = w(p)*(phi(p,i)'*phi(p,j));
            Mij = Mij.*volume;
            sA(index+1:index+ElemS.nElem,p) = Mij;
            index = index + ElemS.nElem;
        end
    end
end
sA = sum(sA,2);
% assemble the matrix
Ms = sparse(iis,jjs,sA,Elem.nDof,Elem.nDof);
Ms = Ms';
ZeroF = sparse(Elem.nDof,Elem.nDof);
Ms = [Ms ZeroF ZeroF ZeroF;ZeroF Ms ZeroF ZeroF;ZeroF ZeroF Ms ZeroF; ZeroF ZeroF ZeroF ZeroF];
Msolid = 1.5 * Phy.rhos * Ms/Sol.delt;
MatlabFlow = MatlabFlow + Msolid;
clear Mij sA

RHS = RHS -1.5*Phy.rhos/Sol.delt * Ms * u1(:) + 2 * Phy.rhos/Sol.delt * Ms * u2(:) - 0.5 * Phy.rhos/Sol.delt * Ms * u3(:);

sA1 = zeros(nnt^2*ElemS.nElem,nQuad);
sA2 = zeros(nnt^2*ElemS.nElem,nQuad);
sA3 = zeros(nnt^2*ElemS.nElem,nQuad);
sA4 = zeros(nnt^2*ElemS.nElem,nQuad);
sA5 = zeros(nnt^2*ElemS.nElem,nQuad);
sA6 = zeros(nnt^2*ElemS.nElem,nQuad);

for p = 1:nQuad
    % Dphi at quadrature points
    J = [xxs*[phix(:,p)], xxs*[phiy(:,p)], xxs*[phiz(:,p)], yys*[phix(:,p)], yys*[phiy(:,p)], yys*[phiz(:,p)], zzs*[phix(:,p)], zzs*[phiy(:,p)], zzs*[phiz(:,p)]];
    if size(J,2)==1
        J = J';
    end
    
    volume =( J(:,1).*(J(:,5).*J(:,9)-J(:,6).*J(:,8)) + ...
                    J(:,4).*(J(:,8).*J(:,3)-J(:,2).*J(:,9)) + ...
                    J(:,7).*(J(:,2).*J(:,6)-J(:,3).*J(:,5)) );
                
    volume = abs(volume);
    
    DphiDx = ((J(:,5).*J(:,9)-J(:,8).*J(:,6))*phix(:,p)'+ ...
        (J(:,6).*J(:,7)-J(:,4).*J(:,9))*phiy(:,p)'+ ...
        (J(:,4).*J(:,8)-J(:,5).*J(:,7))*phiz(:,p)')./repmat(volume,1,ElemS.nnt);
    DphiDy = ((J(:,3).*J(:,8)-J(:,2).*J(:,9))*phix(:,p)'+ ...
        (J(:,1).*J(:,9)-J(:,7).*J(:,3))*phiy(:,p)'+ ...
        (J(:,2).*J(:,7)-J(:,1).*J(:,8))*phiz(:,p)')./repmat(volume,1,ElemS.nnt);
    DphiDz = ((J(:,2).*J(:,6)-J(:,3).*J(:,5))*phix(:,p)'+ ...
        (J(:,3).*J(:,4)-J(:,1).*J(:,6))*phiy(:,p)'+ ...
        (J(:,1).*J(:,5)-J(:,4).*J(:,2))*phiz(:,p)')./repmat(volume,1,ElemS.nnt);
    
        index = 0;
    for i = 1:ElemS.nnt
        for j = 1:ElemS.nnt
            Aij_1 = w(p)*(DphiDx(:,i).*DphiDx(:,j));
            Aij_2 = w(p)*(DphiDy(:,i).*DphiDy(:,j));
            Aij_3 = w(p)*(DphiDz(:,i).*DphiDz(:,j));
            Aij_4 = w(p)*(DphiDx(:,i).*DphiDy(:,j));
            Aij_5 = w(p)*(DphiDy(:,i).*DphiDz(:,j));
            Aij_6 = w(p)*(DphiDz(:,i).*DphiDx(:,j));
            Aij_1 = Aij_1.*volume;
            Aij_2 = Aij_2.*volume;
            Aij_3 = Aij_3.*volume;
            Aij_4 = Aij_4.*volume;
            Aij_5 = Aij_5.*volume;
            Aij_6 = Aij_6.*volume;
            sA1(index+1:index+ElemS.nElem,p) = Aij_1;
            sA2(index+1:index+ElemS.nElem,p) = Aij_2;
            sA3(index+1:index+ElemS.nElem,p) = Aij_3;
            sA4(index+1:index+ElemS.nElem,p) = Aij_4;
            sA5(index+1:index+ElemS.nElem,p) = Aij_5;
            sA6(index+1:index+ElemS.nElem,p) = Aij_6;
            index = index + ElemS.nElem;
        end
    end
end
sA1 = sum(sA1,2);
sA2 = sum(sA2,2);
sA3 = sum(sA3,2);
sA4 = sum(sA4,2);
sA5 = sum(sA5,2);
sA6 = sum(sA6,2);

% assemble the matrix
A1 = sparse(iis,jjs,sA1,Elem.nDof,Elem.nDof);
A2 = sparse(iis,jjs,sA2,Elem.nDof,Elem.nDof);
A3 = sparse(iis,jjs,sA3,Elem.nDof,Elem.nDof);
A4 = sparse(iis,jjs,sA4,Elem.nDof,Elem.nDof);
A5 = sparse(iis,jjs,sA5,Elem.nDof,Elem.nDof);
A6 = sparse(iis,jjs,sA6,Elem.nDof,Elem.nDof);

Ks = Phy.smu*[2*A1+A2+A3 A4' A6 ZeroF; A4 A1+2*A2+A3 A5' ZeroF; A6' A5 A1+A2+2*A3 ZeroF; ZeroF ZeroF ZeroF ZeroF]+Phy.lambda*[A1 A4 A6' ZeroF; A4' A2 A5 ZeroF; A6 A5' A3 ZeroF; ZeroF ZeroF ZeroF ZeroF];

MatlabFlow = MatlabFlow + 0.75 * Sol.delt * Ks;

disp1 = disp(:,:,1);
disp2 = disp(:,:,2);

disp1 = [disp1(:);zeros(Elem.nDof,1)];
disp2 = [disp2(:);zeros(Elem.nDof,1)];

RHS = RHS - 0.75 * Sol.delt * Ks * u1(:) + 0.25 * Sol.delt * Ks * u2(:) - 0.5 * Ks * disp1 - 0.5 * Ks * disp2 ;

clear sA1 sA2 sA3 aA4 sA5 sA6 Aij_1 Aij_2 Aij_3 Aij_4 Aij_5 Aij_6 
clear A1 A2 A3 A4 A5 A6 Dphip


%% LINEAR SOLVER

if flag
    freeNodes = unique([inletNodes;topNodes;bottomNodes;hinged_1;side1Nodes;side2Nodes ]);
else
    freeNodes = find(Elem.dof(:,2)==1 | Elem.dof(:,2)== -1 | (Elem.dof(:,1)==1 & Elem.dof(:,2)<=0) | Elem.dof(:,1)==0);
end

freeNodes = setdiff(1:size(Elem.dof,1),freeNodes);
freeNodes = [freeNodes';freeNodes' + size(Elem.dof,1)];
freeNodes = [freeNodes; [(3*size(Elem.dof,1)+1):(4*size(Elem.dof,1))]'];


if flag
    intSolidNode = 3*size(Elem.dof,1) + setdiff(solidNodes,itfNodes) ;
else
    intSolidNode = 3*size(Elem.dof,1) + find(Elem.dof(:,2)<0);
end

freeNodes = setdiff(freeNodes,intSolidNode);


result = Sol.u(:,:,1);
result = result(:);
result = [result;Sol.p];

linSol = MatlabFlow(freeNodes,freeNodes)\RHS(freeNodes);
result(freeNodes) = result(freeNodes) + linSol;

Sol.u(:,:,1) = reshape(result(1:3*Elem.nDof),[],3);
Sol.p = result((3*Elem.nDof+1):(4*Elem.nDof));

unp1 = Sol.u(:,:,1);
un   = Sol.u(:,:,2);
unm1 = Sol.u(:,:,3);

Sol.u(:,:,4) = Sol.u(:,:,3);
Sol.u(:,:,3) = Sol.u(:,:,2);
Sol.u(:,:,2) = Sol.u(:,:,1);

test = [Elem.dof unp1 disp(:,:,1)];

disp(:,:,3) = disp(:,:,2);
disp(:,:,2) = disp(:,:,1);


if mod(timeStep,100)==0

data = [Elem.dof+disp(:,:,1) Sol.u(:,:,1)  Sol.p disp(:,:,1)];
newElem = [ElemF.elem2dof(:,1:3) ElemF.elem2dof(:,3) ElemF.elem2dof(:,4:6) ElemF.elem2dof(:,6)];

outFileName = ['Filament' int2str(timeStep) '.plt'];

fileID = fopen(outFileName,'w');

FileTitle = 'velocity plot';

fprintf(fileID,' TITLE = \"%s\"\n',FileTitle);

if strcmp(FileTitle,'velocity plot')
        fprintf(fileID,' VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"p\", \"DX\", \"DY\", \"DZ\"\n');
end
fprintf(fileID,' ZONE SOLUTIONTIME=%f, NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK\n',(timeStep-1)*Sol.delt,Elem.nDof,ElemF.nElem);
fprintf(fileID,'%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n',data');

fprintf(fileID,'%d %d %d %d %d %d %d %d\n',newElem');

fclose(fileID);
end

tipDisp = [tipDisp;disp(nodeId,:,1)];

if mod(timeStep,1)==0
    save('tipDisp1.mat','tipDisp');
end
end

