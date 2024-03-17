% MATLAB codes for Finite Element Analysis
% problem9.m
% antonio ferreira 2008
% clear memory
function problem9()
clear all
% E; modulus of elasticity
% I: second moment of area
% L: length of bar
E=1; I=1; EI=E*I;
% generation of coordinates and connectivities
numberElements=80;
nodeCoordinates=linspace(0,1,numberElements+1)';
xx=nodeCoordinates;L=max(nodeCoordinates);
numberNodes=size(nodeCoordinates,1);
xx=nodeCoordinates(:,1);
for i=1:numberElements;
elementNodes(i,1)=i;
elementNodes(i,2)=i+1;
end
% distributed load
P=-1;
% for structure:
% displacements: displacement vector
% force : force vector
% stiffness: stiffness matrix
% GDof: global number of degrees of freedom
GDof=2*numberNodes;
U=zeros(GDof,1);
% stiffess matrix and force vector
[stiffness,force]=...
formStiffnessBernoulliBeam(GDof,numberElements,...
elementNodes,numberNodes,xx,EI,P);
% boundary conditions and solution
% clamped-clamped
%fixedNodeU =[1 2*numberElements+1]’;
%fixedNodeV =[2 2*numberElements+2]’;
% simply supported-simply supported
fixedNodeU =[1 2*numberElements+1]'; fixedNodeV =[]';
% clamped at x=0
%fixedNodeU =[1]’; fixedNodeV =[2]’;
prescribedDof=[fixedNodeU;fixedNodeV];
% solution
displacements=solution(GDof,prescribedDof,stiffness,force);
% output displacements/reactions
outputDisplacementsReactions(displacements,stiffness,...
GDof,prescribedDof)
% drawing deformed shape
U=displacements(1:2:2*numberNodes);
plot(nodeCoordinates,U,'.')


%This code calls function formStiffnessBernoulliBeam.m for the computation of the
%stiffness matrix and the force vector for the Bernoulli beam two-node element.

function [stiffness,force] = formStiffnessBernoulliBeam(GDof,numberElements,elementNodes,numberNodes,xx,EI,P)
force=zeros(GDof,1);
stiffness=zeros(GDof);
% calculation of the system stiffness matrix
% and force vector
for e=1:numberElements;
% elementDof: element degrees of freedom (Dof)
indice=elementNodes(e,:) ;
elementDof=[ 2*(indice(1)-1)+1 2*(indice(2)-1)...
2*(indice(2)-1)+1 2*(indice(2)-1)+2];
% ndof=length(indice);
% length of element
LElem=xx(indice(2))-xx(indice(1)) ;
ll=LElem;
k1=EI/(LElem)^3*[12 6*LElem -12 6*LElem;
6*LElem 4*LElem^2 -6*LElem 2*LElem^2;
-12 -6*LElem 12 -6*LElem ;
6*LElem 2*LElem^2 -6*LElem 4*LElem^2];
f1=[P*LElem/2 P*LElem*LElem/12 P*LElem/2 ...
-P*LElem*LElem/12]';
% equivalent force vector
force(elementDof)=force(elementDof)+f1;
% stiffness matrix
stiffness(elementDof,elementDof)=stiffness(elementDof,elementDof)+k1;
end