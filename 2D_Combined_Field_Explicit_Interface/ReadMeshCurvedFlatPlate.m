function ReadMeshCurvedFlatPlate()

clear all;
fileID = fopen('CurDefPlateEle.txt');
C = textscan(fileID, '%d  %d  %d  %d  %d  %d');
fclose(fileID);
% A1 = double(cell2mat(C(1,1)));
% A2 = double(cell2mat(C(1,2)));
% A3 = double(cell2mat(C(1,3)));
A4 = double(cell2mat(C(1,4)));
A5 = double(cell2mat(C(1,5)));
A6 = double(cell2mat(C(1,6)));
ElemF = [A4 A5 A6];

fileID = fopen('CurDefPlateNod.txt');
C = textscan(fileID, '%d %f %f');
fclose(fileID);
% B1 = double(cell2mat(C(1,1)));
B2 = double(cell2mat(C(1,2)));
B3 = double(cell2mat(C(1,3)));
NodeF = [B2 B3];


% this one of the rows which is required to be modified for the data
% corrections
Nodef = NodeF-repmat([0,4.983],size(NodeF,1),1);
Nodefc = Nodef(:,1)+sqrt(-1)*Nodef(:,2);
[r,ar,ad] = cart2polar(Nodefc);
Nodef = [r ad];

InfNode1 = find(abs(Nodef(:,1)-4.978)<0.001 & Nodef(:,2)>=(-90.001) & Nodef(:,2)<=(-78.499));
InfNode2 = find(abs(Nodef(:,1)-4.988)<0.001 & Nodef(:,2)>=(-90.001) & Nodef(:,2)<=(-78.499));

InfNode = sortrows(round(Nodef([InfNode1;InfNode2],:)*10^3)*10^(-3),[1 2]);

% InfNode = sortrows([InfNode;InfNode(1,1),0;InfNode(size(InfNode,1),1),0],[2 1]);
% innerNode_solid = [InfNode(:,1),zeros(size(InfNode(:,1)))];
innerNode_solid = [];
Nodes = [InfNode;innerNode_solid];
[NodeS(:,1),NodeS(:,2)] = polar2cart(Nodes(:,1),Nodes(:,2));

NodeS = NodeS + repmat([0,4.983],size(NodeS,1),1);

% this line gives the number of elements present in each row the solid.
nt = size(find(abs(Nodes(:,1)-4.978)<=1e-5),1)-1;
% this is total number of nodes present.
nnodes = size(Nodes,1);
% this step defines the number of rows of elements present in each solid
nrows = nnodes/(nt+1)-1;
% this step generates the element numbers in a given row.
ElemNo = (1:nt)';
% this step generates the elements numbers in a given column and need to
% add this step over here.

ElemS = [ElemNo,ElemNo+nt+1,ElemNo+nt+2;ElemNo,ElemNo+1,ElemNo+nt+2];

bd1NodeS = find(abs(Nodes(:,1)-4.978)<=0.001 & Nodes(:,2)>=(-90) & Nodes(:,2)<=(-78.499));
bd2NodeS = find(Nodes(:,2)<=(-78.499) & Nodes(:,2)>=(-78.501) & (Nodes(:,1)>=4.977) & (Nodes(:,1)<=4.989));
bd3NodeS = find( abs(Nodes(:,1)-4.988)<=0.001 & Nodes(:,2)>=(-90) & Nodes(:,2)<=(-78.499));
bd4NodeS = find(Nodes(:,2)<=(-89.999) & Nodes(:,2)>=(-90.001) & (Nodes(:,1)>=4.977)>=0.001 & (Nodes(:,1)<=4.989));

[i,bd1NodeSIndex] = sortrows(Nodes(bd1NodeS,:),2);
[i,bd2NodeSIndex] = sortrows(Nodes(bd2NodeS,:),1);
[i,bd3NodeSIndex] = sortrows(Nodes(bd3NodeS,:),2);
[i,bd4NodeSIndex] = sortrows(Nodes(bd4NodeS,:),1);

gamSbd1 = bd1NodeS(bd1NodeSIndex);
gamSbd2 = bd2NodeS(bd2NodeSIndex);
gamSbd3 = bd3NodeS(bd3NodeSIndex);
gamSbd4 = bd4NodeS(bd4NodeSIndex);

gamSbd1 = [gamSbd1(1:length(gamSbd1)-1) gamSbd1(2:length(gamSbd1))];
gamSbd2 = [gamSbd2(1:length(gamSbd2)-1) gamSbd2(2:length(gamSbd2))];
gamSbd3 = [gamSbd3(1:length(gamSbd3)-1) gamSbd3(2:length(gamSbd3))];
gamSbd4 = [gamSbd4(1:length(gamSbd4)-1) gamSbd4(2:length(gamSbd4))];

n1 = size(gamSbd1,1);
n2 = size(gamSbd2,1);
n3 = size(gamSbd3,1);
n4 = size(gamSbd4,1);

gamS = [gamSbd1;gamSbd2;gamSbd3;gamSbd4];
gamIS = [1 n1;n1+1 n2+n1;n1+n2+1 n1+n2+n3;n1+n2+n3+1 n1+n2+n3+n4];

% this step need to modified..
iS = [1 2 3];
cS = [1 3];

dS = [1 2 3 4];
nS = [];

DS = gamIS(dS',:);
NS = gamIS(nS',:);

bd1NodeF = find(abs(Nodef(:,1)-4.978)<=0.001 & Nodef(:,2)>=(-90) & Nodef(:,2)<=(-78.499));
bd2NodeF = find(Nodef(:,2)<=(-78.499) & Nodef(:,2)>=(-78.501) & (Nodef(:,1)>=4.977) & (Nodef(:,1)<=4.989));
bd3NodeF = find( abs(Nodef(:,1)-4.988)<=0.001 & Nodef(:,2)>=(-90) & Nodef(:,2)<=(-78.499));
bd4NodeF = find(0.005>=NodeF(:,2) & NodeF(:,2)>=(-0.005) & NodeF(:,1)== (-2));
bd5NodeF = find( 10>=NodeF(:,1) & NodeF(:,1)>=(-2) & NodeF(:,2)==1);
bd6NodeF = find(1>=NodeF(:,2) & NodeF(:,2)>=(-1) & NodeF(:,1)== 10);
bd7NodeF = find(10>=NodeF(:,1) & NodeF(:,1)>=(-2) & NodeF(:,2)== -1);
bd8NodeF = find( 1>=NodeF(:,2) & NodeF(:,2)>=(-1) & NodeF(:,1)== (-2));
bd9NodeF = find(0>=NodeF(:,1) & NodeF(:,1)>=(-2) & NodeF(:,2)== -0.005);
bd10NodeF = find(0>=NodeF(:,1) & NodeF(:,1)>=(-2) & NodeF(:,2)== 0.005);

[i,bd1NodeFIndex] = sortrows(Nodef(bd1NodeF,:),2);
[i,bd2NodeFIndex] = sortrows(Nodef(bd2NodeF,:),1);
[i,bd3NodeFIndex] = sortrows(Nodef(bd3NodeF,:),2);
[i,bd4NodeFIndex] = sortrows(Nodef(bd4NodeF,:),1);
[i,bd5NodeFIndex] = sortrows(NodeF(bd5NodeF,:),1);
[i,bd6NodeFIndex] = sortrows(NodeF(bd6NodeF,:),2);
[i,bd7NodeFIndex] = sortrows(NodeF(bd7NodeF,:),1);
[i,bd8NodeFIndex] = sortrows(NodeF(bd8NodeF,:),2);
[i,bd9NodeFIndex] = sortrows(NodeF(bd9NodeF,:),1);
[i,bd10NodeFIndex] = sortrows(NodeF(bd10NodeF,:),1);

gamFbd1 = bd1NodeF(bd1NodeFIndex);
gamFbd2 = bd2NodeF(bd2NodeFIndex);
gamFbd3 = bd3NodeF(bd3NodeFIndex);
gamFbd4 = bd4NodeF(bd4NodeFIndex);
gamFbd5 = bd5NodeF(bd5NodeFIndex);
gamFbd6 = bd6NodeF(bd6NodeFIndex);
gamFbd7 = bd7NodeF(bd7NodeFIndex);
gamFbd8 = bd8NodeF(bd8NodeFIndex);
gamFbd9 = bd9NodeF(bd9NodeFIndex);
gamFbd10 = bd10NodeF(bd10NodeFIndex);

gamFbd1 = [gamFbd1(1:length(gamFbd1)-1) gamFbd1(2:length(gamFbd1))];
gamFbd2 = [gamFbd2(1:length(gamFbd2)-1) gamFbd2(2:length(gamFbd2))];
gamFbd3 = [gamFbd3(1:length(gamFbd3)-1) gamFbd3(2:length(gamFbd3))];
gamFbd4 = [gamFbd4(1:length(gamFbd4)-1) gamFbd4(2:length(gamFbd4))];
gamFbd5 = [gamFbd5(1:length(gamFbd5)-1) gamFbd5(2:length(gamFbd5))];
gamFbd6 = [gamFbd6(1:length(gamFbd6)-1) gamFbd6(2:length(gamFbd6))];
gamFbd7 = [gamFbd7(1:length(gamFbd7)-1) gamFbd7(2:length(gamFbd7))];
gamFbd8 = [gamFbd8(1:length(gamFbd8)-1) gamFbd8(2:length(gamFbd8))];
gamFbd9 = [gamFbd9(1:length(gamFbd9)-1) gamFbd9(2:length(gamFbd9))];
gamFbd10 = [gamFbd10(1:length(gamFbd10)-1) gamFbd10(2:length(gamFbd10))];

gamFbd8 = setdiff(sort(gamFbd8,2),sort(gamFbd4,2),'rows');

n1 = size(gamFbd1,1);
n2 = size(gamFbd2,1);
n3 = size(gamFbd3,1);
% n4 = size(gamFbd4,1); % comment this line
n5 = size(gamFbd5,1);
n6 = size(gamFbd6,1);
n7 = size(gamFbd7,1);
n8 = size(gamFbd8,1);
n9 = size(gamFbd9,1);
n10 = size(gamFbd10,1);

% gamF = [gamFbd1;gamFbd2;gamFbd3;gamFbd4;gamFbd5;gamFbd6;gamFbd7;gamFbd8];
% gamIF = [1 n1;n1+1 n2+n1;n1+n2+1 n1+n2+n3;n1+n2+n3+1 n1+n2+n3+n4;n1+n2+n3+n4+1 n1+n2+n3+n4+n5;...
%     n1+n2+n3+n4+n5+1 n1+n2+n3+n4+n5+n6; n1+n2+n3+n4+n5+n6+1 n1+n2+n3+n4+n5+n6+n7; n1+n2+n3+n4+n5+n6+n7+1 n1+n2+n3+n4+n5+n6+n7+n8];

gamF = [gamFbd1;gamFbd2;gamFbd3;gamFbd5;gamFbd6;gamFbd7;gamFbd8;gamFbd9;gamFbd10];
gamIF = [1 n1;n1+1 n2+n1;n1+n2+1 n1+n2+n3;n1+n2+n3+1 n1+n2+n3+n5;...
    n1+n2+n3+n5+1 n1+n2+n3+n5+n6; n1+n2+n3+n5+n6+1 n1+n2+n3+n5+n6+n7;n1+n2+n3+n5+n6+n7+1 n1+n2+n3+n5+n6+n7+n8;...
    n1+n2+n3+n5+n6+n7+n8+1 n1+n2+n3+n5+n6+n7+n8+n9;n1+n2+n3+n5+n6+n7+n8+n9+1 n1+n2+n3+n5+n6+n7+n8+n9+n10];

% carefull while entering the values because we have a data structure where
% Neummann nodes are entered first and Dirchlet nodes are entered second.
% iF = [1 2 3 5];
% cF = [];
% 
% dF = [4 5 7 8];
% nF = [1 2 3 6];

iF = [1 2 3];
cF = [1 3];

dF = [4 6 7 8 9];
nF = [1 2 3 5];

DF = gamIF(dF',:);
NF = gamIF(nF',:);

save('./FlexPlate/DataF.mat','DF','NF','cF','iF','gamIF','gamF','ElemF','NodeF','InfNode');
save('./FlexPlate/DataS.mat','DS','NS','cS','iS','gamIS','gamS','ElemS','NodeS','InfNode');


function [x,y]= polar2cart (mag, ang_in_deg)
 x = mag .* cos(ang_in_deg*pi/180);
 y = mag .* sin(ang_in_deg*pi/180);
  
 
function [r, ar, ad] = cart2polar(x)
 r = abs(x);
 ar = angle(x);
 ad = ar*180/pi;

