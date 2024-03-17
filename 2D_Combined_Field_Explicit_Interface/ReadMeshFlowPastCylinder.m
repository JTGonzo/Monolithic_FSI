% this function is used for the flow over a cylinder which is a bench mark
% problem for laminar flow

function ReadMeshFlowPastCylinder()

clear all;
clc;
fileID = fopen('BenchMarkFElement.txt');
C = textscan(fileID, '%d  %d  %d  %d  %d  %d');
fclose(fileID);
% A1 = double(cell2mat(C(1,1)));
% A2 = double(cell2mat(C(1,2)));
% A3 = double(cell2mat(C(1,3)));
A4 = double(cell2mat(C(1,4)));
A5 = double(cell2mat(C(1,5)));
A6 = double(cell2mat(C(1,6)));
ElemF = [A4 A5 A6];

fileID = fopen('BenchMarkFNode.txt');
C = textscan(fileID, '%d %f %f');
fclose(fileID);
% B1 = double(cell2mat(C(1,1)));
B2 = double(cell2mat(C(1,2)));
B3 = double(cell2mat(C(1,3)));
NodeF = [B2 B3];

fileID = fopen('BenchMarkSElement.txt');
C = textscan(fileID, '%d  %d  %d  %d  %d  %d');
fclose(fileID);
% A1 = double(cell2mat(C(1,1)));
% A2 = double(cell2mat(C(1,2)));
% A3 = double(cell2mat(C(1,3)));
A4 = double(cell2mat(C(1,4)));
A5 = double(cell2mat(C(1,5)));
A6 = double(cell2mat(C(1,6)));
ElemS = [A4 A5 A6];

fileID = fopen('BenchMarkSNode.txt');
C = textscan(fileID, '%d %f %f');
fclose(fileID);
% B1 = double(cell2mat(C(1,1)));
B2 = double(cell2mat(C(1,2)));
B3 = double(cell2mat(C(1,3)));
NodeS = [B2 B3];

% give the centre of the cicle in the values below
Nodef = NodeF-repmat([0.2,0.2],size(NodeF,1),1);
[r,ad]=cart2polar(Nodef);
Nodef = [r ad];
Nodes = NodeS-repmat([0.2,0.2],size(NodeS,1),1);
[r,ad]=cart2polar(Nodes);
Nodes = [r ad];

% use the radius in the codes below
bd1NodeS = find(abs(Nodes(:,1)-0.05)<=0.001);% & Nodes(:,2)>=(-90) & Nodes(:,2)<=(-78.499));
% bd2NodeS = find(Nodes(:,2)<=(-78.499) & Nodes(:,2)>=(-78.501) & (Nodes(:,1)>=4.977) & (Nodes(:,1)<=4.989));
% bd3NodeS = find( abs(Nodes(:,1)-4.988)<=0.001 & Nodes(:,2)>=(-90) & Nodes(:,2)<=(-78.499));
% bd4NodeS = find(Nodes(:,2)<=(-89.999) & Nodes(:,2)>=(-90.001) & (Nodes(:,1)>=4.977)>=0.001 & (Nodes(:,1)<=4.989));

[i,bd1NodeSIndex] = sortrows(Nodes(bd1NodeS,:),2);
% [i,bd2NodeSIndex] = sortrows(Nodes(bd2NodeS,:),1);
% [i,bd3NodeSIndex] = sortrows(Nodes(bd3NodeS,:),2);
% [i,bd4NodeSIndex] = sortrows(Nodes(bd4NodeS,:),1);

gamSbd1 = bd1NodeS(bd1NodeSIndex);
% gamSbd2 = bd2NodeS(bd2NodeSIndex);
% gamSbd3 = bd3NodeS(bd3NodeSIndex);
% gamSbd4 = bd4NodeS(bd4NodeSIndex);

gamSbd1 = [gamSbd1(1:length(gamSbd1)-1) gamSbd1(2:length(gamSbd1)); gamSbd1(1) gamSbd1(length(gamSbd1))];
% gamSbd2 = [gamSbd2(1:length(gamSbd2)-1) gamSbd2(2:length(gamSbd2))];
% gamSbd3 = [gamSbd3(1:length(gamSbd3)-1) gamSbd3(2:length(gamSbd3))];
% gamSbd4 = [gamSbd4(1:length(gamSbd4)-1) gamSbd4(2:length(gamSbd4))];

n1 = size(gamSbd1,1);
% n2 = size(gamSbd2,1);
% n3 = size(gamSbd3,1);
% n4 = size(gamSbd4,1);

% gamS = [gamSbd1;gamSbd2;gamSbd3;gamSbd4];
% gamIS = [1 n1;n1+1 n2+n1;n1+n2+1 n1+n2+n3;n1+n2+n3+1 n1+n2+n3+n4];

gamS = [gamSbd1];
gamIS = [1 n1];

% this step need to modified..
iS = [1];
cS = [1];

% dS = [1 2 3 4];
dS = [1];
nS = [];

DS = gamIS(dS',:);
NS = gamIS(nS',:);

% use the radius in the lines below
bd1NodeF = find(abs(Nodef(:,1)-0.05)<=0.001);% & Nodef(:,2)>=(-90) & Nodef(:,2)<=(-78.499));
% bd2NodeF = find(Nodef(:,2)<=(-78.499) & Nodef(:,2)>=(-78.501) & (Nodef(:,1)>=4.977) & (Nodef(:,1)<=4.989));
% bd3NodeF = find( abs(Nodef(:,1)-4.988)<=0.001 & Nodef(:,2)>=(-90) & Nodef(:,2)<=(-78.499));
% bd4NodeF = find(Nodef(:,2)<=(-89.999) & Nodef(:,2)>=(-90.001) & (Nodef(:,1)>=4.977)>=0.001 & (Nodef(:,1)<=4.989));
bd5NodeF = find( 2.2>=NodeF(:,1) & NodeF(:,1)>=0 & NodeF(:,2)==0.41);
bd6NodeF = find(0.41>=NodeF(:,2) & NodeF(:,2)>=0 & NodeF(:,1)== 2.2);
bd7NodeF = find(2.2>=NodeF(:,1) & NodeF(:,1)>=0 & NodeF(:,2)== 0);
bd8NodeF = find( 0.41>=NodeF(:,2) & NodeF(:,2)>=0 & NodeF(:,1)== 0);

[i,bd1NodeFIndex] = sortrows(Nodef(bd1NodeF,:),2);
% [i,bd2NodeFIndex] = sortrows(Nodef(bd2NodeF,:),1);
% [i,bd3NodeFIndex] = sortrows(Nodef(bd3NodeF,:),2);
% [i,bd4NodeFIndex] = sortrows(Nodef(bd4NodeF,:),1);
[i,bd5NodeFIndex] = sortrows(NodeF(bd5NodeF,:),1);
[i,bd6NodeFIndex] = sortrows(NodeF(bd6NodeF,:),2);
[i,bd7NodeFIndex] = sortrows(NodeF(bd7NodeF,:),1);
[i,bd8NodeFIndex] = sortrows(NodeF(bd8NodeF,:),2);

gamFbd1 = bd1NodeF(bd1NodeFIndex);
% gamFbd2 = bd2NodeF(bd2NodeFIndex);
% gamFbd3 = bd3NodeF(bd3NodeFIndex);
% gamFbd4 = bd4NodeF(bd4NodeFIndex);
gamFbd5 = bd5NodeF(bd5NodeFIndex);
gamFbd6 = bd6NodeF(bd6NodeFIndex);
gamFbd7 = bd7NodeF(bd7NodeFIndex);
gamFbd8 = bd8NodeF(bd8NodeFIndex);

gamFbd1 = [gamFbd1(1:length(gamFbd1)-1) gamFbd1(2:length(gamFbd1));gamFbd1(1) gamFbd1(length(gamFbd1))];
% gamFbd2 = [gamFbd2(1:length(gamFbd2)-1) gamFbd2(2:length(gamFbd2))];
% gamFbd3 = [gamFbd3(1:length(gamFbd3)-1) gamFbd3(2:length(gamFbd3))];
% gamFbd4 = [gamFbd4(1:length(gamFbd4)-1) gamFbd4(2:length(gamFbd4))];
gamFbd5 = [gamFbd5(1:length(gamFbd5)-1) gamFbd5(2:length(gamFbd5))];
gamFbd6 = [gamFbd6(1:length(gamFbd6)-1) gamFbd6(2:length(gamFbd6))];
gamFbd7 = [gamFbd7(1:length(gamFbd7)-1) gamFbd7(2:length(gamFbd7))];
gamFbd8 = [gamFbd8(1:length(gamFbd8)-1) gamFbd8(2:length(gamFbd8))];

% gamFbd8 = setdiff(sort(gamFbd8,2),sort(gamFbd4,2),'rows');

n1 = size(gamFbd1,1);
% n2 = size(gamFbd2,1); % these lines are not required for the bench mark
% problem
% n3 = size(gamFbd3,1);
% n4 = size(gamFbd4,1); % comment this line
n5 = size(gamFbd5,1);
n6 = size(gamFbd6,1);
n7 = size(gamFbd7,1);
n8 = size(gamFbd8,1);

% gamF = [gamFbd1;gamFbd2;gamFbd3;gamFbd4;gamFbd5;gamFbd6;gamFbd7;gamFbd8];
% gamIF = [1 n1;n1+1 n2+n1;n1+n2+1 n1+n2+n3;n1+n2+n3+1 n1+n2+n3+n4;n1+n2+n3+n4+1 n1+n2+n3+n4+n5;...
%     n1+n2+n3+n4+n5+1 n1+n2+n3+n4+n5+n6; n1+n2+n3+n4+n5+n6+1 n1+n2+n3+n4+n5+n6+n7; n1+n2+n3+n4+n5+n6+n7+1 n1+n2+n3+n4+n5+n6+n7+n8];

% gamF = [gamFbd1;gamFbd2;gamFbd3;gamFbd5;gamFbd6;gamFbd7;gamFbd8];
% gamIF = [1 n1;n1+1 n2+n1;n1+n2+1 n1+n2+n3;n1+n2+n3+1 n1+n2+n3+n5;...
%     n1+n2+n3+n5+1 n1+n2+n3+n5+n6; n1+n2+n3+n5+n6+1 n1+n2+n3+n5+n6+n7;n1+n2+n3+n5+n6+n7+1 n1+n2+n3+n5+n6+n7+n8];

gamF = [gamFbd1;gamFbd5;gamFbd6;gamFbd7;gamFbd8];
gamIF = [1 n1;n1+1 n1+n5;n1+n5+1 n1+n5+n6;n1+n5+n6+1 n1+n5+n6+n7;n1+n5+n6+n7+1 n1+n5+n6+n7+n8];

% carefull while entering the values because we have a data structure where
% Neummann nodes are entered first and Dirchlet nodes are entered second.
% iF = [1 2 3 5];
% cF = [];
% 
% dF = [4 5 7 8];
% nF = [1 2 3 6];

iF = [2];
cF = [2];

dF = [1 2 4 5];
nF = [3];

DF = gamIF(dF',:);
NF = gamIF(nF',:);

save('FlexPlate\BenchMarkDataF.mat','DF','NF','cF','iF','gamIF','gamF','ElemF','NodeF');
save('FlexPlate\BenchMarkDataS.mat','DS','NS','cS','iS','gamIS','gamS','ElemS','NodeS');


function [x,y]= polar2cart (mag, ang_in_deg)
 x = mag .* cos(ang_in_deg*pi/180);
 y = mag .* sin(ang_in_deg*pi/180);
  
 
function [r,ad] = cart2polar(x)
 x = x(:,1)+sqrt(-1)*x(:,2);
 r = abs(x);
 ar = angle(x);
 ad = ar*180/pi;

