clear all;
fileID = fopen('elements.txt');
C = textscan(fileID, '%d  %d  %d  %d  %d  %d');
fclose(fileID);
A1 = double(cell2mat(C(1,1)));
A2 = double(cell2mat(C(1,2)));
A3 = double(cell2mat(C(1,3)));
A4 = double(cell2mat(C(1,4)));
A5 = double(cell2mat(C(1,5)));
A6 = double(cell2mat(C(1,6)));
ElemF = [A4 A5 A6];

fileID = fopen('nodes.txt');
C = textscan(fileID, '%d %f %f');
fclose(fileID);
B1 = double(cell2mat(C(1,1)));
B2 = double(cell2mat(C(1,2)));
B3 = double(cell2mat(C(1,3)));
NodeF = [B2 B3];


% this one of the rows which is required to be modified for the data
% corrections
InfNode = sortrows(intersect(intersect(NodeF(1>=NodeF(:,1),:),NodeF(NodeF(:,1)>=(0),:),'rows'),...
    union(NodeF(NodeF(:,2)==0.005,:),NodeF(NodeF(:,2)== -0.005,:),'rows'),'rows'),[2 1]);
% InfNode = sortrows([InfNode;InfNode(1,1),0;InfNode(size(InfNode,1),1),0],[2 1]);
% innerNode_solid = [InfNode(:,1),zeros(size(InfNode(:,1)))];
innerNode_solid = [];
NodeS = sortrows(unique([InfNode;innerNode_solid],'rows'),[2 1]);


nt = size(unique(NodeS(:,1)),1)-1;
nnodes = size(NodeS,1);
ElemNo = (1:nt)';
ElemS = [ElemNo,ElemNo+nt+1,ElemNo+nt+2;ElemNo,ElemNo+1,ElemNo+nt+2];

bd1NodeS = find(1>=NodeS(:,1) & NodeS(:,1)>=(0) & NodeS(:,2)==0.005);
bd2NodeS = find(0.005>=NodeS(:,2) & NodeS(:,2)>=(-0.005) & NodeS(:,1)== 1);
bd3NodeS = find( 1>=NodeS(:,1) & NodeS(:,1)>=(0) & NodeS(:,2)== -0.005);
bd4NodeS = find(0.005>=NodeS(:,2) & NodeS(:,2)>=(-0.005) & NodeS(:,1)== 0);

[i,bd1NodeSIndex] = sortrows(NodeS(bd1NodeS,:),[1 2]);
[i,bd2NodeSIndex] = sortrows(NodeS(bd2NodeS,:),[1 2]);
[i,bd3NodeSIndex] = sortrows(NodeS(bd3NodeS,:),[1 2]);
[i,bd4NodeSIndex] = sortrows(NodeS(bd4NodeS,:),[1 2]);

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
cS = [];

dS = [1 2 3 4];
nS = [];

DS = gamIS(dS',:);
NS = gamIS(nS',:);

bd1NodeF = find(1>=NodeF(:,1) & NodeF(:,1)>=(0) & NodeF(:,2)==0.005);
bd2NodeF = find(0.005>=NodeF(:,2) & NodeF(:,2)>=(-0.005) & NodeF(:,1)== 1);
bd3NodeF = find( 1>=NodeF(:,1) & NodeF(:,1)>=(0) & NodeF(:,2)== -0.005);
bd4NodeF = find(0.005>=NodeF(:,2) & NodeF(:,2)>=(-0.005) & NodeF(:,1)== (0));
bd5NodeF = find( 10>=NodeF(:,1) & NodeF(:,1)>=(0) & NodeF(:,2)==1);
bd6NodeF = find(1>=NodeF(:,2) & NodeF(:,2)>=(-1) & NodeF(:,1)== 10);
bd7NodeF = find(10>=NodeF(:,1) & NodeF(:,1)>=(0) & NodeF(:,2)== -1);
bd8NodeF = find( 1>=NodeF(:,2) & NodeF(:,2)>=(-1) & NodeF(:,1)== (0));

[i,bd1NodeFIndex] = sortrows(NodeF(bd1NodeF,:),[1 2]);
[i,bd2NodeFIndex] = sortrows(NodeF(bd2NodeF,:),[1 2]);
[i,bd3NodeFIndex] = sortrows(NodeF(bd3NodeF,:),[1 2]);
[i,bd4NodeFIndex] = sortrows(NodeF(bd4NodeF,:),[1 2]);
[i,bd5NodeFIndex] = sortrows(NodeF(bd5NodeF,:),[1 2]);
[i,bd6NodeFIndex] = sortrows(NodeF(bd6NodeF,:),[1 2]);
[i,bd7NodeFIndex] = sortrows(NodeF(bd7NodeF,:),[1 2]);
[i,bd8NodeFIndex] = sortrows(NodeF(bd8NodeF,:),[1 2]);

gamFbd1 = bd1NodeF(bd1NodeFIndex);
gamFbd2 = bd2NodeF(bd2NodeFIndex);
gamFbd3 = bd3NodeF(bd3NodeFIndex);
gamFbd4 = bd4NodeF(bd4NodeFIndex);
gamFbd5 = bd5NodeF(bd5NodeFIndex);
gamFbd6 = bd6NodeF(bd6NodeFIndex);
gamFbd7 = bd7NodeF(bd7NodeFIndex);
gamFbd8 = bd8NodeF(bd8NodeFIndex);

gamFbd1 = [gamFbd1(1:length(gamFbd1)-1) gamFbd1(2:length(gamFbd1))];
gamFbd2 = [gamFbd2(1:length(gamFbd2)-1) gamFbd2(2:length(gamFbd2))];
gamFbd3 = [gamFbd3(1:length(gamFbd3)-1) gamFbd3(2:length(gamFbd3))];
gamFbd4 = [gamFbd4(1:length(gamFbd4)-1) gamFbd4(2:length(gamFbd4))];
gamFbd5 = [gamFbd5(1:length(gamFbd5)-1) gamFbd5(2:length(gamFbd5))];
gamFbd6 = [gamFbd6(1:length(gamFbd6)-1) gamFbd6(2:length(gamFbd6))];
gamFbd7 = [gamFbd7(1:length(gamFbd7)-1) gamFbd7(2:length(gamFbd7))];
gamFbd8 = [gamFbd8(1:length(gamFbd8)-1) gamFbd8(2:length(gamFbd8))];

gamFbd8 = setdiff(sort(gamFbd8,2),sort(gamFbd4,2),'rows');

n1 = size(gamFbd1,1);
n2 = size(gamFbd2,1);
n3 = size(gamFbd3,1);
% n4 = size(gamFbd4,1); % comment this line
n5 = size(gamFbd5,1);
n6 = size(gamFbd6,1);
n7 = size(gamFbd7,1);
n8 = size(gamFbd8,1);

% gamF = [gamFbd1;gamFbd2;gamFbd3;gamFbd4;gamFbd5;gamFbd6;gamFbd7;gamFbd8];
% gamIF = [1 n1;n1+1 n2+n1;n1+n2+1 n1+n2+n3;n1+n2+n3+1 n1+n2+n3+n4;n1+n2+n3+n4+1 n1+n2+n3+n4+n5;...
%     n1+n2+n3+n4+n5+1 n1+n2+n3+n4+n5+n6; n1+n2+n3+n4+n5+n6+1 n1+n2+n3+n4+n5+n6+n7; n1+n2+n3+n4+n5+n6+n7+1 n1+n2+n3+n4+n5+n6+n7+n8];

gamF = [gamFbd1;gamFbd2;gamFbd3;gamFbd5;gamFbd6;gamFbd7;gamFbd8];
gamIF = [1 n1;n1+1 n2+n1;n1+n2+1 n1+n2+n3;n1+n2+n3+1 n1+n2+n3+n5;...
    n1+n2+n3+n5+1 n1+n2+n3+n5+n6; n1+n2+n3+n5+n6+1 n1+n2+n3+n5+n6+n7;n1+n2+n3+n5+n6+n7+1 n1+n2+n3+n5+n6+n7+n8];

% carefull while entering the values because we have a data structure where
% Neummann nodes are entered first and Dirchlet nodes are entered second.
% iF = [1 2 3 5];
% cF = [];
% 
% dF = [4 5 7 8];
% nF = [1 2 3 6];

iF = [1 2 3];
cF = [];

dF = [4 6 7];
nF = [1 2 3 5];

DF = gamIF(dF',:);
NF = gamIF(nF',:);

save('FlexPlate\DataF.mat','DF','NF','cF','iF','gamIF','gamF','ElemF','NodeF','InfNode');
save('FlexPlate\DataS.mat','DS','NS','cS','iS','gamIS','gamS','ElemS','NodeS','InfNode');

