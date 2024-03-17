format long;

fileName1 = 'DataP';
fileName2 = 'DataU';
OutPutFileName = 'XYDiffPlot';
FileTitle = 'XY Plot';
timeStep = 0.002;
FileStep = 0.1;
TimeStep = 0.02;
dataSize = 5; 
numTime = 50;
TotalTime = 15;
Result = [];

load('SolverInputs.mat','ElemS','DeEF','ElemF','Phy');

level = unique(ElemF.dof(intersect(ElemF.itfNode,unique(ElemF.elem2dof(ElemF.itfElem,1:3))),1));
level = sort(level,'descend');
Nodes = find(ElemS.dof(:,2)==-0.005 & ElemS.dof(:,1)==1);

Nodes1 = find(ElemS.dof(:,2)==0);

NumItn = ceil(TotalTime/(FileStep));
f = zeros(length(level)-1,2);
for i =1:length(level)-1
    FileName1 = [fileName1 num2str(TotalTime/FileStep) '.mat'];
    load(FileName1);
    FileName2 = [fileName2 num2str(TotalTime/FileStep) '.mat'];
    load(FileName2);
    f(i,:) = liftDrag(ElemF,Phy,1,DataU,DataP,DeEF,level(i)); % Remember DataU and DataP over here contains only one time step
end
f = [0 0;f];
plot(level,-f(:,1));

% outFileName = [OutPutFileName '.plt'];
% fileID = fopen(outFileName,'w');
% fprintf(fileID,' TITLE = \"%s\"\n',FileTitle);
% if strcmp(FileTitle,'XY Plot')
%         fprintf(fileID,' VARIABLES = \"TIME\", \"TIP DISPLACEMENT\"\n');
% end
% fprintf(fileID,' ZONE I=%d , J=%d, DATAPACKING=POINT\n',[2 size(Result,1)/2]);
% fprintf(fileID,' %f %f\n',Result');
% fclose(fileID);
% figure(10)
% plot(Result(:,1),Result(:,2));
% figure(9)
% plot(nodes1(:,1),nodes1(:,2),'+');