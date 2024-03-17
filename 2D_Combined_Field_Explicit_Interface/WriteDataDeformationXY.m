format long;

fileName = 'DataEta';
OutPutFileName = 'XYDiffPlot';
FileTitle = 'XY Plot';
timeStep = 0.002;
FileStep = 0.1;
TimeStep = 0.02;
dataSize = 5; 
numTime = 50;
TotalTime = 15;
Result1 = [];
Result2 = [];

load('SolverInputs.mat','ElemS');

Nodes = find(ElemS.dof(:,2)==0 & ElemS.dof(:,1)==1);

NumItn = ceil(TotalTime/(FileStep));

for i =1:NumItn
    FileName = [fileName num2str(i) '.mat'];
    load(FileName);
    nodesy = DataEta(Nodes,2,:);
    nodesx = DataEta(Nodes,1,:);
    time = (TimeStep*(1:dataSize)')+dataSize*TimeStep*(i-1);
    Result1 = [Result1;time nodesy(:)];
    Result2 = [Result2;time nodesx(:)];
end

% outFileName = [OutPutFileName '.plt'];
% fileID = fopen(outFileName,'w');
% fprintf(fileID,' TITLE = \"%s\"\n',FileTitle);
% if strcmp(FileTitle,'XY Plot')
%         fprintf(fileID,' VARIABLES = \"TIME\", \"TIP DISPLACEMENT\"\n');
% end
% fprintf(fileID,' ZONE I=%d , J=%d, DATAPACKING=POINT\n',[2 size(Result,1)/2]);
% fprintf(fileID,' %f %f\n',Result');
% fclose(fileID);

plot(Result1(:,1),Result1(:,2));
hold on;
plot(Result2(:,1),Result2(:,2));
hold off;