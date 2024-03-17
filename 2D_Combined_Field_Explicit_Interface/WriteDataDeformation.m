format long;

fileName = 'DataEta';
OutPutFileName = 'XYDiffPlot';
FileTitle = 'XY Plot';
timeStep = 0.002;
FileStep = 0.1;
TimeStep = 0.02;
dataSize = 5; 
numTime = 50;
TotalTime = 10;
Result = [];

load('SolverInputs.mat','ElemS');

Nodes = find(ElemS.dof(:,2)==0 & ElemS.dof(:,1)==1);

NumItn = ceil(TotalTime/(FileStep));
for i =1:NumItn
    FileName = [fileName num2str(i) '.mat'];
    load(FileName);
    nodes = DataEta(Nodes,2,:);
    time = (TimeStep*(1:dataSize)')+dataSize*TimeStep*(i-1);
    Result = [Result;time nodes(:)];
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

plot(Result(:,1),Result(:,2));