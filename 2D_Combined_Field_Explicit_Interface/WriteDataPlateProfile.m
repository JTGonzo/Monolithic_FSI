format long;

fileName = 'DataEta';
OutPutFileName = 'XYPlot';
FileTitle = 'XY Plot Vs Time';

StartTime = 5;
EndTime = 7;

FileStep = 0.1;

timeStep = 0.02;
dataSize = 5; 
numTime = 50;
TotalTime = 10;
ResultAt = [StartTime:0.02:EndTime];
Result = [];

load('SolverInputs.mat','ElemS');

Nodes = find(ElemS.dof(:,2)==0);

load('SolverInputs.mat','ElemF');
for i=1:length(ResultAt)
    time = ResultAt(i);
    fileNum = time/(FileStep);
    FileNum = ceil(fileNum);
    loadFile = [fileName int2str(FileNum) '.mat'];
    load(loadFile);

    DataNum = mod(fileNum,timeStep*dataSize)/timeStep;
    if DataNum == 0
        DataNum = dataSize;
    end
    DataNum = round(DataNum);
    
    dof = ElemS.dof(Nodes,:)+DataEta(Nodes,:,DataNum);
    Result = [dof];
    m = size(Result,1);
    n = size(ElemS.elem2dof,1);
    Result = sortrows(Result,1);
    
    plot(Result(:,1),Result(:,2))
    hold on;
%     outFileName = [OutPutFileName int2str(i) '.plt'];
%     fileID = fopen(outFileName,'w');
%     fprintf(fileID,' TITLE = \"%s\"\n',FileTitle);
%     if strcmp(FileTitle,'XY Plot')
%             fprintf(fileID,' VARIABLES = \"X\", \"Y\"\n');
%     end
%     fprintf(fileID,' ZONE I=%d , J=%d, DATAPACKING=POINT\n',[2 size(Result,1)/2]);
%     fprintf(fileID,' %f %f\n',Result');
%     fclose(fileID);
end

