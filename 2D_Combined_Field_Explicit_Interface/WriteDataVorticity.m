clc;

format long;

fileName = 'DataU';
DispFileName = 'DataME';
OutPutFileName = 'Velocity';
FileTitle = 'vorticity plot';
timeStep = 0.002;
TimeStep = 0.02;
dataSize = 5; 
numTime = 50;
ResultAt = 15;

Result = [];

load('SolverInputs.mat','ElemF');
for i=1:length(ResultAt)
    time = ResultAt(i);
    fileNum = time/(timeStep*numTime);
    FileNum = ceil(fileNum);
    loadFile = [fileName int2str(FileNum) '.mat'];
    DispFile = [DispFileName int2str(FileNum) '.mat'];
    load(loadFile);
    load(DispFile);
    
    DataNum = mod(fileNum,TimeStep*dataSize)/timeStep;
    if DataNum == 0
        DataNum = dataSize;
    end
    DataNum = round(DataNum);
    dof = ElemF.dof+DataME(:,:,DataNum);
    u = DataU(:,:,DataNum);
    Result = [dof u];
    m = size(Result,1);
    n = size(ElemF.elem2dof,1)*4;
        
    outFileName = [OutPutFileName int2str(i) '.plt'];
    fileID = fopen(outFileName,'w');
    fprintf(fileID,' TITLE = \"%s\"\n',FileTitle);
    if strcmp(FileTitle,'vorticity plot')
            fprintf(fileID,' VARIABLES = \"X\", \"Y\", \"Ux\", \"Uy\"\n');
    end
    fprintf(fileID,' ZONE NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',m,n);
    fprintf(fileID,'%f %f %f %f\n',Result');
    
    fprintf(fileID,'%d %d %d\n',[ElemF.elem2dof(:,[1 6 5])]');
    fprintf(fileID,'%d %d %d\n',[ElemF.elem2dof(:,[6 2 4])]');
    fprintf(fileID,'%d %d %d\n',[ElemF.elem2dof(:,[6 4 5])]');
    fprintf(fileID,'%d %d %d\n',[ElemF.elem2dof(:,[5 4 3])]');

    fclose(fileID); 
end

