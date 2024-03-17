clc;

format long;

fileName = 'DataU';
DispFileName = 'DataME';
OutPutFileName = 'Vorticity';
FileTitle = 'vorticity plot';
timeStep = 0.002;
TimeStep = 0.02;
dataSize = 5; 
numTime = 50;
ResultAt = [15];

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
    
    DeEF = deformedMeshData(dof,ElemF);
    
    b1=fU_gradp(DeEF,ElemF,u(:,1));
    b2=fU_gradp(DeEF,ElemF,u(:,2));
    [~,fM,~]=initMatrix2old(DeEF,ElemF);
    vort=fM\(b2(:,1)-b1(:,2));
    newVort = vort;
    newNodes=dof;
    newElem = [ElemF.elem2dof(:,[1 6 5]);ElemF.elem2dof(:,[6 2 4]);ElemF.elem2dof(:,[6 4 5]);ElemF.elem2dof(:,[5 4 3])];
%     for al=1:2
%         [newNodes, newElem, newVort] = refinement(newNodes, newElem, newVort);
%     end
    
    Result = [newNodes newVort];
    m = size(Result,1);
    n = size(newElem,1);
        
    outFileName = [OutPutFileName int2str(i) '.plt'];
    fileID = fopen(outFileName,'w');
    fprintf(fileID,' TITLE = \"%s\"\n',FileTitle);
    if strcmp(FileTitle,'vorticity plot')
            fprintf(fileID,' VARIABLES = \"X\", \"Y\", \"omega_Z\"\n');
    end
    fprintf(fileID,' ZONE NODES=%d , ELEMENTS=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',m,n);
    fprintf(fileID,'%f %f %f\n',Result');
    
    fprintf(fileID,'%d %d %d\n',newElem');
    
    fclose(fileID); 
end

