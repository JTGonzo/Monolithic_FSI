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

Y = fft(Result(:,2));
n = length(Y);
power = abs(Y).^2;
nyquist = 1/2;
freq = (1:n/2)/(n/2)*nyquist;
plot(freq,power);