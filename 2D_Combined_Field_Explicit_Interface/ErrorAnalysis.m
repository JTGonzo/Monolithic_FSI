function ErrorAnalysis()
clc;
fileName = {'eta01';'eta005';'eta0025';'eta00125';'eta000625'};
fileNameU = {'U01' 'U005' 'U0025' 'U00125' 'U000625'};
fileNameP = {'P01' 'P005' 'P0025' 'P00125' 'P000625'};
fileNameV = {'V01' 'V005' 'V0025' 'V00125' 'V000625'};

time = [0.1 0.05,0.025,0.0125];%,0.00625];
Eta = [];
no = 4;
for j=1:no
    str = [cell2mat(fileName(j)) '.mat'];
    strU = [cell2mat(fileNameU(j)) '.mat'];
    strP = [cell2mat(fileNameP(j)) '.mat'];
    strV = [cell2mat(fileNameV(j)) '.mat'];

    load(str);
    load(strU);
    load(strP);
    load(strV);

    Eta(:,:,j)= eta(:,:,1);
    u(:,:,j)=U;
    p(:,:,j)=P;
    v(:,:,j)=V;

    load('ElemS.mat','ElemS');
    load('ElemF.mat','ElemF');

%     load('eta00125.mat');
%     load('U00125.mat');
%     load('P00125.mat');
%     load('V00125.mat');
     
    load('eta000625.mat');
    load('U000625.mat');
    load('P000625.mat');
    load('V000625.mat');

%     load('eta0001.mat');
%     load('U0001.mat');
%     load('P0001.mat');
% %     load('V0001.mat');

%     load('etaAct.mat');
%     load('UAct.mat');
%     load('PAct.mat');
%     load('VAct.mat');

end
    eta = repmat(eta(:,:,1),[1,1,no]);
    U = repmat(U,[1,1,no]);
    P = repmat(P,[1,1,no]);
    V = repmat(V,[1,1,no]);


    E = Eta(ElemS.itfNode,:,:) - eta(ElemS.itfNode,:,:);
    ui = u(ElemF.itfNode,:,:) - U(ElemF.itfNode,:,:);
    vi = v(ElemS.itfNode,:,:) - V(ElemS.itfNode,:,:);

%     P = P(ElemF.itfNode,:,:) - P(ElemF.itfNode,:,:);
    E1 = Eta-eta;
    E1U = u-U;
    E1V = v-V;
%     E1 = Eta-eta;
    L2 = sqrt(sum(E1.^2,1));
    Linf = max(abs(E1));
    
    L2U = sqrt(sum(E1U.^2,1));
    LinfU = max(abs(E1U));
    
    L2V = sqrt(sum(E1V.^2,1));
    LinfV = max(abs(E1V));
        
    disp('log(L2) error is');
    Err = log(L2)
    Einf = log(Linf)
    
    ErrU = log(L2U)
    EinfU = log(LinfU)
    
    ErrV = log(L2V)
    EinfV = log(LinfV)
    
    L2Order = [];
    L2UOrder = [];
    L2VOrder = [];
    
    LInfOrder = [];
    LInfUOrder = [];
    LInfVOrder = [];
    for i=1:size(Err,3)-1
        orderL2 = (Err(:,:,i)-Err(:,:,i+1))/log(2);
        orderLinf = (Einf(:,:,i)-Einf(:,:,i+1))/log(2);
        
        orderL2U = (ErrU(:,:,i)-ErrU(:,:,i+1))/log(2);
        orderLinfU = (EinfU(:,:,i)-EinfU(:,:,i+1))/log(2);
        
        orderL2V = (ErrV(:,:,i)-ErrV(:,:,i+1))/log(2);
        orderLinfV = (EinfV(:,:,i)-EinfV(:,:,i+1))/log(2);
        
        disp('Order of the L2 Solution is');
        L2Order = [L2Order;orderL2];
        L2UOrder = [L2UOrder;orderL2U];
        L2VOrder = [L2VOrder;orderL2V];
        disp('Order of the Linf Solution is');
        LInfOrder = [LInfOrder;orderLinf];
        LInfUOrder = [LInfUOrder;orderLinfU];
        LInfVOrder = [LInfVOrder;orderLinfV];
    end
    L2Order
    LInfOrder
    
    L2UOrder
    LInfUOrder
    
    L2VOrder
    LInfVOrder
    
%     L2Bd = sqrt(sum(E.^2,2));
%     [r,d] = cart2polar(ElemS.dof(ElemS.itfNode,:)-repmat([1.5 -0.5],size(ElemS.itfNode,1),1));
%     for i=1:size(L2Bd,3)
%         plot(d,L2Bd(:,1,i))
%         hold on;
%     end
    
    figure(2)
    data = reshape(L2(1,1,:),[],1);
    loglog(time,data)
    
    xlabel('Time Step, \Delta T');
    ylabel('L^2 Error');

    h_xlabel = get(gca,'XLabel');
    set(h_xlabel,'FontSize',14); 

    h_ylabel = get(gca,'YLabel');
    set(h_ylabel,'FontSize',14);
    
%     hold on;
%     data = reshape(Linf(1,1,:),[],1);
%     plot(time,data)
    
%     figure(3)
    hold on;
    data = reshape(L2U(1,1,:),[],1);
    loglog(time,data)
    
%     xlabel('Time Step, \Delta T');
%     ylabel('L^2 Error');
% 
%     h_xlabel = get(gca,'XLabel');
%     set(h_xlabel,'FontSize',14); 
% 
%     h_ylabel = get(gca,'YLabel');
%     set(h_ylabel,'FontSize',14);
    
% %     figure(4)
%     hold on;
%     data = reshape(L2V(1,1,:),[],1);
%     loglog(time,data)
    
    hold off;
    
%     xlabel('Time Step, \Delta T');
%     ylabel('L^2 Error');
% 
%     h_xlabel = get(gca,'XLabel');
%     set(h_xlabel,'FontSize',14); 
% 
%     h_ylabel = get(gca,'YLabel');
%     set(h_ylabel,'FontSize',14);
    
function [x,y]= polar2cart (mag, ang_in_deg)
 x = mag .* cos(ang_in_deg*pi/180);
 y = mag .* sin(ang_in_deg*pi/180);
  
 
function [r ad] = cart2polar(x)
 j = sqrt(-1);
 x = x(:,1)+j*x(:,2);
 r = abs(x);
 ar = angle(x);
 ad = ar*180/pi;