clc;
clear all;

load('ElemS.mat');
load('Eta.mat');
p = find(ElemS.dof(:,1)==0.6 & ElemS.dof(:,2)==0.2);

x = DataEta(p,1,:);
y = DataEta(p,2,:);

x = reshape(x,[],1);
y = reshape(y,[],1);

n = 1:size(y,1);
n = n';
figure(1)
max(y(4000:4200))
min(y(4000:4200))
plot(n*0.005,y,'k','LineWidth',2)
xlabel('Time, Sec');
ylabel('Y-displacement, m');

h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',30); 

h_ylabel = get(gca,'YLabel');
set(h_ylabel,'FontSize',30);

set(gca,'FontSize',24)

figure(2)
max(x(4000:4200))
min(x(4000:4200))
plot(n*0.005,x,'k','LineWidth',2)

xlabel('Time, Sec');
ylabel('X-displacement, m');

h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',30); 

h_ylabel = get(gca,'YLabel');
set(h_ylabel,'FontSize',30);

set(gca,'FontSize',24)
