load('SOL.mat');
LiftDrag = Sol.ldc;
n = size(LiftDrag,1);
n = [1:1:n]';
x = [n LiftDrag(:,1)];
x = [0 0;x];
figure(1)
max(x(4000:4200,2))
min(x(4000:4200,2))
plot(x(:,1)*0.005,x(:,2),'k','LineWidth',2);
xlabel('Time, Sec');
ylabel('Drag Force, N');

h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',30); 

h_ylabel = get(gca,'YLabel');
set(h_ylabel,'FontSize',30);

set(gca,'FontSize',24)


y = [n LiftDrag(:,2)];
y = [0 0;y];
figure(2)
max(y(4000:4200,2))
min(y(4000:4200,2))
plot(y(:,1)*0.005,y(:,2),'k','LineWidth',2);
xlabel('Time, Sec');
ylabel('Lift Force, N');

h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',30); 

h_ylabel = get(gca,'YLabel');
set(h_ylabel,'FontSize',30);

set(gca,'FontSize',24)
