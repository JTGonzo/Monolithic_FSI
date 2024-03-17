close all
clear variables;
set(0,'defaultlinelinewidth',2,'defaultpatchlinewidth',1)
set(0,'defaulttextfontsize',14,'defaultaxesfontsize',14)


rand('state',100);

a=1.75;
elem2dof=[1,2,3];
%dof0(:,:,1)=[1,0;2,0;1,1].*0.8;
dof0(:,:,1)=[1,0;2,0;1,1].*0.8+[0,a;0,a;0,a];
dof0(:,:,2)=[1,0;2,0;1,1].*0.8+[0,2*a;0,2*a;0,2*a];

%dof1(:,:,1)=[-1,0;0,0.2;-0.4,0.7];
dof1(:,:,1)=[-1,0;0,0;-0.4,0.7]+[0,a;0,a;0,a]+randn(3,2)*0.1+[-0.5,0;-0.5,0;-0.5,0];
dof1(:,:,2)=[-1,0;0,0;-0.4,0.7]+[0,2*a;0,2*a;0,2*a]+randn(3,2)*0.1+[-0.2,0;-0.2,0;-0.2,0];

dof4(1:2,:,:)=dof1(2:3,:,:);
%dof4(3,:,1)=[0.35,0.8];
dof4(3,:,1)=[0.35-0.5,0.75]+[0,a];
dof4(3,:,2)=[0.35-0.25,0.8]+[0,2*a];


dof2(1,:,:)=dof1(3,:,:);
dof2(3,:,:)=dof1(1,:,:);
%dof2(2,:,1)=[-1.2,0.7];
dof2(2,:,1)=[-1.2-0.5,0.68]+[0,a];
dof2(2,:,2)=[-1.2-0.3,0.75]+[0,2*a];

n=20;
dof3=zeros(n+1,2,2);
dof3(1,:,:)=dof2(1,:,:);
for i=1:2
    a1=dof2(1,:,i);
    a2=dof2(2,:,i);
    a3=dof2(3,:,i);
    r1=norm(a2-a1);
    r2=norm(a3-a1);
    t1=atan2(a2(2)-a1(2),a2(1)-a1(1));
    if t1<0, t1=t1+2*pi; end
    t2=atan2(a3(2)-a1(2),a3(1)-a1(1));
    if t2<0, t2=t2+2*pi; end
    r=linspace(r1,r2,n)';
    t=linspace(t1,t2,n)';
    aa=[r.*cos(t),r.*sin(t)]+repmat(a1,n,1);
    dof3(2:end,:,i)=aa;
end

dof0(:,2,2)=dof0(:,2,2)+1;
dof1(:,2,2)=dof1(:,2,2)+1;
dof2(:,2,2)=dof2(:,2,2)+1;
dof3(:,2,2)=dof3(:,2,2)+1;
dof4(:,2,2)=dof4(:,2,2)+1;

dof0(:,:,3)=(dof0(:,:,1)+dof0(:,:,2))/2;
dof1(:,:,3)=(dof1(:,:,1)+dof1(:,:,2))/2;
dof2(:,:,3)=(dof2(:,:,1)+dof2(:,:,2))/2;
dof3(:,:,3)=(dof3(:,:,1)+dof3(:,:,2))/2;
dof4(:,:,3)=(dof4(:,:,1)+dof4(:,:,2))/2;

for i=1:3
    cen0(i,:)=[2/3,1/6,1/6]*dof0(:,:,i);
    cen1(i,:)=[2/3,1/6,1/6]*dof1(:,:,i);
    cen2(i,:)=[2/3,1/6,1/6]*dof2(:,:,i);    
    cen4(i,:)=[1/6,2/3,1/6]*dof4(1:3,:,i);    
end


trisurf(elem2dof,dof0(:,1,2),dof0(:,2,2),zeros(size(dof0,1),1)); hold on;
for i=1:3
%     trisurf(elem2dof,dof0(:,1,i),dof0(:,2,i),zeros(size(dof0,1),1)); hold on;
    trisurf(elem2dof,dof1(:,1,i),dof1(:,2,i),zeros(size(dof1,1),1)); hold on;
    fill(dof3(:,1,i),dof3(:,2,i),'g');
    fill(dof4(:,1,i),dof4(:,2,i),'g');
%     cen0(i,:)=[2/3,1/6,1/6]*dof0(:,:,i);
%     cen1(i,:)=[2/3,1/6,1/6]*dof1(:,:,i);
%     cen2(i,:)=[2/3,1/6,1/6]*dof2(1:3,:,i);    
%     cen4(i,:)=[1/6,2/3,1/6]*dof4(1:3,:,i);    
end

for i=1:5
    if i==4
        x = cen1(:,2)';
        y = cen1(:,1)';
    elseif i==5
        x=squeeze(dof2(3,2,:))';
        y=squeeze(dof2(3,1,:))';
    else
        x=squeeze(dof1(i,2,:))';
        y=squeeze(dof1(i,1,:))';
    end
    cs = spline(x(1:2),[y(1:2)]);
    xx = linspace(min(x),max(x),101);
    switch i
        case {1,2}
            plot(y,x,'r.',ppval(cs,xx),xx,'r-'); hold on;
        case {3,5}
            plot(y,x,'r.',ppval(cs,xx),xx,'r:'); hold on;
        case 4
            plot(y,x,'b.'); hold on; %,ppval(cs,xx),xx,'b:'); hold on;
    end
end

plot(cen0(2,1),cen0(2,2),'b.'); hold on;
plot(cen2(:,1),cen2(:,2),'b.'); hold on;
plot(cen4(:,1),cen4(:,2),'b.'); hold on;

for i=1:3
    text(dof0(i,1,2),dof0(i,2,2),'(0,0)'); hold on;
end


text(-2.25,2*a+0.3,'$t^n$','Interpreter','latex','FontSize',15); 
text(-2.25,a+0.3,'$t^{n-1}$','Interpreter','latex','FontSize',15); 


text(dof1(1,1,2)-0.35,dof1(1,2,2),'$a_{i_{(j;0)}}^{n}$','Interpreter','latex','FontSize',12);
text(dof1(2,1,2)+0.1,dof1(2,2,2),'$a_{i_{(j;1)}}^{n}$','Interpreter','latex','FontSize',12);
text(dof1(3,1,2),dof1(3,2,2)+0.15,'$a_{i_{(j;2)}}^{n}$','Interpreter','latex','FontSize',12);

text((dof1(1,1,2)-0.35+dof1(2,1,2)+0.1)/2,...
    (dof1(1,2,2)+dof1(2,2,2))/2,'$a_{i_{(j;5)}}^{n}$','Interpreter','latex','FontSize',12);
text((dof1(2,1,2)+0.1+dof1(3,1,2))/2,...
    (dof1(2,2,2)+dof1(3,2,2)+0.15)/2,'$a_{i_{(j;3)}}^{n}$','Interpreter','latex','FontSize',12);
text((dof1(3,1,2)+dof1(1,1,2)-0.35)/2,...
    (dof1(3,2,2)+0.15+dof1(1,2,2))/2,'$a_{i_{(j;4)}}^{n}$','Interpreter','latex','FontSize',12);

plot(dof1(1,1,2),dof1(1,2,2),'r*');
plot(dof1(2,1,2),dof1(2,2,2),'r*');
plot(dof1(3,1,2),dof1(3,2,2),'r*');

plot((dof1(1,1,2)+dof1(2,1,2))/2,...
    (dof1(1,2,2)+dof1(2,2,2))/2,'r*');
plot((dof1(2,1,2)+dof1(3,1,2))/2,...
    (dof1(2,2,2)+dof1(3,2,2))/2,'r*');
plot((dof1(3,1,2)+dof1(1,1,2))/2,...
    (dof1(3,2,2)+dof1(1,2,2))/2,'r*');

text(cen1(1,1)+0.1,cen1(1,2),'$x^{n-1}$','Interpreter','latex','FontSize',15);
text(cen1(2,1)+0.1,cen1(2,2),'$x^{n}$','Interpreter','latex','FontSize',15);
text(cen0(2,1)+0.1,cen0(2,2),'$\hat{x}$','Interpreter','latex','FontSize',15);
view(2); 
axis equal; axis off;



hold on;

dof0old=dof0;
dof1old=dof1;
dof2old=dof2;
dof4old=dof4;

dof0(:,1,:)=dof0(:,1,:)+4;
dof1(:,1,:)=dof1(:,1,:)+4;
dof2(:,1,:)=dof2(:,1,:)+4;
dof4(:,1,:)=dof4(:,1,:)+4;


for i=1:2
    cen0(i,:)=[2/3,1/6,1/6]*dof0(:,:,i);
    cen1(i,:)=[2/3,1/6,1/6]*dof1(:,:,i);
    cen2(i,:)=[2/3,1/6,1/6]*dof2(1:3,:,i);    
    cen4(i,:)=[1/6,2/3,1/6]*dof4(1:3,:,i);    
end

trisurf(elem2dof,dof0(:,1,2),dof0(:,2,2),zeros(size(dof0,1),1),'FaceColor','y'); hold on;
for i=1:2
    trisurf(elem2dof,dof1(:,1,i),dof1(:,2,i),zeros(size(dof1,1),1),'FaceColor','y'); hold on;
    cen0(i,:)=[2/3,1/6,1/6]*dof0(:,:,i);
    cen1(i,:)=[2/3,1/6,1/6]*dof1(:,:,i);
end

for i=1:3
    text(dof0(i,1,2),dof0(i,2,2),'(0,0)'); hold on;
end

x=[cen1(:,1);cen0(2,1)];
y=[cen1(:,2);cen0(2,2)];
plot(x,y,'b.'); hold on; %,ppval(cs,xx),xx,'b:'); hold on;

text(cen0(2,1)+0.1,cen0(2,2),'$\hat{x}$','Interpreter','latex','FontSize',15);

ax=x([2,4]-1);
ay=y([2,4]-1);
line(ax,ay);
text(mean(ax),mean(ay),'$x^{n-1}=\Psi^{f;n-1}_j(\hat{x})$','Interpreter','latex');%,'FontSize',12);
ax=x([3,4]-1);
ay=y([3,4]-1);
line(ax+0.4*rand(2,1),ay+0.4*rand(2,1));
text(mean(ax),mean(ay),'$x^{n}=\Psi^{f;n}_j(\hat{x})$','Interpreter','latex');%,'FontSize',12);

ax=x([3,2]-1);
ay=y([3,2]-1);
line(ax+0.4*rand(2,1),ay+0.4*rand(2,1));
text(mean(ax)-2,mean(ay),'$x^{n-1}=\Phi^n(x^{n};t^{n-1})$','Interpreter','latex');%,'FontSize',12);
text(mean(ax)-2,mean(ay)-1,'$=\Psi^{f;n-1}_j\circ (\Psi^{f;n}_j)^{-1}(x^n)$','Interpreter','latex');%,'FontSize',12);

view(2); 
axis equal; axis off;



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% figure;

dof0=dof0old;
dof1=dof1old;
dof2=dof2old;
dof4=dof4old;

dof0(:,1,:)=dof0(:,1,:)+8;
dof1(:,1,:)=dof1(:,1,:)+8;
dof2(:,1,:)=dof2(:,1,:)+8;
dof4(:,1,:)=dof4(:,1,:)+8;

% dof0(:,2,2)=dof0(:,2,2)+1;
% dof1(:,2,2)=dof1(:,2,2)+1;
% dof2(:,2,2)=dof2(:,2,2)+1;
% dof3(:,2,2)=dof3(:,2,2)+1;
% dof4(:,2,2)=dof4(:,2,2)+1;
% % dof4(:,2,2)=dof4(:,2,2)+1;


dof1(:,:,3)=(dof1(:,:,1)+dof1(:,:,2))/2;
dof2(:,:,3)=(dof2(:,:,1)+dof2(:,:,2))/2;
dof3(:,:,3)=(dof3(:,:,1)+dof3(:,:,2))/2;
dof4(:,:,3)=(dof4(:,:,1)+dof4(:,:,2))/2;

for i=1:3
    %cen0(i,:)=[2/3,1/6,1/6]*dof0(:,:,i);
    cen1(i,:)=[2/3,1/6,1/6]*dof1(:,:,i);
    cen2(i,:)=[2/3,1/6,1/6]*dof2(:,:,i);    
    cen4(i,:)=[1/6,2/3,1/6]*dof4(1:3,:,i);    
end

%trisurf(elem2dof,dof0(:,1,2),dof0(:,2,2),zeros(size(dof0,1),1)); hold on;
for i=1:3
    % trisurf(elem2dof,dof0(:,1,i),dof0(:,2,i),zeros(size(dof0,1),1)); hold on;
    trisurf(elem2dof,dof1(:,1,i),dof1(:,2,i),zeros(size(dof1,1),1),'FaceColor','y'); hold on;
    %fill(dof3(:,1,i),dof3(:,2,i),'g');
    %fill(dof4(:,1,i),dof4(:,2,i),'g');
    % cen0(i,:)=[2/3,1/6,1/6]*dof0(:,:,i);
    cen1(i,:)=[2/3,1/6,1/6]*dof1(:,:,i);
    %cen2(i,:)=[2/3,1/6,1/6]*dof2(1:3,:,i);    
    % cen4(i,:)=[1/6,2/3,1/6]*dof4(1:3,:,i);    
end


view(2); 
axis equal; axis off;

for i=1:5
    if i==4
        x = cen1(:,2)';
        y = cen1(:,1)';
    elseif i==5
        x=squeeze(dof2(3,2,:))';
        y=squeeze(dof2(3,1,:))';
    else
        x=squeeze(dof1(i,2,:))';
        y=squeeze(dof1(i,1,:))';
    end
    cs = spline(x,[y]);
    xx = linspace(min(x),max(x),101);
    switch i
        case {1,2}
            plot(y,x,'r.',ppval(cs,xx),xx,'r-'); hold on;
        case {3,5}
            plot(y,x,'r.',ppval(cs,xx),xx,'r:'); hold on;
        case 4
            plot(y,x,'b.'); hold on; %,ppval(cs,xx),xx,'b:'); hold on;
    end
end

plot(cen0(2,1),cen0(2,2),'b.'); hold on;
plot(cen2(:,1),cen2(:,2),'b.'); hold on;
plot(cen4(:,1),cen4(:,2),'b.'); hold on;


text(-2.25,3*a+0.3,'$t^n$','Interpreter','latex','FontSize',15); 
text(-2.25,2*a+0.3,'$t^{n-\frac{1}{2}}$','Interpreter','latex','FontSize',15); 
text(-2.25,a+0.3,'$t^{n-1}$','Interpreter','latex','FontSize',15); 
% 
% 
% text(dof1(1,1,2)-0.35,dof1(1,2,2),'$a_{i_{(j,0)}}^{n}$','Interpreter','latex','FontSize',12);
% text(dof1(2,1,2)+0.1,dof1(2,2,2),'$a_{i_{(j,1)}}^{n}$','Interpreter','latex','FontSize',12);
% text(dof1(3,1,2),dof1(3,2,2)+0.15,'$a_{i_{(j,2)}}^{n}$','Interpreter','latex','FontSize',12);
% 
% text((dof1(1,1,2)-0.35+dof1(2,1,2)+0.1)/2,...
%     (dof1(1,2,2)+dof1(2,2,2))/2,'$a_{i_{(j,5)}}^{n}$','Interpreter','latex','FontSize',12);
% text((dof1(2,1,2)+0.1+dof1(3,1,2))/2,...
%     (dof1(2,2,2)+dof1(3,2,2)+0.15)/2,'$a_{i_{(j,3)}}^{n}$','Interpreter','latex','FontSize',12);
% text((dof1(3,1,2)+dof1(1,1,2)-0.35)/2,...
%     (dof1(3,2,2)+0.15+dof1(1,2,2))/2,'$a_{i_{(j,4)}}^{n}$','Interpreter','latex','FontSize',12);

plot(dof1(1,1,2),dof1(1,2,2),'r*');
plot(dof1(2,1,2),dof1(2,2,2),'r*');
plot(dof1(3,1,2),dof1(3,2,2),'r*');

plot((dof1(1,1,2)+dof1(2,1,2))/2,...
    (dof1(1,2,2)+dof1(2,2,2))/2,'r*');
plot((dof1(2,1,2)+dof1(3,1,2))/2,...
    (dof1(2,2,2)+dof1(3,2,2))/2,'r*');
plot((dof1(3,1,2)+dof1(1,1,2))/2,...
    (dof1(3,2,2)+dof1(1,2,2))/2,'r*');

text(cen1(1,1)+0.1,cen1(1,2),'$x^{n-1}$','Interpreter','latex','FontSize',15);
text(cen1(2,1)+0.1,cen1(2,2),'$x^{n}$','Interpreter','latex','FontSize',15);
text(cen0(2,1)+0.1,cen0(2,2),'$\hat{x}$','Interpreter','latex','FontSize',15);
view(2); 
axis equal; axis off;


% dof1(:,1,:)=dof1(:,1,:)+4;
% dof3(:,1,:)=dof3(:,1,:)+4;
% dof4(:,1,:)=dof4(:,1,:)+4;
% 
% %trisurf(elem2dof,dof0(:,1,2),dof0(:,2,2),zeros(size(dof0,1),1)); hold on;
% for i=1:3
%     % trisurf(elem2dof,dof0(:,1,i),dof0(:,2,i),zeros(size(dof0,1),1)); hold on;
%     trisurf(elem2dof,dof1(:,1,i),dof1(:,2,i),zeros(size(dof1,1),1)); hold on;
%     fill(dof3(:,1,i),dof3(:,2,i),'g');
%     fill(dof4(:,1,i),dof4(:,2,i),'g');
%     % cen0(i,:)=[2/3,1/6,1/6]*dof0(:,:,i);
%     cen1(i,:)=[2/3,1/6,1/6]*dof1(:,:,i);
%     cen2(i,:)=[2/3,1/6,1/6]*dof2(1:3,:,i);    
%     % cen4(i,:)=[1/6,2/3,1/6]*dof4(1:3,:,i);    
% end


for k=1:3
    p=[dof1(1,1,k),dof1(1,2,k)+0.75];
    b0=[dof1(1,1,k),dof1(1,2,k)]; % a_0
    b1=[dof1(2,1,k),dof1(2,2,k)];
    b2=[dof1(3,1,k),dof1(3,2,k)];
    b5=[(dof1(1,1,k)+dof1(2,1,k))/2,...
        (dof1(1,2,k)+dof1(2,2,k))/2];
    b3=[(dof1(2,1,k)+dof1(3,1,k))/2,...
        (dof1(2,2,k)+dof1(3,2,k))/2];
    b4=[(dof1(3,1,k)+dof1(1,1,k))/2,...
        (dof1(3,2,k)+dof1(1,2,k))/2];
    
    for i=1:3
        if i==1
            pp=[p;b0];
            x=linspace(pp(1,1),pp(2,1),1000);
            y=linspace(pp(1,2),pp(2,2),1000);
        elseif i==2
            pp=[p;b5;b1];
            xx=pp(:,1); yy=pp(:,2);
            A=[xx.^2,xx,[1;1;1]];
            c=A\yy;
            x=linspace(min(xx),max(xx),1000);
            y=c(1)*x.^2+c(2)*x+c(3);
        else
            pp=[p;b4;b2];
            xx=pp(:,1); yy=pp(:,2);
            A=[xx.^2,xx,[1;1;1]];
            c=A\yy;
            x=linspace(min(xx),max(xx),1000);
            y=c(1)*x.^2+c(2)*x+c(3);
        end
        plot(x,y,'b:'); hold on;
    end
    
    % text(dof1(1,1,k)-0.35,dof1(1,2,k),'$a_{i_{(j,0)}}^{n}$','Interpreter','latex','FontSize',12);
    % text(dof1(2,1,k)+0.1,dof1(2,2,k),'$a_{i_{(j,1)}}^{n}$','Interpreter','latex','FontSize',12);
    % text(dof1(3,1,k),dof1(3,2,k)+0.15,'$a_{i_{(j,2)}}^{n}$','Interpreter','latex','FontSize',12);
    %
    % text((dof1(1,1,k)-0.35+dof1(2,1,k)+0.1)/2,...
    %     (dof1(1,2,k)+dof1(2,2,k))/2,'$a_{i_{(j,5)}}^{n}$','Interpreter','latex','FontSize',12);
    % text((dof1(2,1,k)+0.1+dof1(3,1,k))/2,...
    %     (dof1(2,2,k)+dof1(3,2,k)+0.15)/2,'$a_{i_{(j,3)}}^{n}$','Interpreter','latex','FontSize',12);
    % text((dof1(3,1,k)+dof1(1,1,k)-0.35)/2,...
    %     (dof1(3,2,k)+0.15+dof1(1,2,k))/2,'$a_{i_{(j,4)}}^{n}$','Interpreter','latex','FontSize',12);
    
end