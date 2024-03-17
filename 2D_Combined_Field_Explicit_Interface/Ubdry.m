function g=Ubdry(xy,str,t,Phy)
if (strcmp(str,'FlexFlatPlateF')==1) ...
     ||(strcmp(str,'FlexFlatPlateTh02F')==1)...
            ||(strcmp(str,'FlexFlatPlateTh005F')==1)
    if t<=0
        f=0;
    elseif t<=.5
        f=Phy.Uo*(1-cos(pi/(2*5e-1)*t));
    else
        f=Phy.Uo;
    end
 elseif strcmp(str(1:4),'FS4f')==1
     if t<=0
         f = 0;
     elseif t<=2
        f = Phy.Uo*(1-cos(pi/(2)*t/2));
     else
         f = Phy.Uo;
     end
%     f = Phy.Uo;
elseif strcmp(str(1:4),'FPCf')==1
     if t<=0
         f = 0;
     elseif t<=2
        f = Phy.Uo*(1-cos(pi/(2)*t))/2;
     else
         f = Phy.Uo;
     end
else
    disp('not implement Ubdry');
end
        

x=xy(:,1); y=xy(:,2); N=length(x);
u=zeros(N,1); 

if strcmp(str(1:4),'FS4f')==1
    n1=find(abs(x+0)<1e-8);
    % u(n1)=2*f*(y(n1)).*(0.02-y(n1));
    u(n1)=f*(1+2*y(n1)).*(1-y(n1));
elseif strcmp(str(1:4),'FPCf')==1
    n1=find(abs(x+0)<1e-8);
    % n1=[n1;find(abs(x-2.5)<1e-8)];
    u(n1)=1.5*4/0.1681*f*(y(n1)).*(0.41-y(n1));
elseif strcmp(str(1:4),'HVf1')==1
    f=2*sin(2*pi*t);
    n1=find(abs(x+0)<1e-8);
    u(n1)=f*(0.7-y(n1)).*(0.7+y(n1));    
elseif strcmp(str,'FlexFlatPlateF')==1 ...
        ||(strcmp(str,'FlexFlatPlateTh02F')==1)...
            ||(strcmp(str,'FlexFlatPlateTh005F')==1) 
%     n1 = find (abs(x+2)<=1e-8);
%     u(n1)=f*(5-y(n1)).*(5+y(n1))/25;
%     n1 = find (abs(y+0.005)<=1e-8);
%     u(n1)= f*(5-y(n1)).*(5+y(n1))/25;
%     n1 = find (abs(y-0.005)<=1e-8);
%     u(n1)= f*(5-y(n1)).*(5+y(n1))/25;
    n1 = find (abs(x+2)<=1e-8);
    u(n1)= f;
    n1 = find (abs(y+Phy.thk/2)<=1e-8);
    u(n1)= f;
    n1 = find (abs(y-Phy.thk/2)<=1e-8);
    u(n1)= f;
    n1 = find (abs(y-5)<=1e-8);
    u(n1)= f;
    n1 = find (abs(y+5)<=1e-8);
    u(n1)= f;  
elseif strcmp(str,'FlexCurvedPlateF')==1
    n1 = find (abs(x+2)<=1e-8);
    u(n1)=f*(1-y(n1)).*(1+y(n1));
    n1 = find (abs(y+0.005)<=1e-8);
    u(n1)= f*(1-y(n1)).*(1+y(n1));
    n1 = find (abs(y-0.005)<=1e-8);
    u(n1)= f*(1-y(n1)).*(1+y(n1));
elseif strcmp(str,'BenchMarkF')==1
    n1 = find (abs(x+0)<=1e-8);
    u(n1)=4*f*(0.41-y(n1)).*y(n1)/0.41^2;
else
    disp('Ubdry is undefined');
    pause;
end

g=[u,zeros(N,1)];