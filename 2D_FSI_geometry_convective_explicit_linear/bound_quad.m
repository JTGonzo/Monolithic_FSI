function [x,w,phi] = bound_quad(np,a,b)

n=np-1;
phi   = [];

if np<=1
    x=0;w=2;
    return
end

x=zeros(np,1);
w=zeros(np,1);
jac=zeros(np);

k=[1:n];
v=(k)./(sqrt(4*(k.^2)-1));
jac=jac+diag(v,1)+diag(v,-1);

[w1,x]=eig(jac);
norm2=sqrt(diag(w1'*w1));    
w1=(2*w1(1,:)'.^2)./norm2;  
x=diag(x);		    

[x,ip]=sort(x);

for i=1:np
    w(i)=w1(ip(i));
end

if nargin == 3
    bma=(b-a)*.5;
    bpa=(b+a)*.5;
    x=bma*x+bpa;
    w=w*bma;
end

x = x';
w = w';

X = [x; 0*x];

x = X(1,:);
y = X(2,:);

phi(1,:) = (1-x-y).*(1-2*x-2*y);
phi(2,:) = x.*(-1+2*x);
phi(3,:) = y.*(-1+2*y);
phi(4,:) = 4*x.*(1-x-y);
phi(5,:) = 4*x.*y;
phi(6,:) = 4*y.*(1-x-y);

BDOF = [1 2 4];

phi   = phi(BDOF,:);

end