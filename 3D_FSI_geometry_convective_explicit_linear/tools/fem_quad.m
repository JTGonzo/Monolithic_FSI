function [phi, csi, eta, w] = fem_quad( )

ax = [0 1 0];
ay = [0 0 1];
area = 0.5;
             
w = [9/40 0.132394152788506*ones(1,3) 0.125939180544827*ones(1,3)];
x = sum(ax)/3;
y = sum(ay)/3;
p1 = 0.059715871789770;
p2 = 0.470142064105115;

x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
p1 = 0.797426985353087;
p2 = 0.101286507323456;

x = [x, p1*ax(1)+p2*ax(2)+p2*ax(3) p2*ax(1)+p1*ax(2)+p2*ax(3) p2*ax(1)+p2*ax(2)+p1*ax(3)];
y = [y, p1*ay(1)+p2*ay(2)+p2*ay(3) p2*ay(1)+p1*ay(2)+p2*ay(3) p2*ay(1)+p2*ay(2)+p1*ay(3)];
        
w = w*area;        
quad_points = [x;y];

csi = quad_points(1,:);
eta = quad_points(2,:);

phi   = [];

X = [csi; eta; 0*eta];

x = X(1,:);
y = X(2,:);
z = X(3,:);
        
phi(1,:) = (1-x-y-z).*(1-2*x-2*y-2*z);
phi(2,:) = x.*(2*x-1);
phi(3,:) = y.*(2*y-1);
phi(4,:) = z.*(2*z-1);
phi(5,:) = 4*x.*(1-x-y-z);
phi(6,:) = 4*x.*y;
phi(7,:) = 4*y.*(1-x-y-z);
phi(8,:) = 4*z.*(1-x-y-z);
phi(9,:) = 4*x.*z;
phi(10,:) = 4*y.*z;

BDOF = [1 2 3 5 6 7];

phi   = phi(BDOF,:);

end
