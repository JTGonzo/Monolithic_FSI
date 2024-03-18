// medium
mesh_s1 = 0.75;
mesh_s2 = 0.2;
mesh_s3 = 0.015;

// coarse
mesh_s1 = 0.8;
mesh_s2 = 0.4;
mesh_s3 = 0.075;

H  = 4;
L  = 6;
l  = 0.1; 
h = 1;

Point(1) = {0, 0, 0, mesh_s2};
Point(2) = {0, H, 0, mesh_s2};
Point(3) = {L, H, 0, mesh_s2};
Point(4) = {L, 0,    0, mesh_s2};

Point(5) = {L/2 + l/2, 0, 0, mesh_s3};
Point(6) = {L/2 + l/2, h, 0, mesh_s3};
Point(7) = {L/2 - l/2, h, 0, mesh_s3};
Point(8) = {L/2 - l/2, 0, 0, mesh_s3};


Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};


Line(9) = {5, 8};

Line Loop(10) = {8, 1, 2, 3, 4, 5, 6, 7};
Line Loop(11) = {9, -5, -6, -7};

Plane Surface(12) = {10};
Plane Surface(13) = {11};

