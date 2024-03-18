Merge "meshFS.msh";

Physical Line(2) = {1};// Inlet Fluid
Physical Line(3) = {3};// Outlet Fluid
Physical Line(4) = {2, 4, 8};// lateral Fluid
Physical Line(5) = {7, 6, 5};// Beam
//Physical Line(6) = {9}; // Obstacle-beam interface

//Physical Surface(2) = {13};// Solid
Physical Surface(1) = {12};// Fluid

Physical Point(1) = {5, 8}; // rings
