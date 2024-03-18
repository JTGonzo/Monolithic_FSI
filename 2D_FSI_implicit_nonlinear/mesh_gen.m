%   This file is part of redbKIT.
%   Copyright (c) 2016, Ecole Polytechnique Federale de Lausanne (EPFL)
%   Author: Federico Negri <federico.negri@epfl.ch>

clear all; close all; clc;

[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');

dim      =  2;

%% Load meshes
[mshS.vertices, mshS.boundaries, mshS.elements, mshS.rings] = msh_to_Mmesh('mesh/mesh_Solid', dim);
vertices = mshS.vertices;
boundaries = mshS.boundaries;
elements = mshS.elements; 
rings = mshS.rings;
save('flap_S.mat','vertices', 'boundaries', 'elements', 'rings')


[mshF.vertices, mshF.boundaries, mshF.elements, mshF.rings] = msh_to_Mmesh('mesh/mesh_Fluid', dim);
vertices = mshF.vertices;
boundaries = mshF.boundaries;
elements = mshF.elements; 
rings = mshF.rings;
save('flap_F.mat','vertices', 'boundaries', 'elements', 'rings')
