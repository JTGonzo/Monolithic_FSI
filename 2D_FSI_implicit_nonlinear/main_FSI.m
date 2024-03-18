% Raw_FSI_P2
clear all; close all; clc;
[~,~,~] = mkdir('Figures');
[~,~,~] = mkdir('Results');
addpath('C_Files/')
addpath('import_export_time')
dim      =  2;

load('flap_S.mat')
meshSolid.boundaries = boundaries;
meshSolid.elements = elements;
meshSolid.vertices = vertices;
meshSolid.rings = rings;

load('flap_F.mat')
meshFluid.boundaries = boundaries;
meshFluid.elements = elements;
meshFluid.vertices = vertices;
meshFluid.rings = rings;

fem_F = {'P2', 'P1'};
fem_S = 'P2';
data_file_F = 'NS_data';
data_file_S = 'CSM_data';

t = [];
param = [];
vtk_filename = 'Figures/PerpFlap_';%[];
use_SUPG = false;
quad_order   = 4;
MESH.dim  = dim;

eval(data_file_F);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)
    
     eval(['DATA.Fluid.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
    
end
DATA.Fluid.param = param;

clear data_fields_name data

eval(data_file_S);
data_fields_name = fieldnames(data);

for i = 1 : length(data_fields_name)
    
    eval(['DATA.Solid.',data_fields_name{i}, '=', 'data.',data_fields_name{i},';']);
    
end
DATA.Solid.param = param;


[ MESH.Fluid, FE_SPACE_v, FE_SPACE_p] = build_F( dim, meshFluid.elements, meshFluid.vertices, meshFluid.boundaries, fem_F{1}, quad_order, DATA.Fluid, 'CFD', meshFluid.rings );

[ MESH.Solid, FE_SPACE_s ] = build_S( dim, meshSolid.elements, meshSolid.vertices, meshSolid.boundaries, fem_S, quad_order, DATA.Solid, 'CSM', meshSolid.rings );

[ FE_SPACE_g] = build_G( MESH.Fluid, fem_F{1}, dim, quad_order );%

MESH.Fluid.internal_dof_c{MESH.dim+1} = 1:FE_SPACE_p.numDof;

%% Generates mappings from solid to fluid interface dofs and viceversa
[MESH] = FSI_InterfaceMap(DATA, MESH);

for k = 1 : dim
    MESH.ndof_interface{k} = length(MESH.Interface_FSmap{k});
end

%% Find interface indices in the internal numbering
tmp = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
for i = 1 : dim
    tmp([FE_SPACE_v.numDofScalar*(i-1)+MESH.Fluid.dof_interface{i}]) = 1;
end

% fluid interface DoFs wrt internal numbering
MESH.Fluid.Gamma     = find(tmp(MESH.Fluid.internal_dof));
MESH.Fluid.II        = setdiff(1:length(MESH.Fluid.internal_dof),MESH.Fluid.Gamma);
MESH.Fluid.II_global = setdiff(MESH.Fluid.internal_dof,MESH.Fluid.internal_dof(MESH.Fluid.Gamma));

tmp = zeros(FE_SPACE_s.numDof,1);
for i = 1 : dim
    tmp([FE_SPACE_s.numDofScalar*(i-1)+MESH.Solid.dof_interface{i}]) = 1;
end

MESH.Solid.Gamma     = find(tmp(MESH.Solid.internal_dof));
MESH.Solid.II        = setdiff(1:length(MESH.Solid.internal_dof),MESH.Solid.Gamma);
MESH.Solid.II_global = setdiff(MESH.Solid.internal_dof,MESH.Solid.internal_dof(MESH.Solid.Gamma));

MESH.internal_dof    = [MESH.Fluid.II MESH.Fluid.Gamma' ...
            length(MESH.Fluid.internal_dof)+1:(length(MESH.Fluid.internal_dof)+length(MESH.Solid.II))];
Fluid_ReferenceNodes = MESH.Fluid.nodes(1:dim,:);

%% Time Setting
BDF_orderF = DATA.Fluid.time.BDF_order;
t0        = DATA.Fluid.time.t0;
dt        = DATA.Fluid.time.dt;
tf        = DATA.Fluid.time.tf;
t         = DATA.Fluid.time.t0;

X_n  = zeros(length(MESH.Fluid.II) + length(MESH.Fluid.Gamma) + length(MESH.Solid.II),1);
k_t  = 0;

%% Initalize Fluid Time Advance
TimeAdvanceF = BDF_TimeAdvance( BDF_orderF );

% read initial condition
v0  = [];
for k = 1 : FE_SPACE_v.numComponents
    v0  = [v0; DATA.Fluid.u0{k}(  MESH.Fluid.nodes(1,:), MESH.Fluid.nodes(2,:), t0, param )'];
end
    
p0 = zeros(FE_SPACE_p.numDof,1);

u = [v0; p0];
X_n(1:length(MESH.Fluid.internal_dof)) = u(MESH.Fluid.internal_dof);

% export initial condition (if it's the case)
if ~isempty(vtk_filename)
    CFD_export_solution(MESH.dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), ...
        MESH.Fluid.vertices, MESH.Fluid.elements, MESH.Fluid.numNodes, [vtk_filename,'Fluid'], 0);
end

TimeAdvanceF.Initialize( v0 );
for bd = 2 : BDF_orderF
    TimeAdvanceF.Append( v0 );
end

%% Initalize Solid Time Advance
TimeAdvanceS = Newmark_TimeAdvance( DATA.Solid.time.beta, DATA.Solid.time.gamma, dt );

% read displacement and velocity initial condition
u0  = [];
du0 = [];
for k = 1 : FE_SPACE_s.numComponents
    u0  = [u0; DATA.Solid.u0{k}(  MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), t0, param )'];
    du0 = [du0; DATA.Solid.du0{k}( MESH.Solid.nodes(1,:), MESH.Solid.nodes(2,:), t0, param )'];
end

d2u0 = 0*du0;

% export initial condition (if it's the case)
if ~isempty(vtk_filename)
    CSM_export_solution(MESH.dim, u0, MESH.Solid.vertices, MESH.Solid.elements, MESH.Solid.numNodes, [vtk_filename,'Solid'], 0);
end

TimeAdvanceS.Initialize( u0, du0, d2u0 );
Coef_MassS = TimeAdvanceS.MassCoefficient( );

%% Initalize Geometry Time Advance
% time discretization of the ALE follows the fluid one
ALE_velocity = zeros(MESH.Fluid.numNodes*dim, 1);
d_Fn         = zeros(MESH.Fluid.numNodes*dim, 1);

% only used in the fully implicit case
TimeAdvanceG = BDF_TimeAdvance( BDF_orderF );
TimeAdvanceG.Initialize( d_Fn );
for bd = 2 : BDF_orderF
    TimeAdvanceG.Append( d_Fn );
end

% create mesh motion DATA structure: young and poisson coefficients for the
% solid extension have to be set in the solid datafile
DATA.Geometry                  = DATA.Solid;
DATA.Geometry.Material_Model   = 'SEMMT';
DATA.Geometry.Stiffening_power = 0.8;

% assemble matrix
MeshMotionAssembler = CSM_Assembler( MESH.Fluid, DATA.Geometry, FE_SPACE_g );
Solid_Extension.matrix  = MeshMotionAssembler.compute_jacobian( zeros(FE_SPACE_g.numDof, 1) );

% solid-extension internal DoFs
internal_dofs_HE = [];
for k = 1 : dim
    internal_dofs_HE       = [ internal_dofs_HE (k-1)*FE_SPACE_v.numDofScalar+setdiff(1:FE_SPACE_v.numDofScalar, ...
        [MESH.Fluid.dof_interface{k}; MESH.ALE_dirichlet{k}])];
end
Solid_Extension.internal_dofs = internal_dofs_HE;

% compute LU factorization and store it
[Solid_Extension.L , Solid_Extension.U,...
    Solid_Extension.perm , q ]   = lu(Solid_Extension.matrix(internal_dofs_HE,internal_dofs_HE), 'vector');
Solid_Extension.invp             = 0*q ;
Solid_Extension.invp(q)          = 1:length(q);

%% Interface transfer matrices: solid to fluid and viceversa

% empty matrices for zero matrix blocks in the monolithic system
Z_FS  = sparse(length(MESH.Fluid.II), length(MESH.Solid.II));
Z_SF  = sparse(length(MESH.Solid.II), length(MESH.Fluid.II));

% transfer matrix from solid to fluid
IdGamma_SF = [];
for k = 1 : dim
    IdGamma_SF_tmp = sparse(MESH.ndof_interface{k}, MESH.ndof_interface{k});
    IdGamma_SF_tmp(MESH.Interface_FSmap{k}, :) = speye(MESH.ndof_interface{k},MESH.ndof_interface{k});
    IdGamma_SF = blkdiag(IdGamma_SF, IdGamma_SF_tmp);
end

% transfer matrix from fluid to solid
IdGamma_FS = [];
for k = 1 : dim
    IdGamma_FS_tmp = sparse(MESH.ndof_interface{k}, MESH.ndof_interface{k});
    IdGamma_FS_tmp(MESH.Interface_SFmap{k}, :) = speye(MESH.ndof_interface{k},MESH.ndof_interface{k});
    IdGamma_FS = blkdiag(IdGamma_FS, IdGamma_FS_tmp);
end

%% Create Solid Assembler Object
SolidModel = CSM_Assembler( MESH.Solid, DATA.Solid, FE_SPACE_s );

% Assemble mass matrix
M_s    =  SolidModel.compute_mass();
M_s    =  M_s * DATA.Solid.Density;

A_s = SolidModel.compute_jacobian( zeros(FE_SPACE_s.numDof, 1) );

A_robin = SolidModel.assemble_ElasticRobinBC();

%% Initialize Linear Solver
tol        = DATA.Solid.NonLinearSolver.tol;
maxIter    = DATA.Solid.NonLinearSolver.maxit;

%% PreProcessing for Drag and Lift Computation

compute_AerodynamicForces = true;

if compute_AerodynamicForces
    AeroF_x(k_t+1)  = 0;
    AeroF_y(k_t+1)  = 0;
    AeroF_z(k_t+1)  = 0;
    dofs_drag    = [];
    
    for j = 1 : length(DATA.Fluid.Output.DragLift.flag)
        Dirichlet_side         = find(MESH.Fluid.boundaries(MESH.Fluid.bc_flag_row,:) == DATA.Fluid.Output.DragLift.flag(j));
        Dirichlet_side         = unique(Dirichlet_side);
        Dirichlet_dof          = MESH.Fluid.boundaries(1:MESH.Fluid.numBoundaryDof,Dirichlet_side);
        dofs_drag              = [dofs_drag; Dirichlet_dof(:)];
    end
    dofs_drag = unique(dofs_drag);
    
    fileDragLift = fopen(DATA.Fluid.Output.DragLift.filename, 'w+');
    fprintf(fileDragLift, 'Time          F_x          F_y          F_z');
    fprintf(fileDragLift, '\n%1.4e  %1.4e  %1.4e  %1.4e', t, AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
end

R_P = zeros(FE_SPACE_s.numDof,1);
J_P = sparse(FE_SPACE_s.numDof,FE_SPACE_s.numDof);

X_nk = X_n;

%% Time Loop
fprintf('\n **** Starting temporal loop ****\n');
while ( t < tf )
    
    iter_time = tic;
    
    t       = t   + dt;
    k_t     = k_t + 1;
    
    fprintf('\n=========================================================================')
    fprintf('\n==========  t0 = %2.4f  t = %2.4f  tf = %2.4f\n',t0,t,tf);
    
    v_BDF = TimeAdvanceF.RhsContribute( );
    u_BDF = [v_BDF; zeros(FE_SPACE_p.numDof,1)];
    alphaF = TimeAdvanceF.GetCoefficientDerivative();

    % get lifting for F-velocity and S_displacement
    [~, ~, v_D]                 = CFD_ApplyBC([], [], FE_SPACE_v, FE_SPACE_p, MESH.Fluid, DATA.Fluid, t);            
    [~, ~, DisplacementDir_np1] = CSM_ApplyBC([], [], FE_SPACE_s, MESH.Solid, DATA.Solid, t);
            
     Csi    = TimeAdvanceS.RhsContribute( );
     F_ext  = SolidModel.compute_volumetric_forces( t );
     F_S    = F_ext + M_s * Csi;
            
     % Get displacment d^alpha on the interface
     d_alpha = TimeAdvanceS.M_dU - dt * DATA.Solid.time.gamma * Csi ...
                + dt* (1-DATA.Solid.time.gamma)*TimeAdvanceS.M_d2U;
     F_L     = d_alpha(MESH.Solid.internal_dof(MESH.Solid.Gamma));
            
     alpha   = DATA.Solid.time.gamma / (dt * DATA.Solid.time.beta);
            
     % Prepare for Newton's method
     X_nk = X_n;
     dX   = 0*X_n;
            
     % fluid current iteration
     u_nk = 0*u;
     u_nk(MESH.Fluid.internal_dof)  = X_nk(1:length(MESH.Fluid.internal_dof));
     u_nk(MESH.Fluid.Dirichlet_dof) = v_D;
            
     % solid current iteration
     d_nk = TimeAdvanceS.M_U;
     d_nk(MESH.Solid.Dirichlet_dof) = DisplacementDir_np1;
            
     ALE_rhsBDF = TimeAdvanceG.RhsContribute( );
            
     k = 1;
     norm_k   = tol + 1;
     
            % Start Newton's method
            while( k <= maxIter && (norm_k > tol || isnan(norm_k) ) )
                
                % Assemble Solid Tangent matrix and internal forces vector
                dA = SolidModel.compute_jacobian(  d_nk  );
                GS = SolidModel.compute_internal_forces( d_nk );
                
                dG_STR    = Coef_MassS * M_s + dA + A_robin + J_P;
                G_S       = Coef_MassS * M_s * d_nk + GS + A_robin * d_nk + R_P + J_P * d_nk - F_S;
                
                % Apply Solid boundary conditions
                [dG_STR, G_S] = CSM_ApplyBC(dG_STR, -G_S, FE_SPACE_s, MESH.Solid, DATA.Solid, t, 1);
                
                % Update Fluid Matrices
                FluidModel = CFD_Assembler( MESH.Fluid, DATA.Fluid, FE_SPACE_v, FE_SPACE_p );
                
                F_gravity = FluidModel.compute_external_forces(t);

                [A_Stokes] = FluidModel.compute_Stokes_matrix();
                
                Mv = FluidModel.compute_mass_velocity();
                Mp = FluidModel.compute_mass_pressure();
                M  = blkdiag(DATA.Fluid.density * Mv, 0*Mp);
                                
                [C1, C2]    = FluidModel.compute_convective_matrix_ALE( u_nk(1:FE_SPACE_v.numDof),  ALE_velocity);
                
                F_NS = 1/dt * M * (alphaF*u_nk - u_BDF) + A_Stokes * u_nk + C1 * u_nk - F_gravity;
                C_NS = alphaF/dt * M + A_Stokes + C1 + C2;
                
                % Apply Fluid boundary conditions
                [dG_NS, G_NS]   =  CFD_ApplyBC(C_NS, -F_NS, FE_SPACE_v, FE_SPACE_p, MESH.Fluid, DATA.Fluid, t, 1, u);
                
                % Interface Solid Stiffness expressed in the fluid numbering
                S_GG     = IdGamma_SF * (dG_STR(MESH.Solid.Gamma,MESH.Solid.Gamma) * IdGamma_FS);
                
                % Interface/Internal Solid Stiffness expressed in fluid/solid numbering
                S_GI     = IdGamma_SF *  dG_STR(MESH.Solid.Gamma,MESH.Solid.II);
                F_L2     = F_L - IdGamma_FS*X_nk(MESH.Fluid.Gamma) + alpha*d_nk(MESH.Solid.internal_dof(MESH.Solid.Gamma));
                
                % Form Monolothic System
                dG_FSI    = [dG_NS(MESH.Fluid.II,MESH.Fluid.II)                      dG_NS(MESH.Fluid.II,MESH.Fluid.Gamma)                    Z_FS ;...
                    dG_NS(MESH.Fluid.Gamma,MESH.Fluid.II)      dG_NS(MESH.Fluid.Gamma,MESH.Fluid.Gamma)+1/alpha*S_GG               S_GI  ;...
                    Z_SF                                                    1/alpha*dG_STR(MESH.Solid.II,MESH.Solid.Gamma)*IdGamma_FS   dG_STR(MESH.Solid.II,MESH.Solid.II)   ];
                
                G_FSI    = [G_NS(MESH.Fluid.II); ...
                    G_NS(MESH.Fluid.Gamma)+IdGamma_SF*G_S(MESH.Solid.Gamma)+1/alpha*(IdGamma_SF*(dG_STR(MESH.Solid.Gamma,MESH.Solid.Gamma)*F_L2)); ...
                    G_S(MESH.Solid.II)+1/alpha*dG_STR(MESH.Solid.II,MESH.Solid.Gamma)*F_L2];
                
                dX(MESH.internal_dof) = dG_FSI \ G_FSI;
                
                % Update solution
                norm_k   = norm(dX)/norm(X_nk);
                X_nk     = X_nk + dX;
                
                % Update current velocity and pressure
                u_nk(MESH.Fluid.internal_dof)= X_nk(1:length(MESH.Fluid.internal_dof));
                u_nk(MESH.Fluid.Dirichlet_dof) = v_D;
                
                % Update current solid displacement
                d_nk(MESH.Solid.II_global)= X_nk(1+length(MESH.Fluid.internal_dof):end);
                d_nk(MESH.Solid.Dirichlet_dof)= DisplacementDir_np1;
                d_nk(MESH.Solid.internal_dof(MESH.Solid.Gamma))= 1/alpha * ( IdGamma_FS*X_nk(MESH.Fluid.Gamma) -  F_L);
                
                % Deform Fluid mesh by Solid-Extension Mesh Motion technique
                d_F = FSI_SolidExtension(MESH, d_nk, Solid_Extension);
                Fluid_def_nodes    = Fluid_ReferenceNodes + d_F;
                Fluid_def_vertices = Fluid_def_nodes(1:dim, 1:MESH.Fluid.numVertices);
                
                d_F = reshape(d_F',dim*MESH.Fluid.numNodes,1);
                
                % Compute Fluid mesh velocity by BDF formula
                ALE_velocity  =  1/dt * ( alphaF*d_F - ALE_rhsBDF );
                
                % Update Fluid MESH jacobian etc
                MESH.Fluid.vertices = Fluid_def_vertices;
                MESH.Fluid.nodes    = Fluid_def_nodes;
                [MESH.Fluid.jac, MESH.Fluid.invjac, MESH.Fluid.h] = geotrasf(dim, MESH.Fluid.vertices, MESH.Fluid.elements);
                
                fprintf('\n   Iteration  k= %d;  norm(dX)/norm(Xk) = %1.2e, normRES = %1.2e \n',k,full(norm_k), full(norm(G_FSI)))
                k = k + 1;
                
            end
            
    norm_n   = norm(X_nk - X_n) / norm(X_n);
    X_n      = X_nk;
    
        %% Export solid displacement on reference mesh
    Displacement_np1                           = zeros(FE_SPACE_s.numDof,1);
    Displacement_np1(MESH.Solid.II_global)     = X_n(1+length(MESH.Fluid.internal_dof):end);
    Displacement_np1(MESH.Solid.Dirichlet_dof) = DisplacementDir_np1;
    Displacement_np1(MESH.Solid.internal_dof(MESH.Solid.Gamma))   = 1/alpha * ( IdGamma_FS*X_n(MESH.Fluid.Gamma) -  F_L);
    
    % Export to VTK
    if ~isempty(vtk_filename)
        CSM_export_solution(MESH.dim, Displacement_np1, MESH.Solid.vertices, ...
            MESH.Solid.elements, MESH.Solid.numNodes, [vtk_filename, 'Solid'], k_t);
    end
    
    % update time advance
    TimeAdvanceS.Update( Displacement_np1 );
    
    % Update geometry time advance
    TimeAdvanceG.Append( d_F );
    d_Fn          =  d_F;
    
    %% Export Fluid velocity and pressure on deformed mesh
    u(MESH.Fluid.internal_dof)  = X_n(1:length(MESH.Fluid.internal_dof));
    u(MESH.Fluid.Dirichlet_dof) = v_D;
        
    % Export to VTK
    if ~isempty(vtk_filename)
        CFD_export_solution(dim, u(1:FE_SPACE_v.numDof), u(1+FE_SPACE_v.numDof:end), ...
            MESH.Fluid.vertices, MESH.Fluid.elements, MESH.Fluid.numNodes, [vtk_filename,'Fluid'], k_t);
    end
    
    % Update fluid time advance
    TimeAdvanceF.Append( u(1:FE_SPACE_v.numDof) );
    
        %% Compute Aerodynamic Forces
    if compute_AerodynamicForces
        
        C_NS = 0*C_NS;
        F_NS = -F_NS;
        
        Z = zeros(FE_SPACE_v.numDofScalar,1);
        Z(dofs_drag)= 1;
        
        W = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
        W(1:FE_SPACE_v.numDofScalar) = Z;
        AeroF_x(k_t+1) = DATA.Fluid.Output.DragLift.factor*(W'*(-C_NS*u + F_NS));
        
        W = zeros(FE_SPACE_v.numDof+FE_SPACE_p.numDof,1);
        W(FE_SPACE_v.numDofScalar+[1:FE_SPACE_v.numDofScalar])  = Z;
        AeroF_y(k_t+1)  = DATA.Fluid.Output.DragLift.factor*(W'*(-C_NS*u  + F_NS));
        
        AeroF_z(k_t+1) = 0.0;
        
        fprintf('\n *** F_x = %e, F_y = %e, F_z = %e *** \n',  AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
        fprintf(fileDragLift, '\n%1.4e  %1.4e  %1.4e  %1.4e', t, AeroF_x(k_t+1), AeroF_y(k_t+1), AeroF_z(k_t+1));
    end
    
    X = [u; Displacement_np1];
    
end
  
 