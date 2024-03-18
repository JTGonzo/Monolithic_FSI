classdef CFD_Assembler < handle
    
    properties (GetAccess = public, SetAccess = protected)
        M_MESH;
        M_DATA;
        M_FE_SPACE_v;
        M_FE_SPACE_p;
        M_totSize;
        M_density;
        M_dynamic_viscosity;
        M_gravity;
    end
   
    methods
        
        %==========================================================================
        %% Constructor
        function obj = CFD_Assembler( MESH, DATA, FE_SPACE_v, FE_SPACE_p )
            
            obj.M_MESH        = MESH;
            obj.M_DATA        = DATA;
            obj.M_FE_SPACE_v  = FE_SPACE_v;
            obj.M_FE_SPACE_p  = FE_SPACE_p;
            obj.M_totSize     = FE_SPACE_v.numDof + FE_SPACE_p.numDof;
            obj = SetFluidParameters( obj );
            
            if isfield(obj.M_DATA, 'gravity')
                obj.M_gravity = obj.M_DATA.gravity;
            else
                obj.M_gravity = [0 0 0];
            end
            
        end
        
        %==========================================================================
        %% SetMaterialParameters
        function obj = SetFluidParameters( obj )
            
             if isfield(obj.M_DATA, 'density')
                 obj.M_density = obj.M_DATA.density;
             else
                 warning('\nCFD_ASSEMBLER class: density not provided, set to 1 by default \n')
                 obj.M_density = 1;
             end
             
             if isfield(obj.M_DATA, 'dynamic_viscosity')
                 obj.M_dynamic_viscosity = obj.M_DATA.dynamic_viscosity;
             else
                 warning('\nCFD_ASSEMBLER class: kinematic viscosity not provided\n set to 1 by default \n')
                 obj.M_dynamic_viscosity = 1;
             end
             
        end
        
        %==========================================================================
        %% Compute_external_forces
        function F_ext = compute_external_forces(obj, t)
            
            if nargin < 2 || isempty(t)
                t = [];
            end
            
            % Computations of all quadrature nodes in the elements
            coord_ref = obj.M_MESH.chi;

            x = zeros(obj.M_MESH.numElem,obj.M_FE_SPACE_v.numQuadNodes); y = x;
            for j = 1 : 3
                 i = obj.M_MESH.elements(j,:);
                 vtemp = obj.M_MESH.vertices(1,i);
                 x = x + vtemp'*coord_ref(j,:);
                 vtemp = obj.M_MESH.vertices(2,i);
                 y = y + vtemp'*coord_ref(j,:);
            end
                    
                    % Evaluation of external forces in the quadrature nodes
            for k = 1 : obj.M_MESH.dim
                f{k}  = obj.M_DATA.force{k}(x,y,t,obj.M_DATA.param);
            end
                                
            % C_OMP assembly, returns matrices in sparse vector format
            F_ext = [];
            for k = 1 : obj.M_MESH.dim
                
                [rowF, coefF] = CSM_assembler_ExtForces(f{k}, obj.M_MESH.elements, obj.M_FE_SPACE_v.numElemDof, ...
                    obj.M_FE_SPACE_v.quad_weights, obj.M_MESH.jac, obj.M_FE_SPACE_v.phi);
                
                % Build sparse matrix and vector
                F_ext    = [F_ext; sparse(rowF, 1, coefF, obj.M_MESH.numNodes, 1)];
                
            end
            
            F_ext = [F_ext; zeros(obj.M_FE_SPACE_p.numDof,1)];
        end
        
        %==========================================================================
        %% compute_Stokes_matrix
        function A = compute_Stokes_matrix(obj)
            
            % C_OMP assembly, returns matrices in sparse vector format
            [rowA, colA, coefA] = ...
                CFD_assembler_C_omp('Stokes', obj.M_dynamic_viscosity, obj.M_MESH.dim, obj.M_MESH.elements, ...
                obj.M_FE_SPACE_v.numElemDof, obj.M_FE_SPACE_p.numElemDof, obj.M_FE_SPACE_v.numDof, ...
                obj.M_FE_SPACE_v.quad_weights, obj.M_MESH.invjac, obj.M_MESH.jac, ...
                obj.M_FE_SPACE_v.phi, obj.M_FE_SPACE_v.dphi_ref, obj.M_FE_SPACE_p.phi);
            
            % Build sparse matrix
            A   = sparse(rowA, colA, coefA, obj.M_totSize, obj.M_totSize);

        end
                        
        %==========================================================================
        %% compute_convective_matrix_ALE
        function [C1, C2] = compute_convective_matrix_ALE(obj, U_h, ALE_velocity)
            
            if nargin < 2 || isempty(U_h)
                U_h = zeros(obj.M_totSize,1);
            end
            
            if nargin < 3 || isempty(ALE_velocity)
                ALE_velocity = zeros(obj.M_totSize,1);
            end
            
            convective_velocity = U_h - ALE_velocity;
            
            % C_OMP assembly, returns matrices in sparse vector format
            [rowA, colA, coefA, rowB, colB, coefB] = ...
                CFD_assembler_C_omp('convectiveALE', 1.0, obj.M_MESH.dim, obj.M_MESH.elements, ...
                obj.M_FE_SPACE_v.numElemDof, obj.M_FE_SPACE_v.numDof, ...
                obj.M_FE_SPACE_v.quad_weights, obj.M_MESH.invjac, obj.M_MESH.jac, ...
                obj.M_FE_SPACE_v.phi, obj.M_FE_SPACE_v.dphi_ref, U_h, convective_velocity);
            
            % Build sparse matrix
            C1   = sparse(rowA, colA, coefA, obj.M_totSize, obj.M_totSize);
            C2   = sparse(rowB, colB, coefB, obj.M_totSize, obj.M_totSize);

            C1   = obj.M_density * C1;
            C2   = obj.M_density * C2;
        end
        
        %==========================================================================
        %% compute_mass_velocity
        function [M] = compute_mass(obj, FE_SPACE)
            
            % C_OMP assembly, returns matrices in sparse vector format
            [rowM, colM, coefM] = Mass_assembler_C_omp(obj.M_MESH.dim, obj.M_MESH.elements, FE_SPACE.numElemDof, ...
                FE_SPACE.quad_weights, obj.M_MESH.jac, FE_SPACE.phi);
            
            % Build sparse matrix
            M_scalar   = sparse(rowM, colM, coefM, FE_SPACE.numDofScalar, FE_SPACE.numDofScalar);
            M          = [];
            for k = 1 : FE_SPACE.numComponents
                M = blkdiag(M, M_scalar);
            end
            
        end
        
        %==========================================================================
        %% compute_mass_velocity
        function [Mv] = compute_mass_velocity(obj)
            Mv = compute_mass(obj, obj.M_FE_SPACE_v);
        end
        
        %==========================================================================
        %% compute_mass_pressure
        function [Mp] = compute_mass_pressure(obj)
            Mp = compute_mass(obj, obj.M_FE_SPACE_p);
        end
               
    end
    
end