classdef CSM_Assembler < handle
    
    properties (GetAccess = public, SetAccess = protected)
        M_MESH;
        M_DATA;
        M_FE_SPACE;
        M_subdomain;
        M_MaterialModel;
        M_MaterialParam;
    end
   
    methods
        
        %==========================================================================
        %% Constructor
        function obj = CSM_Assembler( MESH, DATA, FE_SPACE )
            
            obj.M_MESH      = MESH;
            obj.M_DATA      = DATA;
            obj.M_FE_SPACE  = FE_SPACE;
            obj.M_MaterialModel = DATA.Material_Model;
            obj = SetMaterialParameters(obj);
                       
        end
        
        %==========================================================================
        %% SetMaterialParameters
        function obj = SetMaterialParameters( obj )
            
            switch obj.M_MaterialModel
                case {'Linear', 'StVenantKirchhoff', 'NeoHookean'}
                    obj.M_MaterialParam = [obj.M_DATA.Young obj.M_DATA.Poisson];
                    
                case {'RaghavanVorp'}
                    obj.M_MaterialParam = [obj.M_DATA.Alpha obj.M_DATA.Beta obj.M_DATA.Bulk];
                     
                case 'SEMMT'
                    obj.M_MaterialParam = [obj.M_DATA.Young obj.M_DATA.Poisson obj.M_DATA.Stiffening_power];
            end
            
        end
        
        %==========================================================================
        %% Compute Volumetric Forces
        function F_ext = compute_volumetric_forces( obj, t )
            
            if nargin < 2 || isempty(t)
                t = [];
            end
            
            % Computations of all quadrature nodes in the elements
            coord_ref = obj.M_MESH.chi;
                    
            x = zeros(obj.M_MESH.numElem,obj.M_FE_SPACE.numQuadNodes); y = x;
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
                
                [rowF, coefF] = CSM_assembler_ExtForces(f{k}, obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                    obj.M_FE_SPACE.quad_weights, obj.M_MESH.jac, obj.M_FE_SPACE.phi);
                
                % Build sparse matrix and vector
                F_ext    = [F_ext; sparse(rowF, 1, coefF, obj.M_MESH.numNodes, 1)];
                
            end
            
        end
        
        %==========================================================================
        %% Compute mass matrix
        function [M] = compute_mass( obj )
            
            % C_OMP assembly, returns matrices in sparse vector format
            [rowM, colM, coefM] = Mass_assembler_C_omp(obj.M_MESH.dim, obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                obj.M_FE_SPACE.quad_weights, obj.M_MESH.jac, obj.M_FE_SPACE.phi);
            
            % Build sparse matrix
            M_scalar   = sparse(rowM, colM, coefM, obj.M_MESH.numNodes, obj.M_MESH.numNodes);
            M          = [];
            for k = 1 : obj.M_FE_SPACE.numComponents
                M = blkdiag(M, M_scalar);
            end
            
        end
                
        %==========================================================================
        %% Compute internal forces
        function [F_in] = compute_internal_forces(obj, U_h)
            
            % C_OMP assembly, returns matrices in sparse vector format
            [rowG, coefG] = ...
                CSM_assembler_C_omp(obj.M_MESH.dim, [obj.M_MaterialModel,'_forces'], obj.M_MaterialParam, full( U_h ), ...
                obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                obj.M_FE_SPACE.quad_weights, obj.M_MESH.invjac, obj.M_MESH.jac, obj.M_FE_SPACE.phi, obj.M_FE_SPACE.dphi_ref);
            
            % Build sparse matrix and vector
            F_in    = sparse(rowG, 1, coefG, obj.M_MESH.numNodes*obj.M_MESH.dim, 1);
            
        end
        
        %==========================================================================
        %% Compute internal forces Jacobian
        function [dF_in] = compute_jacobian(obj, U_h)
            
            if nargin < 2 || isempty(U_h)
                U_h = zeros(obj.M_MESH.dim*obj.M_MESH.numNodes,1);
            end

            % C_OMP assembly, returns matrices in sparse vector format
            [rowdG, coldG, coefdG] = ...
                CSM_assembler_C_omp(obj.M_MESH.dim, [obj.M_MaterialModel,'_jacobian'], obj.M_MaterialParam, full( U_h ), ...
                obj.M_MESH.elements, obj.M_FE_SPACE.numElemDof, ...
                obj.M_FE_SPACE.quad_weights, obj.M_MESH.invjac, obj.M_MESH.jac, obj.M_FE_SPACE.phi, obj.M_FE_SPACE.dphi_ref);
            
            % Build sparse matrix and vector
            dF_in   = sparse(rowdG, coldG, coefdG, obj.M_MESH.numNodes*obj.M_MESH.dim, obj.M_MESH.numNodes*obj.M_MESH.dim);
        end
        
        %==========================================================================
        %% Assemble Robin Condition: Pn + K d = 0 on \Gamma_Robin, with K = ElasticCoefRobin
        function [A] = assemble_ElasticRobinBC(obj)
            
            A = sparse(obj.M_MESH.numNodes*obj.M_MESH.dim, obj.M_MESH.numNodes*obj.M_MESH.dim);
            
        end
        
    end
    
end