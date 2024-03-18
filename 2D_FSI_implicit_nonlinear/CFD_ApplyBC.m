function [A_in, F_in, u_D] =  CFD_ApplyBC(A, F, FE_SPACE, FE_SPACE_p, MESH, DATA, t, zero_Dirichlet, U_k)

if nargin < 7
    t = [];
end

if isempty(A)
    A = sparse(FE_SPACE.numDof + FE_SPACE_p.numDof, FE_SPACE.numDof + FE_SPACE_p.numDof);
end

if isempty(F)
    F = sparse(FE_SPACE.numDof + FE_SPACE_p.numDof, 1);
end

if nargin < 8
    zero_Dirichlet = 0;
end


if nargin < 9
    U_k = zeros(FE_SPACE.numDof + FE_SPACE_p.numDof,1);
end

param = DATA.param;

u_D = [];
        
        %% Dirichlet condition
        for k = 1 : 2
            x = MESH.nodes(1,MESH.Dirichlet_dof_c{k});
            y = MESH.nodes(2,MESH.Dirichlet_dof_c{k});
            u_Dirichlet{k} = DATA.bcDir{k}(x,y,t,param);
            u_D = [u_D; u_Dirichlet{k}'];
        end
               
u_D  = u_D * (1 - zero_Dirichlet);

F_in = F(MESH.internal_dof) - A(MESH.internal_dof,MESH.Dirichlet_dof)*u_D;

A_in = A(MESH.internal_dof,MESH.internal_dof);

end
