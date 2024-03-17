function nodes2 = Qua4toTri3(nodes1)
% Purpose: To convert QUA4 mesh to TRI3 mesh
% Input:
%      nodes1 - nodal connectivity matrix of QUA4 elements 
% Output:
%      nodes2 - nodal connectivity matrix for TRI3 elements

% Coded by :    Siva Srinivas Kolukula, PhD      
%               Indian Tsunami Early Warning Centre (ITEWC)
%               Advisory Services and Satellite Oceanography Group (ASG)
%               Indian National Centre for Ocean Information Services (INCOIS)
%               Hyderabad, INDIA
% E-mail   :    allwayzitzme@gmail.com                                        
% web-link :    https://sites.google.com/site/kolukulasivasrinivas/   
%% Input check 
if ~any(size(nodes1)==4)
    error('Input the nodal connectivity data for Qua4 elements')
end
nel1 = size(nodes1,1) ;          % total number of QUA4 elements 
%% Initilaize TRI3 data 
nel2 = 2*nel1 ;
nodes2 = zeros(nel2,3) ;
%
nodes2(1:2:end,:) = nodes1(:,1:3) ;
nodes2(2:2:end,:) = nodes1(:,[1 3 4]) ;
