function [node,elem,vorticity] = refinement(node,elem,vorticity)
%**************************************************************************
% Construct data structure
%--------------------------------------------------------------------------
totalEdge = [elem(:,[2,3]); elem(:,[3,1]); elem(:,[1,2])];
totalEdge = sort(totalEdge,2);
[edge, i1, ij] = unique(totalEdge,'rows','first'); 
% edge=totalEdge(i1,:), totalEdge=edge(ij,:);
% i1 is the last time it appears
N = size(node,1); NT = size(elem,1); NS = size(edge,1);
% ii(ij(3*NT:-1:1)) = 3*NT:-1:1; ii=ii'; % totalEdge(ii,:)-edge
% ii is the first time it appears
% use i1 and ii to determine the triangle the edge belongs to
% t1 and t2 are the triangle the edge belongs to
% k1 = ceil(i1/NT); t1 = i1 - NT*(k1-1);
% k2 = ceil(ii/NT); t2 = ii - NT*(k2-1);
% ix = (i1 ~= ii);
% edgeIx = 1:NS;
% elem2edge=accumarray([[t1(ix),k1(ix)];[t2,k2]],[edgeIx(ix),edgeIx],[NT 3]);

elem2edge = reshape(ij,NT,3);

%----------------------------------------------------------------------
% Add new nodes
%----------------------------------------------------------------------
% new node additions
node(N+1:N+NS,:) = (node(edge(:,1),:)+node(edge(:,2),:))/2;
vorticity(N+1:N+NS,:)= (vorticity(edge(:,1),:)+vorticity(edge(:,2),:))/2;

edge2newNode = (N+1:N+NS)';
%------------------------------------------------------------------
% refine each triangle into four triangles by regular refinement
%     3
%    / \
%   5 - 4
%  / \ / \
% 1 - 6 - 2
%------------------------------------------------------------------

% for creating the new element one need is the details of the edges which
% are present in a element.
t=1:NT;
p(t,1:3) = elem(t,1:3);
p(t,4:6) = edge2newNode(elem2edge(t,1:3));
elem(t,:) = p(t,[1,6,5]); %[p(t,1), p(t,6), p(t,5)];
elem(NT+1:2*NT,:) = p(t,[2,4,6]); %[p(t,2), p(t,4), p(t,6)];
elem(2*NT+1:3*NT,:) = p(t,[3,5,4]); %[p(t,3), p(t,5), p(t,4)];
elem(3*NT+1:4*NT,:) = p(t,4:6); %[p(t,4), p(t,5), p(t,6)];
% % Update boundary edges
% %--------------------------------------------------------------------------
% N = max(newNode(:,1)); if (N==0), return; end
% node2newNode = sparse(newNode(:,[2,3]),newNode(:,[3,2]),...
%     [newNode(:,1),newNode(:,1)],N,N);
% newgam=zeros(length(gam)*2,2);
% newgamI=gamI; 
% 
% % [gamI.N;gamI.D] has been ordered successively
% a=[gamI.N;gamI.D];
% s1=size(gamI.N,1); s2=size(gamI.D,1);
% b=zeros(size(a));
% b(1,1)=1;
% for i=1:size(a,1)
%     ind=a(i,1):a(i,2);
%     newgamAdd=refinebd(gam(ind,:),node2newNode);
%     if i>=2
%         b(i,1)=b(i-1,2)+1;
%     end
%     b(i,2)=b(i,1)+size(newgamAdd,1)-1;
%     newgam(b(i,1):b(i,2),:)=newgamAdd;
% end
% if s1>=1, newgamI.N=b(1:s1,:); end
% if s2>=1, newgamI.D=b(s1+1:s1+s2,:); end
% 
% % re-map when there is a curved bdry
% % gamI.i stands for interface, 
% % gamI.c stands for curved bdry. (odering in [bdN,bdD])
% if ~isempty(gamI.c)
%     a=[newgamI.N;newgamI.D];
%     for k=1:length(newgamI.c)
%         bd=a(newgamI.c(k),1):a(newgamI.c(k),2);
%         for i=1:length(bd)
%             p=node(newgam(bd(i),:),:);
%             node(newgam(bd(i),:),:)=projEllipse(p,newgamI.ellp(k,:));
%         end
%     end
% end