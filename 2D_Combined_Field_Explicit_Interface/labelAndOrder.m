function [elem,hmin,hmax] = labelAndOrder(node,elem)
%**************************************************************************
% make it clockwisely oriented
%--------------------------------------------------------------------------
p1p2 = node(elem(:,2),:)-node(elem(:,1),:); % edge 12
p2p3 = node(elem(:,3),:)-node(elem(:,2),:); % edge 23
p3p1 = node(elem(:,1),:)-node(elem(:,3),:); % edge 31
detT = p1p2(:,2).*p3p1(:,1)-p1p2(:,1).*p3p1(:,2);
% equal to p1p2(:,1).*p2p3(:,2)-p1p2(:,2).*p2p3(:,1);
elem((detT<0),:) = elem((detT<0),[1,3,2]);
%--------------------------------------------------------------------------
% Compute length of each edge
%--------------------------------------------------------------------------
edgeLength(:,1) = sum(p2p3.^2,2);
edgeLength(:,2) = sum(p3p1.^2,2);
edgeLength(:,3) = sum(p1p2.^2,2);
hmin=sqrt(min(min(edgeLength)));
hmax=sqrt(max(max(edgeLength)));
%--------------------------------------------------------------------------
% Switch indices according the edge length
%--------------------------------------------------------------------------
[temp,I] = max(edgeLength,[],2);
elem((I==2),[1 2 3]) = elem((I==2), [2 3 1]);
elem((I==3),[1 2 3]) = elem((I==3), [3 1 2]);