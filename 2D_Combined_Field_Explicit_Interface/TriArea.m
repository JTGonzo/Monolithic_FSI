function detT=TriArea(node,elem)
p1p2 = node(elem(:,2),:)-node(elem(:,1),:); % edge 12
p3p1 = node(elem(:,1),:)-node(elem(:,3),:); % edge 31
% after reordering the nodes on each triangle, it is positive
detT = p1p2(:,2).*p3p1(:,1)-p1p2(:,1).*p3p1(:,2); % det only. No abs
detT=abs(detT);
