function Elem=fixedMeshData(deg,meshString,nGlobalRefine)
% The result contans
% [gamI,elem2dof,dof,bdNeuElem,bdDirElem,inElem,itfElem,hm,...
%  nElem,nDof,nVertex,NeuBDN,DirBDN,freeN,nFreeN,nDirBDN,p1p2,p3p1,detT]

% gam was used to generate bdNeuElem, bdDirElem and itfElem,
% other than that, gam is not used in other places.
% gamI contain gamI.D and gamI.N are the index of the Dir and Neu edge
% We assume there is no element having two boundary edges


switch deg
    case 1
        indLBN=[3,2];
        cindLBN = [1 2 3];
    case 2
        indLBN=[3,4,2];
        cindLBN = [1 2 3 5 6];
    case 3
        indLBN=[3,7,4,2];
        cindLBN = [1 2 3 5 6 8 9];
    case 4
        indLBN=[3,10,4,7,2];
    case 5
        indLBN=[3,13,10,7,4,2];
    otherwise
        disp('not implemented yet');
        return;
end
Elem.deg=deg;
Elem.nnt=(deg+1)*(deg+2)/2; % is the number of nodes on each triangle
Elem.indLBN=indLBN;
Elem.cindLBN=cindLBN;
Elem.meshString=meshString;
Elem.nGlobalRefine=nGlobalRefine;
Elem.Quad=quadrature(deg); % assuming deg=sdeg
% Need to get a better idea about quadrature....

load(['./PerssonMesh/',meshString]);

gamI.A=[gamI.N;gamI.D]; % first N, then D.......gamI is the part of the meshstring which we have loaded...
% apart from gamI we have elem,gam,node
% gamI is a structure which again contains N, D

elem = double(elem); % this step has been added to include to handle the int32 format of the elem and sparse only handles double.

[elem,hm(1),hm(2)] = labelAndOrder(node,elem);

% % comment out the following 12 lines for the debugging
% To facilitate the refinement of boundary gam, 
% I will re-order the gam so that [gamI.N;gamI.D] are consecutively.
% for some historical reason, gamI.N goes first
% gamI.i is the index of the interface gam in the list of [gamI.N;gamI.D]

%gamI and gam have their role in mesh refinement process and interface elements...

a=[gamI.N;gamI.D] ; % N stands Neumman and D stands for Dirichlet
s1=size(gamI.N,1) ;
s2=size(gamI.D,1) ;
b=a(:,2)-a(:,1)+1 ;
b = cumsum(b); % cumsum means cummulative sum/addition
b =[[1;b(1:end-1)+1],b];
aa=[]; for i=1:size(a,1), aa=[aa,a(i,1):a(i,2)]; end

gam=gam(aa,:); % make gam in the order of the gamI.N and gamI.D

gamI.N=b(1:s1,:);
gamI.D=b(s1+1:s1+s2,:);

for i=1:nGlobalRefine
    [node,elem,gam,gamI] = uniformrefine(node,elem,gam,gamI);
end

hm=hm./(2^nGlobalRefine);

%--------------------------------------------------------------------------
% generate the index for the edge and element on or touching the boundary
%--------------------------------------------------------------------------

% the following is the same as the function initFE in refine_init
% we assume there is no element has two boundary edges
elem1=sort(elem(:,[2,3]),2);
elem2=sort(elem(:,[3,1]),2);
elem3=sort(elem(:,[1,2]),2);

orderedEdge=[elem1;elem2;elem3];
[edge,ii,ij]= unique(orderedEdge,'rows','first');
% [edge,jj]= unique(totalEdge,'rows','last');
% bdEdge=totalEdge(intersect(ii,jj),:);
% 

NN=size(node,1); NT=size(elem,1); NE=size(edge,1);

%% This part has been added on 29th dec 2012
%
%% Code starts Here
Node2Node = sparse(edge,[edge(:,2) edge(:,1)],ones(size(edge)),NN,NN);
Node2NodeN = sparse(edge,[edge(:,2) edge(:,1)],[edge(:,2) edge(:,1)],NN,NN);
%% Code Ends Here

% node2Edge=sparse(edge(:,1),edge(:,2),1:size(edge,1),NN,NN)+sparse(edge(:,2),edge(:,1),1:size(edge,1),NN,NN);
% 
% % Loop is the only way to access each element separately...
% elem2Edge=[];
% for p=1:size(elem,1)
%     elem2Edge = [elem2Edge;node2Edge(elem(p,2),elem(p,3)),node2Edge(elem(p,3),elem(p,1)),node2Edge(elem(p,1),elem(p,2)) ];
% end

elem2Edge = reshape(ij,NT,3);

% edge2Elem=[];
% 
% [tt,ss]=ismember(edge,elem1,'rows');
% [pp,qq]=find(tt==1);
% edge2Elem=[edge2Elem,ss];
% 
% [tt,ss]=ismember(edge,elem2,'rows');
% [pp,qq]=find(tt==1);
% edge2Elem=[edge2Elem,ss];
% 
% [tt,ss]=ismember(edge,elem3,'rows');
% [pp,qq]=find(tt==1);
% edge2Elem=[edge2Elem,ss];
% 
% edge2Elem = sort(edge2Elem,2,'descend');
% edge2Elem=edge2Elem(:,[1,2]);

% [row,col]=find((sum(edge2Elem(:,[1,2])>0,2)==1)==1);
% bdElem=edge2Elem(row,1);


% are all triangle

% the relative odering of interfaceElem in fluid and solid are the same

% bdEdge2Elem=edge2Elem(bdElem,1);

bdEdge = gam;
bdEdge = sort(bdEdge,2);
bdElem = find(sum(ismember(orderedEdge,bdEdge,'rows'),2)>=1);
bdElem(bdElem > NT & bdElem <=2*NT)=bdElem(bdElem > NT & bdElem <= 2*NT)-NT;
bdElem(bdElem>2*NT)=bdElem(bdElem>2*NT)-2*NT;
bdElem = unique(bdElem);

BdElem = sort(elem(bdElem,:),2);
BdEdges = [BdElem(:,[1 2]);BdElem(:,[2 3]);BdElem(:,[1 3])];
[var1,var2,BdEdges2BdElem] = intersect(bdEdge,BdEdges,'rows'); % var1 and var2 are not of importance
BdEdges2BdElem = [var2,BdEdges2BdElem];
BdEdges2BdElem = sortrows(BdEdges2BdElem,1);
BdEdges2BdElem = BdEdges2BdElem(:,2);
BdEdges2BdElem(BdEdges2BdElem > length(bdElem) & BdEdges2BdElem <=2*length(bdElem))...
    =BdEdges2BdElem(BdEdges2BdElem > length(bdElem) & BdEdges2BdElem <= 2*length(bdElem))-length(bdElem);
BdEdges2BdElem(BdEdges2BdElem >2*length(bdElem))=BdEdges2BdElem(BdEdges2BdElem>2*length(bdElem))-2*length(bdElem);
bdEdge2bdElem = BdEdges2BdElem;

inElem=setdiff((1:NT)',bdElem); % the index of the inner element 

% Nitf=zeros(length(gamI.i),1);
if ~isempty(gamI.i)
    a=[gamI.N;gamI.D];
    iedge=[];
    for j=1:length(gamI.i)
%         Nitf(j)=a(gamI.i(j),2)-a(gamI.i(j),1)+1; % # of edges for jth inf
        iedge=[iedge;sort(gam(a(gamI.i(j),1):a(gamI.i(j),2),:),2)]; % interface edge
    end
    itfelem = find(sum(ismember(orderedEdge,iedge,'rows'),2)==1);
    itfelem(itfelem > NT & itfelem <=2*NT)=itfelem(itfelem > NT & itfelem <= 2*NT)-NT;
    itfelem(itfelem>2*NT)=itfelem(itfelem>2*NT)-2*NT;
  
    [itfElem,ii,jj] = unique(itfelem);
    citfElem = itfelem(setdiff(1:size(itfelem),ii')); % citfElem defines the corner elements
else
    itfElem=[];
end

% temp1=cumsum(Nitf);
% temp2=cumsum([1;Nitf]);
% itfElemInfo=[temp2(1:end-1),temp1];

% icElem is the boundary element for the ball and the flag (contains
% both interface and curved boundary elements)
% icElem=itfElem;
% c_i=setdiff(gamI.c(:),gamI.i(:));

if ~isempty(gamI.c)
    a=[gamI.N;gamI.D];
    cedge=[];
    for j=1:length(gamI.c)
        cedge=[cedge;sort(gam(a(gamI.c(j),1):a(gamI.c(j),2),:),2)];
    end
    icElem = find(sum(ismember(orderedEdge,cedge,'rows'),2)+sum(ismember(orderedEdge,cedge(:,[2 1]),'rows'),2)==1);
    icElem(icElem > NT & icElem <=2*NT)=icElem(icElem > NT & icElem <= 2*NT)-NT;
    icElem(icElem>2*NT)=icElem(icElem>2*NT)-2*NT;
    icElem = unique(icElem);
else
    icElem = [];
end

%--------------------------------------------------------------------------
% add the edge index by N so that it's the index for the mid pts


%ix=(i1~=ii); %what is the role of these ix and these terms....
%edgeIx=1:NE;


% elem2edge(:,1) is the refinement edge
%elem2edge=accumarray([[t1(ix),k1(ix)];[t2,k2]],[edgeIx(ix),edgeIx],[NT,3]);

switch deg
    case 1
        % The local ordering is
        % 3
        % 1 2
        elem2dof=elem;
        dof=node;
    case 2        
        % The local node ordering is 
        %  3
        %  5  4
        %  1  6  2
        % The quadratic dof index in a element = edge index + N
        
        elem2dof=[elem,elem2Edge+NN];
        
        % the coordinates of the quadratic dof
        dof=[node;0.5*(node(edge(:,1),:)+node(edge(:,2),:))];
    case 3
        elem2dof=elem;
        elem2dof=[elem2dof,elem2Edge+NN,elem2Edge+(NE+NN)];
        % 1:NT is elem2elem
        elem2dof=[elem2dof,(1:NT)'+(2*NE+NN)];
        dof=node;
        % add points on the edge
        p1=node(edge(:,1),:);
        p2=node(edge(:,2),:);
        dof=[dof;2/3*p1+1/3*p2;1/3*p1+2/3*p2];
        % add points in the center
        p1=node(elem(:,1),:); p2=node(elem(:,2),:); p3=node(elem(:,3),:);
        dof=[dof;1/3*p1+1/3*p2+1/3*p3];

        % Now, elem2dof([1,2,3,10]) have the right index     
        % but not the rest. We might need to switch 
        % elem2dof(j,4) and elem2dof(j,7), 5,8 and 6,9. 
        % our goal is to get local ordering 
        %  3
        %  5   7  
        %  8   10   4
        %  1   6    9   2
        
        oldOrder=[4,7; 5,8; 6,9];
        newOrder=[7,4; 8,5; 9,6];
        for i=1:NT
            % r is the index of the pair that need to be switched
            r=antiParal(dof(elem2dof(i,[4,7,5,8,6,9]),:),...
                        dof(elem2dof(i,[2,3,3,1,1,2]),:));
            elem2dof(i,oldOrder(r,:))=elem2dof(i,newOrder(r,:));
            
%             temp1=elem2dof(i,:);
%             temp2=dof(temp1,:);
%             clf;
%             for j=1:10 
%                 plot(temp2(j,1),temp2(j,2),'o'); hold on; 
%                 text(temp2(j,1),temp2(j,2),num2str(j)); 
%             end 
%             pause();
        end

    case 4
        elem2dof=elem;
        elem2dof=[elem2dof,elem2Edge+NN,elem2Edge+(NE+NN),...
            elem2Edge+(2*NE+NN)];
        % 1:NT is elem2elem
        elem2dof=[elem2dof,(1:NT)'+(3*NE+NN),(1:NT)'+(NT+3*NE+NN),...
            (1:NT)'+(2*NT+3*NE+NN)];
        
        dof=node;
        % add points on the edge
        p1=node(edge(:,1),:);
        p2=node(edge(:,2),:);
        dof=[dof;0.5*p1+0.5*p2;0.75*p1+0.25*p2;0.25*p1+0.75*p2];
        % add points in the center
        p1=node(elem(:,1),:); p2=node(elem(:,2),:); p3=node(elem(:,3),:);
        dof=[dof;0.50*p1+0.25*p2+0.25*p3;...
                 0.50*p2+0.25*p3+0.25*p1;...
                 0.50*p3+0.25*p1+0.25*p2];

        % Now, elem2dof([1,2,3,4,5,6,13,14,15]) have the right index     
        % but not the rest. We might need to switch 
        % elem2dof(j,7) and elem2dof(j,10), 8,11 and 9,12. 
        % our goal is to get local ordering 
        %  3
        %  8   10  
        %  5   15   4
        %  11  13   14   7
        %  1   9    6    12  2
        
        oldOrder=[7,10;8,11;9,12];
        newOrder=[10,7;11,8;12,9];
        for i=1:NT
            % r is the index of the pair that need to be switched
            r=antiParal(dof(elem2dof(i,[7,10,8,11,9,12]),:),...
                        dof(elem2dof(i,[2,3,3,1,1,2]),:));
            elem2dof(i,oldOrder(r,:))=elem2dof(i,newOrder(r,:));
            
%             temp1=elem2dof(i,:);
%             temp2=dof(temp1,:);
%             clf;
%             for j=1:15 
%                 plot(temp2(j,1),temp2(j,2),'o'); hold on; 
%                 text(temp2(j,1),temp2(j,2),num2str(j)); 
%             end 
%             pause();
        end
        
    case 5
        elem2dof=elem;
        elem2dof=[elem2dof,elem2Edge+NN,elem2Edge+(NE+NN),...
            elem2Edge+(2*NE+NN),elem2Edge+(3*NE+NN)];
        % 1:NT is elem2elem
        elem2dof=[elem2dof,(1:NT)'+(4*NE+NN),(1:NT)'+(NT+4*NE+NN),...
            (1:NT)'+(2*NT+4*NE+NN),(1:NT)'+(3*NT+4*NE+NN),...
            (1:NT)'+(4*NT+4*NE+NN),(1:NT)'+(5*NT+4*NE+NN)];
        dof=node;
        % add points on the edge
        p1=node(edge(:,1),:);
        p2=node(edge(:,2),:);
        dof=[dof;0.8*p1+0.2*p2;0.6*p1+0.4*p2;0.4*p1+0.6*p2;0.2*p1+0.8*p2];
        % add points in the center
        p1=node(elem(:,1),:); p2=node(elem(:,2),:); p3=node(elem(:,3),:);
        dof=[dof;0.6*p1+0.2*p2+0.2*p3;...
                 0.2*p1+0.6*p2+0.2*p3;...
                 0.2*p1+0.2*p2+0.6*p3;...
                 0.2*p1+0.4*p2+0.4*p3;...
                 0.4*p1+0.2*p2+0.4*p3;...
                 0.4*p1+0.4*p2+0.2*p3];

        % Now, elem2dof([1,2,3,16:21]) have the right index     
        % but not the rest. We might need to switch 
        % elem2dof(j,[4,7]) and elem2dof(j,[13,10]), 
        % [5,8],[14,11] and [6,9],[15,12]. 
        % our goal is to get local ordering 
        %  3
        %  5   13  
        %  8   18   10
        %  11  20   19   7
        %  14  16   21   17  4
        %  1    6    9   12  15  2 
        
        oldOrder=[4,7,13,10; 5,8,14,11; 6,9,15,12];
        newOrder=[13,10,4,7; 14,11,5,8; 15,12,6,9];
        for i=1:NT
            % r is the index of the pair that need to be switched
            r=antiParal(dof(elem2dof(i,[7,10,8,11,9,12]),:),...
                        dof(elem2dof(i,[2,3,3,1,1,2]),:));
            elem2dof(i,oldOrder(r,:))=elem2dof(i,newOrder(r,:));
            
%             temp1=elem2dof(i,:);
%             temp2=dof(temp1,:);
%             clf;
%             for j=1:21 
%                 plot(temp2(j,1),temp2(j,2),'o'); hold on; 
%                 text(temp2(j,1),temp2(j,2),num2str(j)); 
%             end 
%             pause();

        end
    otherwise
        disp('not implemented yet');
        return;
end
%--------------------------------------------------------------------------
% rebuild the ordering of element touching the boundary so that the 2th 
% and 3rd pts are the pts on the boundary edge. 
% don't switch elem so that mesh retains the structure of Chen Long 
% Assuming no element has two boundary edges     
% when elem2dof is finished, for deg=4, elem2dof(bdEId,4,7,10) are the 
% points on the boundary edge. (bdEID,1) is the pt oppos. to it. 
% The rest are ordered counterclockwisely.
switch deg
    case 1
        reOrder=[1,2,3; 2,3,1; 3,1,2];
    case 2
        reOrder=[1,2,3,4,5,6;     2,3,1,5,6,4;     3,1,2,6,4,5];
        
    case 3
        reOrder=[1,2,3,4,5,6,7,8,9,10;...
                 2,3,1,5,6,4,8,9,7,10;...
                 3,1,2,6,4,5,9,7,8,10];
    case 4
        reOrder=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15;...
                 2,3,1,5,6,4,8,9,7,11,12,10,14,15,13;...
                 3,1,2,6,4,5,9,7,8,12,10,11,15,13,14];
    case 5
        reOrder=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21;...
                 2,3,1,5,6,4,8,9,7,11,12,10,14,15,13,17,18,16,20,21,19;...
                 3,1,2,6,4,5,9,7,8,12,10,11,15,13,14,18,16,17,21,19,20];
    otherwise
        disp('not available');
        return;
end

% % bdNode2bdEdge(i,j) is the index of edge ij or ji. i or j is index of node
% bdNode2bdEdge=sparse(edge(bdEdge,[1,2]),edge(bdEdge,[2,1]),...
%     repmat(bdEdge,1,2),NN,NN);

bdNeuElem=[]; bdNeuN=[]; a=[]; % NeuElem is the Neuman type bdry elem
for i=1:size(gamI.N,1) 
    a=[a,gamI.N(i,1):gamI.N(i,2)];
end
if ~isempty(a)
    [ii,bdNeuN,jj] = intersect(edge,sort(gam(a,:),2),'rows');
    if deg == 1
        if size(a,1)>1            
            bdNeuN = unique(gam(a,:));
        elseif size(a,1)==1
            bdNeuN = unique(gam(a,:)');
        end
    elseif deg ==2
        bdNeuN = bdNeuN+NN;
        if size(a,1)>1            
            bdNeuN = [bdNeuN;unique(gam(a,:))];
        elseif size(a,1)==1
            bdNeuN = [bdNeuN;unique(gam(a,:)')];
        end
    elseif deg ==3
        bdNeuN1 = bdNeuN+NN;
        bdNeuN2 = bdNeuN+NN+NE;
        if size(a,1)>1            
            bdNeuN = [bdNeuN1;bdNeuN2;unique(gam(a,:))];
        elseif size(a,1)==1
            bdNeuN = [bdNeuN1;bdNeuN2;unique(gam(a,:)')];
        end
    else
        disp('deg > 3 not yet implemnted');
        return;
    end
    bdNeuElem = find(sum(ismember(orderedEdge,gam(a,:),'rows'),2)+sum(ismember(orderedEdge,gam(a,[2 1]),'rows'),2)==1);
    bdNeuElem(bdNeuElem > NT & bdNeuElem <=2*NT)=bdNeuElem(bdNeuElem > NT & bdNeuElem <= 2*NT)-NT;
    bdNeuElem(bdNeuElem>2*NT)=bdNeuElem(bdNeuElem>2*NT)-2*NT;
    bdNeuElem = unique(bdNeuElem);
end

bdDirElem=[]; bdDirN = []; a=[]; % DirElem is the Dirichlet type bdry elem
for i=1:size(gamI.D,1) 
    a=[a,gamI.D(i,1):gamI.D(i,2)];
end
if ~isempty(a)
    [ii,bdDirN,jj] = intersect(edge,sort(gam(a,:),2),'rows');
    if deg ==1    
        if size(a,1)>1            
            bdDirN = unique(gam(a,:));
        elseif size(a,1)==1
            bdDirN = unique(gam(a,:)');
        end
    elseif deg == 2
        bdDirN = bdDirN+NN;
        if size(a,1)>1            
            bdDirN = [bdDirN;unique(gam(a,:))];
        elseif size(a,1)==1
            bdDirN = [bdDirN;unique(gam(a,:)')];
        end
    elseif deg ==3
        bdDirN1 = bdDirN+NN;
        bdDirN2 = bdDirN+NN+NE;
        if size(a,1)>1            
            bdDirN = [bdDirN1;bdDirN2;unique(gam(a,:))];
        elseif size(a,1)==1
            bdDirN = [bdDirN1;bdDirN2;unique(gam(a,:)')];
        end
    else
        disp('deg > 3 not yet implemnted');
        return;
    end
    bdDirElem = find(sum(ismember(orderedEdge,gam(a,:),'rows'),2)+sum(ismember(orderedEdge,gam(a,[2 1]),'rows'),2)==1);
    bdDirElem(bdDirElem > NT & bdDirElem <=2*NT)=bdDirElem(bdDirElem > NT & bdDirElem <= 2*NT)-NT;
    bdDirElem(bdDirElem>2*NT)=bdDirElem(bdDirElem>2*NT)-2*NT;
    bdDirElem = unique(bdDirElem);
end
    
for i=1:length(bdElem)
    
    eId=bdElem(i);
       
    [k,localInd]=setdiff(elem(eId,:),unique(bdEdge));
%     % check
%     a1=dof(elem2dof(eId,:),:);
%     for k=1:15, plot(a1(k,1)+sqrt(-1)*a1(k,2),'*'); hold on; pause(); end
    if isempty(k)
        [k,localInd]=setdiff(elem(eId,:),unique(iedge));
        if isempty(k)
            [ii,localInd,kk]=intersect([sort(elem(eId,[2 3]),2);sort(elem(eId,[1 3]),2);sort(elem(eId,[1 2]),2);],iedge,'rows');
            if length(localInd)==2               
                localInd=setdiff([1 2 3],localInd);
            end
        else
            [ii,localInd,kk]=intersect([sort(elem(eId,[2 3]),2);sort(elem(eId,[1 3]),2);sort(elem(eId,[1 2]),2);],iedge,'rows');
        end
        if isempty (localInd)
            [ii,localInd,kk]=intersect([sort(elem(eId,[2 3]),2);sort(elem(eId,[1 3]),2);sort(elem(eId,[1 2]),2);],bdEdge,'rows');
            if length(localInd)==2
                localInd=setdiff([1 2 3],localInd);
            end                
        end
    end
    if ~isempty(localInd)
        if (localInd>3 | localInd<1)
            pause;
        end
        if length(localInd)>1
            pause;
        end
        elem2dof(eId,:)=elem2dof(eId,reOrder(localInd,:)); % switch elem2dof
    end
%     % check
%     a1=dof(elem2dof(eId,:),:);
%     for k=1:15, plot(a1(k,1)+sqrt(-1)*a1(k,2),'*'); hold on;
%         text(a1(k,1),a1(k,2),num2str(k));
%     end
%     sum(a1([2,3,4,7,10],:).^2,2)-0.25
%     pause();

end

if ~isempty(gamI.c)
    % update the interface bdry node to make it on the curve.
    for k=1:length(gamI.c)
        ellp=gamI.ellp(k,:);
        b=[gamI.N;gamI.D];
        a=b(gamI.c(k),1):b(gamI.c(k),2);
        tempBDElem=bdElem(bdEdge2bdElem(a),:);

        p=zeros(3,2);
        for i=1:length(tempBDElem)
            nodeId=elem2dof(tempBDElem(i),:);
            p=dof(nodeId([2;3]),:); % pt 2 and pt 3
            theta=atan2((p(:,2)-ellp(2))./ellp(4), (p(:,1)-ellp(1))./ellp(3));
            theta0=theta(1); % pt 2
            dtheta=theta(2)-theta(1); % pt 3 - pt 2
            if dtheta>pi % in case pt2 and pt3 are on two sides of negative real axis
                theta(2)=theta(2)-2*pi;
                dtheta=theta(2)-theta(1); % pt 3 - pt 2
            elseif dtheta<-pi
                theta(1)=theta(1)-2*pi;
                dtheta=theta(2)-theta(1); % pt 3 - pt 2
            end

            switch deg
                case 1
                    
                case 2
                    theta=theta0+0.5*dtheta;
                    dof(nodeId(4),:)=ellp(1:2)...
                        +[ellp(3)*cos(theta),ellp(4)*sin(theta)];
                    
                case 3
                    theta=theta0+[1/3;2/3]*dtheta;
                    dof(nodeId([4,7]),:)=repmat(ellp(1:2),2,1)...
                        +[ellp(3)*cos(theta),ellp(4)*sin(theta)];
                    dof(nodeId(10),:)=0.5*dof(nodeId(8),:)+0.5*dof(nodeId(4),:);
                case 4
                    theta=theta0+[0.25;0.5;0.75]*dtheta;
                    dof(nodeId([7,4,10]),:)=repmat(ellp(1:2),3,1)...
                        +[ellp(3)*cos(theta),ellp(4)*sin(theta)];

                    dof(nodeId(13),:)=2/3*dof(nodeId(11),:)+1/3*dof(nodeId(7),:);
                    dof(nodeId(14),:)=1/3*dof(nodeId(11),:)+2/3*dof(nodeId(7),:);
                    dof(nodeId(15),:)=0.5*dof(nodeId(5),:)+0.5*dof(nodeId(4),:);
                case 5
                    theta=theta0+[0.2;0.4;0.6;0.8]*dtheta;
                    dof(nodeId([4,7,10,13]),:)=repmat(ellp(1:2),4,1)...
                        +[ellp(3)*cos(theta),ellp(4)*sin(theta)];

                    dof(nodeId(16),:)=3/4*dof(nodeId(14),:)+1/4*dof(nodeId(4),:);
                    dof(nodeId(17),:)=1/4*dof(nodeId(14),:)+3/4*dof(nodeId(4),:);
                    dof(nodeId(18),:)=0.5*dof(nodeId(8),:)+0.5*dof(nodeId(10),:);

                    dof(nodeId(19),:)=1/3*dof(nodeId(11),:)+2/3*dof(nodeId(7),:);
                    dof(nodeId(20),:)=2/3*dof(nodeId(11),:)+1/3*dof(nodeId(7),:);
                    dof(nodeId(21),:)=0.5*dof(nodeId(14),:)+0.5*dof(nodeId(4),:);
                otherwise
                    disp('not implemented');
            end
        end
    end
end

itfNode=unique(elem2dof(setdiff(itfElem,citfElem),Elem.indLBN)); 
citfNode=unique(elem2dof(citfElem,cindLBN)); 
if length(citfElem)==1
    itfNode = unique([itfNode;citfNode']);

else
    itfNode = unique([itfNode;citfNode]);
end

% showmesh(node,elem);
% for i=1:length(itfElem)
%     p=dof(elem2dof(itfElem(i),:),:);
%     hold on;
%     plot(p(:,1),p(:,2),'r*');
%     hold on;
%     for j=1:size(p,1)
%         hold on
%         text(p(j,1),p(j,2),num2str(j));
%     end
% end

% Elem.gamI=gamI;
Elem.hm=hm;
Elem.elem2dof=elem2dof;
Elem.dof=round(dof*10^(4))*10^(-4);
Elem.bdNeuElem=bdNeuElem;
Elem.bdDirElem=bdDirElem;
Elem.bdElem=[bdElem];
Elem.inElem=inElem;
Elem.itfElem=itfElem;
Elem.citfElem=citfElem;
% Elem.itfElemInfo=itfElemInfo;
Elem.icElem=icElem;
Elem.itfEdge = iedge;
Elem.itfNode = itfNode;

%% added on Dec 29 2012
Elem.Node2Node = Node2Node;
%% added on Dec 30 2012
Elem.Node2NodeN = Node2NodeN;

% Remark added on (June 24, 2011):
% please note that bdNeuElem also contains the interface element
% open boundary is only part of it.
Elem.openbdryElem=setdiff(Elem.bdNeuElem,Elem.itfElem);



Elem.bdNeuN= bdNeuN; % Neu include end points
            
Elem.bdDirN= bdDirN; % Dir include end points
Elem.freeN=setdiff((1:size(dof,1))',Elem.bdDirN); % excluding the Dir dofs

Elem.nDof=size(dof,1);
Elem.nVertex=size(node,1);
Elem.nElem=size(elem2dof,1);   
Elem.nBDElem=length(Elem.bdElem);
Elem.nFreeN=length(Elem.freeN);



% % p1p2 and p3p1 are calculated based on dof and elem2dof
% Elem.p1p2 = dof(elem2dof(:,2),:)-dof(elem2dof(:,1),:); % edge 12
% Elem.p3p1 = dof(elem2dof(:,1),:)-dof(elem2dof(:,3),:); % edge 31
% % after reordering the nodes on each triangle, it is positive
% Elem.detT = Elem.p1p2(:,2).*Elem.p3p1(:,1)-Elem.p1p2(:,1).*Elem.p3p1(:,2); % det only. No abs
% if min(Elem.detT)<0 | abs(length(Elem.bdElem)+length(inElem)-Elem.nElem)~=0
%     disp('det<0 or wrong elem in fixedMeshData'); pause; 
% end

%**************************************************************************
function [node,elem,newgam,newgamI] = uniformrefine(node,elem,gam,gamI)
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
% this column of the new node contains node numbers of the new nodes
newNode(:,1) = (N+1:N+NS)';
% This column contains the node numbers of the left corner node
newNode(:,2) = edge(:,1);
% this column contains the node numbers of the right corner node
newNode(:,3) = edge(:,2);
% this is a straight forward conversion or relation between the new node
% and the edge.
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
% Update boundary edges
%--------------------------------------------------------------------------
N = max(newNode(:,1)); if (N==0), return; end
node2newNode = sparse(newNode(:,[2,3]),newNode(:,[3,2]),...
    [newNode(:,1),newNode(:,1)],N,N);
newgam=zeros(length(gam)*2,2);
newgamI=gamI; 

% [gamI.N;gamI.D] has been ordered successively
a=[gamI.N;gamI.D];
s1=size(gamI.N,1); s2=size(gamI.D,1);
b=zeros(size(a));
b(1,1)=1;
for i=1:size(a,1)
    ind=a(i,1):a(i,2);
    newgamAdd=refinebd(gam(ind,:),node2newNode);
    if i>=2
        b(i,1)=b(i-1,2)+1;
    end
    b(i,2)=b(i,1)+size(newgamAdd,1)-1;
    newgam(b(i,1):b(i,2),:)=newgamAdd;
end
if s1>=1, newgamI.N=b(1:s1,:); end
if s2>=1, newgamI.D=b(s1+1:s1+s2,:); end

% re-map when there is a curved bdry
% gamI.i stands for interface, 
% gamI.c stands for curved bdry. (odering in [bdN,bdD])
if ~isempty(gamI.c)
    a=[newgamI.N;newgamI.D];
    for k=1:length(newgamI.c)
        bd=a(newgamI.c(k),1):a(newgamI.c(k),2);
        for i=1:length(bd)
            p=node(newgam(bd(i),:),:);
            node(newgam(bd(i),:),:)=projEllipse(p,newgamI.ellp(k,:));
        end
    end
end
%--------------------------------------------------------------------------
function newp=projEllipse(p,ellp) % project onto the ellipse
theta=atan2((p(:,2)-ellp(2))./ellp(4), (p(:,1)-ellp(1))./ellp(3));
newp=[ellp(3)*cos(theta)+ellp(1),ellp(4)*sin(theta)+ellp(2)];
%--------------------------------------------------------------------------
function bdEdge = refinebd(bdEdge,node2newNode)
if isempty(bdEdge), return; end
NB = size(bdEdge,1);
bdEdge(2*NB,:) = [0 0];
for k = 1:NB
    i = bdEdge(k,1);
    j = bdEdge(k,2);
    p = node2newNode(i,j); % node p is added to edge ij
    if p >0
        bdEdge(k,:) = [i,p];
        bdEdge(NB+k,:) = [p,j];
    end
end
bdEdge;
%**************************************************************************
function r=antiParal(L1,L2)
% L1(i,[1,2]) are the x y cord of ith points in L1. i=1:6
r=sum(((L1(1:2:5,:)-L1(2:2:6,:))./(L2(1:2:5,:)-L2(2:2:6,:)+1e-8)),2)<0;
% the pair that is anti-parallel
r=find(r>0);