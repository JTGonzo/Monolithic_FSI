function dof=completeDof(dof,Elem)
% assume the bdry dof in dof is good
% bdry elements are treated as isoparametric finite element
% note that at the beginning, all the dofs on bdry is correct        

bdNode=unique([Elem.bdNeuN;Elem.bdDirN]); %it can have overlap
elem2dof=Elem.elem2dof;
if Elem.deg>=3, bdElem=Elem.bdElem; end

switch Elem.deg
    case 1
        return;
    case 2        
        % The local node ordering is 
        %  3
        %  5  4
        %  1  6  2
        temp=dof(bdNode,:); % store the dofs on bdry which are correct
        dof(elem2dof(:,4),:)=0.5*(dof(elem2dof(:,2),:)+dof(elem2dof(:,3),:));
        dof(elem2dof(:,5),:)=0.5*(dof(elem2dof(:,3),:)+dof(elem2dof(:,1),:));
        dof(elem2dof(:,6),:)=0.5*(dof(elem2dof(:,1),:)+dof(elem2dof(:,2),:));
        dof(bdNode,:)=temp;
        % the curved Dir boundary don't need to change
    case 3
        % The local ordering is
        %  3
        %  5   7  
        %  8   10   4
        %  1   6    9   2
        
        temp=dof(bdNode,:); % store the dofs on bdry which are correct
        dof(elem2dof(:,4),:)=(2/3)*dof(elem2dof(:,2),:)+(1/3)*dof(elem2dof(:,3),:);
        dof(elem2dof(:,5),:)=(2/3)*dof(elem2dof(:,3),:)+(1/3)*dof(elem2dof(:,1),:);
        dof(elem2dof(:,6),:)=(2/3)*dof(elem2dof(:,1),:)+(1/3)*dof(elem2dof(:,2),:);
        
        dof(elem2dof(:,7),:)=(1/3)*dof(elem2dof(:,2),:)+(2/3)*dof(elem2dof(:,3),:);
        dof(elem2dof(:,8),:)=(1/3)*dof(elem2dof(:,3),:)+(2/3)*dof(elem2dof(:,1),:);
        dof(elem2dof(:,9),:)=(1/3)*dof(elem2dof(:,1),:)+(2/3)*dof(elem2dof(:,2),:);

        dof(elem2dof(:,10),:)=0.5*(dof(elem2dof(:,8),:)+dof(elem2dof(:,4),:));

        dof(bdNode,:)=temp;

        % for the bdry element, only need to take care of the inner dof
        dof(elem2dof(bdElem,10),:)=0.5*dof(elem2dof(bdElem,8),:)+0.5*dof(elem2dof(bdElem,4),:);        
    case 4
        % The local ordering 
        %  3
        %  8   10  
        %  5   15   4
        %  11  13   14   7
        %  1   9    6    12  2
        temp=dof(bdNode,:); % store the dofs on bdry which are correct
        dof(elem2dof(:,4),:)=0.5*(dof(elem2dof(:,2),:)+dof(elem2dof(:,3),:));
        dof(elem2dof(:,5),:)=0.5*(dof(elem2dof(:,3),:)+dof(elem2dof(:,1),:));
        dof(elem2dof(:,6),:)=0.5*(dof(elem2dof(:,1),:)+dof(elem2dof(:,2),:));
        dof(elem2dof(:,7),:)=0.5*(dof(elem2dof(:,2),:)+dof(elem2dof(:,4),:));
        dof(elem2dof(:,8),:)=0.5*(dof(elem2dof(:,3),:)+dof(elem2dof(:,5),:));
        dof(elem2dof(:,9),:)=0.5*(dof(elem2dof(:,1),:)+dof(elem2dof(:,6),:));
        dof(elem2dof(:,10),:)=0.5*(dof(elem2dof(:,4),:)+dof(elem2dof(:,3),:));
        dof(elem2dof(:,11),:)=0.5*(dof(elem2dof(:,5),:)+dof(elem2dof(:,1),:));
        dof(elem2dof(:,12),:)=0.5*(dof(elem2dof(:,6),:)+dof(elem2dof(:,2),:));
        dof(elem2dof(:,13),:)=0.5*(dof(elem2dof(:,5),:)+dof(elem2dof(:,6),:));
        dof(elem2dof(:,14),:)=0.5*(dof(elem2dof(:,6),:)+dof(elem2dof(:,4),:));
        dof(elem2dof(:,15),:)=0.5*(dof(elem2dof(:,4),:)+dof(elem2dof(:,5),:));
        dof(bdNode,:)=temp;
        
        % for the bdry element, only need to take care of the inner dof
        dof(elem2dof(bdElem,13),:)=2/3*dof(elem2dof(bdElem,11),:)+1/3*dof(elem2dof(bdElem,7),:);
        dof(elem2dof(bdElem,14),:)=1/3*dof(elem2dof(bdElem,11),:)+2/3*dof(elem2dof(bdElem,7),:);
        dof(elem2dof(bdElem,15),:)=0.5*(dof(elem2dof(bdElem,5),:)+dof(elem2dof(bdElem,4),:));
    case 5
        % The local ordering 
        %  3
        %  5   13  
        %  8   18   10
        %  11  20   19   7
        %  14  16   21   17  4
        %  1    6    9   12  15  2 
        temp=dof(bdNode,:); % store the dofs on bdry which are correct
        dof(elem2dof(:,4),:)=(4/5)*dof(elem2dof(:,2),:)+(1/5)*dof(elem2dof(:,3),:);
        dof(elem2dof(:,5),:)=(4/5)*dof(elem2dof(:,3),:)+(1/5)*dof(elem2dof(:,1),:);
        dof(elem2dof(:,6),:)=(4/5)*dof(elem2dof(:,1),:)+(1/5)*dof(elem2dof(:,2),:);
        
        dof(elem2dof(:,7),:)=(3/5)*dof(elem2dof(:,2),:)+(2/5)*dof(elem2dof(:,3),:);
        dof(elem2dof(:,8),:)=(3/5)*dof(elem2dof(:,3),:)+(2/5)*dof(elem2dof(:,1),:);
        dof(elem2dof(:,9),:)=(3/5)*dof(elem2dof(:,1),:)+(2/5)*dof(elem2dof(:,2),:);
        
        dof(elem2dof(:,10),:)=(2/5)*dof(elem2dof(:,2),:)+(3/5)*dof(elem2dof(:,3),:);
        dof(elem2dof(:,11),:)=(2/5)*dof(elem2dof(:,3),:)+(3/5)*dof(elem2dof(:,1),:);
        dof(elem2dof(:,12),:)=(2/5)*dof(elem2dof(:,1),:)+(3/5)*dof(elem2dof(:,2),:);
        
        dof(elem2dof(:,13),:)=(1/5)*dof(elem2dof(:,2),:)+(4/5)*dof(elem2dof(:,3),:);
        dof(elem2dof(:,14),:)=(1/5)*dof(elem2dof(:,3),:)+(4/5)*dof(elem2dof(:,1),:);
        dof(elem2dof(:,15),:)=(1/5)*dof(elem2dof(:,1),:)+(4/5)*dof(elem2dof(:,2),:);
        
        dof(elem2dof(:,16),:)=0.5*(dof(elem2dof(:,11),:)+dof(elem2dof(:,9),:));
        dof(elem2dof(:,17),:)=0.5*(dof(elem2dof(:,12),:)+dof(elem2dof(:,7),:));
        dof(elem2dof(:,18),:)=0.5*(dof(elem2dof(:,10),:)+dof(elem2dof(:,8),:));

        dof(elem2dof(:,19),:)=0.5*(dof(elem2dof(:,15),:)+dof(elem2dof(:,5),:));
        dof(elem2dof(:,20),:)=0.5*(dof(elem2dof(:,13),:)+dof(elem2dof(:,6),:));
        dof(elem2dof(:,21),:)=0.5*(dof(elem2dof(:,14),:)+dof(elem2dof(:,4),:));
        dof(bdNode,:)=temp;

        % for the bdry element, only need to take care of the inner dof
        dof(elem2dof(bdElem,16),:)=3/4*dof(elem2dof(bdElem,14),:)+1/4*dof(elem2dof(bdElem,4),:);
        dof(elem2dof(bdElem,17),:)=1/4*dof(elem2dof(bdElem,14),:)+3/4*dof(elem2dof(bdElem,4),:);
        dof(elem2dof(bdElem,18),:)=0.5*(dof(elem2dof(bdElem,8),:)+dof(elem2dof(bdElem,10),:));

        dof(elem2dof(bdElem,19),:)=1/3*dof(elem2dof(bdElem,11),:)+2/3*dof(elem2dof(bdElem,7),:);
        dof(elem2dof(bdElem,20),:)=2/3*dof(elem2dof(bdElem,11),:)+1/3*dof(elem2dof(bdElem,7),:);
        dof(elem2dof(bdElem,21),:)=0.5*(dof(elem2dof(bdElem,14),:)+dof(elem2dof(bdElem,4),:));

        
    otherwise
        disp('not implemented yet');
        pause;
end