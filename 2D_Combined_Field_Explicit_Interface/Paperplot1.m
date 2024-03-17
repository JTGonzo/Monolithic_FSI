meshString={'FPCf2';'FPCs2'};
plotFrame=[0,1.5,0,0.41];


k=1;
for i=[1,3,4,5,6,8]
    str=['flelueGCLT',num2str(i),meshString{1}(1:4),'.mat'];
    load(str);
    subplot(2,3,k);
    k=k+1;

    FPCplot(DeEF,ElemF,ElemS,Sol,vort,plotFrame,[]);
    %             FPCplot(DeEF{1},ElemF,ElemS,Sol,Sol.p,'Pressure');
    colorbar off;
    drawnow;

end

load(str)

figure;
trisurf(ElemF.elem2dof(:,1:3),ElemF.dof(:,1),ElemF.dof(:,2)+0.3,zeros(size(ElemF.dof,1),1))
hold on;
trisurf(ElemS.elem2dof(:,1:3),ElemS.dof(:,1),ElemS.dof(:,2),zeros(size(ElemS.dof,1),1))
view(2); axis equal; axis off;

figure;
dof=DeEF.dof;
trisurf(ElemF.elem2dof(:,1:3),dof(:,1),dof(:,2)+0.3,zeros(size(dof,1),1))
hold on;
trisurf(ElemS.elem2dof(:,1:3),ElemS.dof(:,1),ElemS.dof(:,2),zeros(size(ElemS.dof,1),1))
view(2); axis equal; axis off;