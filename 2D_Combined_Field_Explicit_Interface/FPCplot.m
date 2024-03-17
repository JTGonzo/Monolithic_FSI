function FPCplot(DeEF,ElemF,ElemS,Sol,p,frame,nlev,str)
nv=ElemF.nVertex;
% showsolution(DeEF.dof(1:nv,:),ElemF.elem2dof(:,1:3),p(1:nv));
% trisurf(ElemF.elem2dof(:,1:3),DeEF.dof(1:nv,1),DeEF.dof(1:nv,2),p(1:nv),...
%     'FaceColor', 'interp', 'EdgeColor', 'interp');
% view(2);

lev=linspace(min(p),max(p),nlev);
a=min(abs(lev));
lev=lev(abs(lev)>2*a);

% lev1=linspace(min(p),min(p)/20,50);
% lev2=linspace(max(p)/20,max(p),50);
tricontour(DeEF.dof(1:nv,:),ElemF.elem2dof(:,1:3),p(1:nv),lev);


% for i=1:length(ElemF.barBDInd)
%     a=DeEF.dof(ElemF.barBDInd{i},:);
%     hold on; fill(a(:,1),a(:,2),'r');
% end
axis equal;
axis(frame);
text(1.3,0.2,num2str(Sol.t,4));
% title([str,',  t=',num2str(Sol.t,4)]);
if ~isempty(str)
    colorbar;
end

% figure(2);
% X=DeEF.dof(:,1);
% Y=DeEF.dof(:,2);
% F=p;
% LIN=ElemF.barBDInd;
% gx=unique(X);
% gy=unique(Y)';
% gx=linspace(0,2.5,3000);
% gy=linspace(0,0.41,800)';
% [Xr,Yr,Fr] = griddata(X,Y,F,gx,gy); %note the transpose at the end
% IN = inpolygon(Xr,Yr,X(LIN),Y(LIN));
% Fr(IN) = NaN;
% contourf(Xr,Yr,Fr);
% %contourf(Xr,Yr,Fr,200)
% shading flat
% colorbar

% trisurf(TRI,X,Y,F,'FaceColor','interp','EdgeColor','interp')
% view(2)
% colorbar