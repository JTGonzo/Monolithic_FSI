function AnalyseData()

fileName = 'Data';
time = 0.002;
for j=1:81
    str = [fileName num2str(j) '.mat'];
    load(str);
    load('SolverInputs.mat','DeEF','ElemF','ElemS','Sol');

%     dof = cart2polar(ElemS.dof-repmat([0,4.983],size(ElemS.dof,1),1));
    dof = ElemS.dof;
    point = find(dof(:,1)==1 & dof(:,2)==0);

    for i=1:50
        pt = Data(:,:,i);
        pt = pt(ElemS.itfNode==point,:)+ElemS.dof(point,:);
        plot(time,pt(:,2),'+r');
        hold on;
        time = time+Sol.delt;
    end
end
hold off;

find(DeEF.detT<0)


function [x,y]= polar2cart (mag, ang_in_deg)
 x = mag .* cos(ang_in_deg*pi/180);
 y = mag .* sin(ang_in_deg*pi/180);
  
 
function pt = cart2polar(x)
 j = sqrt(-1);
 x = x(:,1)+j*x(:,2);
 r = abs(x);
 ar = angle(x);
 ad = ar*180/pi;
 pt = [r ad];
