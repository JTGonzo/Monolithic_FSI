clc;
clear all;

global GM ALE % Global matrix and ALE matrix
% u-eta formulation
% use the conservative formulation

tic
% If the first letter is capital, it is a structure.
% all input arg's are odered alphabatically. So is output arg's.
% Structures are always listed in the end

% suffix s is used to define the solid structure
% sffix f is used to define the fluid medium...
delt= [0.1];% 0.05 .025 .0125 .00625 1e-3];%0.1;0.05;0.025;0.0125;8e-3;4e-3;2e-3;1e-3;%5e-5;
titleat1 = {'eta01at1.mat' 'eta005at1.mat' 'eta0025at1.mat' 'eta00125at1.mat' 'eta000625at1.mat' 'etaActat1.mat'};
title = {'eta01.mat' 'eta005.mat' 'eta0025.mat' 'eta00125.mat' 'eta000625.mat' 'etaAct.mat'};
titleU = {'U01.mat' 'U005.mat' 'U0025.mat' 'U00125.mat' 'U000625.mat' 'UAct.mat'};
titleP = {'P01.mat' 'P005.mat' 'P0025.mat' 'P00125.mat' 'P000625.mat' 'PAct.mat'};
titleV = {'V01.mat' 'V005.mat' 'V0025.mat' 'V00125.mat' 'V000625.mat' 'VAct.mat'};


for loop = 1:length(delt)

ContinueProb = 0;

if ContinueProb
    % This block shoould contain the code which will load the necessary
    % data for continuation of a scheme.
else
    
    meshString={'FPCf_NonConservative';'FPCs_NonConservative'};

    if strcmp(meshString{1}(1:4),'FPCf')==1 % string comparison...
        Sol.delt = 0.002;
        Sol.tend=25; %8 for benchmark test, 0.4 for accuracy check
        Sol.movie=0;
        Sol.SVK=1; % =1; SVK with linear approximation (semi-implicit)

        % physical parameters and computational domain
        Phy.frho=1000;
        Phy.fmu=1e-3;

        Phy.srho=10000;
        Phy.smu=5e5;
        Phy.snu = 0.4;
        E = 2*(1+Phy.snu)*Phy.smu;
        lam = Phy.snu*E/((1+Phy.snu)*(1-2*Phy.snu));
        Phy.slam=2*Phy.snu*Phy.smu/(1-2*Phy.snu);%500;
        Phy.g=0;
        Phy.theta=0;
        Phy.Uo = 1;
        Phy.outbdryx=2.5;

        Sol.ad=zeros(floor(Sol.tend/Sol.delt),2);
        Sol.ldc=zeros(floor(Sol.tend/Sol.delt),2);
        nn=1;
        plotDt=Sol.delt;
        plotFrame=[0,1.5,0,0.41];
        
        Sol.Conservative = 0;
        Sol.Temam = 1;
        Sol.nGlobalRefine=0; 
        Sol.deg=2; 
        Sol.Order = 2;

    elseif (strcmp(meshString{1},'FS4f')==1)    
        Sol.delt = delt(loop);
        Sol.tend=10;
        Sol.movie=0;
        Sol.SVK=1; % =0; Linear elasticity

        Phy.frho=1;
        Phy.fmu=1;%1e-2;
        Phy.srho=1;
        PoisonR = 0.3;
        E = 100;
        
        Phy.smu = E/(2*(1+PoisonR));
        Phy.slam= E*PoisonR/((1+PoisonR)*(1-2*PoisonR));
        Phy.g=0;
        Phy.theta=pi/2;
        Phy.Uo = 1;
        
        plotDt=Sol.delt;
        plotFrame=[0,6.5,-0.5,1];

        Sol.nGlobalRefine=0;
        Sol.deg=2;  
        Sol.Order =2;

    elseif (strcmp(meshString{1},'BenchMarkF')==1)    
        Sol.delt=0.01;0.1;0.05;0.025;0.0125;8e-3;4e-3;2e-3;1e-3;%5e-5;
        Sol.tend=8;
        Sol.movie=1;
        Sol.SVK=0; % =0; Linear elasticity

        Phy.frho=1;
        Phy.fmu=10^(-3);%1e-2;
        Phy.srho=1;
        Phy.smu=10^50;
        Phy.slam=10^50;
        Phy.g=0;
        Phy.theta=0;
        Phy.Uo = 1.5;

        plotDt=Sol.delt;
        plotFrame=[0,2.2,0,0.41];

        Sol.nGlobalRefine=1;
        Sol.deg=2;  
        Sol.Order =1;

    elseif (strcmp(meshString{1},'FlexFlatPlateF')==1)...
            ||(strcmp(meshString{1},'FlexFlatPlateTh02F')==1)...
            ||(strcmp(meshString{1},'FlexFlatPlateTh005F')==1)       
        Sol.delt=0.001;0.1;0.05;0.025;0.0125;8e-3;4e-3;2e-3;1e-3;%5e-5;
        Sol.tend=15;
        Sol.movie=0;
        Sol.SVK=0; % =0; Linear elasticity
        
        Re = 1000;
        Kb = 0.0001;
        PoisonR = 0.3;
        h = 0.01; % plate thickness in length scale, this is computation model dependent
        L = 1; % Length of the plate in length scale computation model dependent
        Phy.frho=1;
        Phy.Uo = 1;
        Phy.fmu= Phy.Uo*L/Re;
        Phy.srho=2.5;
        Phy.thk = 0.01;
                
        E = Kb*12*(1-PoisonR^2)*Phy.frho*Phy.Uo^2*L^3/h^3;
        
        Phy.smu = E/(2*(1+PoisonR));
        Phy.slam= E*PoisonR/((1+PoisonR)*(1-2*PoisonR));
        Phy.g=0;
        Phy.theta=0;
        Phy.thk = 0.01;


        plotDt=Sol.delt;
        plotFrame=[-2,20.5,-5,5];

        Sol.nGlobalRefine=0;
        Sol.deg=2;
        Sol.MinAngle = 20; % in degrees
        Sol.state = 'deformed';
        Sol.Order =2;
        
    elseif (strcmp(meshString{1},'FlexCurvedPlateF')==1)    
        Sol.delt=0.002;0.1;0.05;0.025;0.0125;8e-3;4e-3;2e-3;1e-3;%5e-5;
        Sol.tend=10;
        Sol.movie=0;
        Sol.SVK=0; % =0; Linear elasticity

        Phy.frho=1;
        Phy.fmu=1e-3;1e-2;1*10^(-3);%1e-2;
        Phy.srho=10;
        Phy.smu=1680/4;7.617*10^10;
        Phy.slam=2520/4;9.695*10^10;
        Phy.g=0;
        Phy.theta=0;
        Phy.Uo = 1;

        plotDt=Sol.delt;
        plotFrame=[-2,10,-1,1];

        Sol.nGlobalRefine=1;
        Sol.deg=2;
        Sol.MinAngle = 20; % in degrees

    elseif strcmp(meshString{1},'HVf1')==1
        Sol.delt=0.2;
        Sol.tend=10;
        Sol.movie=0;
        Sol.SVK=0; % =0; Linear elasticity

        Phy.frho=1;
        Phy.fmu=1;%1e-2;
        Phy.srho=1;
        Phy.smu=100;
        Phy.slam=2000;
        Phy.g=0;
        Phy.theta=0;

        plotDt=Sol.delt;
        plotFrame=[0,6,-0.7,0.7];

        Sol.nGlobalRefine=0;
        Sol.deg=3;  
    end

    % Elem will store pure geometric and topologic inform. of mesh and node
    % DeE will store those data that will change when dof changes

    % this line is being modified and the original line is commented below.. In
    % this line DeES has no application in this function..

    % [ElemF,ElemS,Sol,DeES]=fixedDataue(meshString,Phy,Sol);

    [ElemF,ElemS,Sol]=fixedDataue(meshString,Phy,Sol);
    save('ElemS.mat','ElemS');
    save('ElemF.mat','ElemF');
%     ALE.dof = ElemF.dof;

    %--------------------------------------------------------------------------
    % initialize 
    Sol.t=0;
    Sol.me=zeros(ElemF.nDof,2,3); 

    % load('U.mat');
    % Sol.ActU = U;
    Sol.U=zeros(ElemF.nDof,2,2);
    Sol.MU =zeros(ElemF.nDof,2); % Mass matrix times u on the old mesh
    Sol.sMU =zeros(ElemS.nDof,2); % Mass matrix times u on the old mesh

    Sol.eta=zeros(ElemS.nDof,2); % displacement
    
    Sol.dof = ElemF.dof;
    if ((strcmp(meshString{1},'FlexFlatPlateF')==1)...
            ||(strcmp(meshString{1},'FlexFlatPlateTh02F')==1)...
            ||(strcmp(meshString{1},'FlexFlatPlateTh005F')==1))...
            && (strcmp(Sol.state,'deformed'))
        R = 4.983;
        no = unique(ElemS.dof(:,2));
        no = sort(no,'descend');
        if min(no)~= -Phy.thk/2 && max(no)~= Phy.thk/2
            disp('Error');
            return;
        end
        for i=1:length(no)
            r = R-no(i);
            Sol.eta(ElemS.dof(:,2)==no(i),:)=-ElemS.dof(ElemS.dof(:,2)==no(i),:)+r*[sin(ElemS.dof(ElemS.dof(:,2)==no(i),1)/r) (R/r-cos(ElemS.dof(ElemS.dof(:,2)==no(i),1)/r))];
        end
    end
    
    Sol.etat=zeros(ElemS.nDof,2,2); % velocity deta/dt
    Sol.etatt = zeros(ElemS.nDof,2,2); % accelaration 
    Sol.etaN=Sol.eta+Sol.delt*(Sol.etat(:,:,2)); % displacement
    Sol.sL = 0;
    
    Sol.error = [];
    DataU = [];
    DataEta = [];
    DataEtat = [];
    DataP = [];
    DataDeF = [];
    DataVort = [];


    % aviobj = avifile('example.avi','compression','none','fps',10,'quality',75);
    DeEF=deformedMeshData(ElemF.dof,ElemF); 
    moviecount=1;
    nn = 1;
    iterationStart = 1;

    end

    DataCount = 1;
    Counter = 1;
    
    for nit=iterationStart:round(Sol.tend/Sol.delt) % nit = number of iterations...
        t=nit*Sol.delt;

        Sol.t=t;
        Sol.nit = nit;     
        nit
        if Sol.Order ==2
            if nit ==1                
                Sol.fM=0;
                Sol.eta(:,:,2:3)=zeros(ElemS.nDof,2,2);
                Sol.MU(:,:,2)=Sol.MU(:,:,1);
                
                Sol.etat=zeros(ElemS.nDof,2,3); % velocity deta/dt
                Sol.etatt=zeros(ElemS.nDof,2,2); % accelaration
                Sol.sMU(:,:,1:3)=zeros(ElemS.nDof,2,3);
                Sol.gg = 0;
            end  
            
            % Select the solver either conservative or non conservative
            if Sol.Conservative ==1                
                [DeEF,Sol]=feue2SecondOrder(ALE,DeEF,ElemF,ElemS,Phy,Sol,meshString); %%% most important step and function might be..
            elseif Sol.Conservative ==0
               [DeEF,Sol] = feue2SecondOrderNonCons(ALE,DeEF,ElemF,ElemS,Phy,Sol,meshString);
            end
            
        elseif Sol.Order ==1
            [DeEF,Sol]=feue2(ALE,DeEF,ElemF,ElemS,Phy,Sol,meshString); %%% most important step and function might be..
        else
            disp('Order greater than 2 not yet implemented');
            return;
        end
         if (strcmp(meshString{1},'FlexFlatPlateF')==1)...
                ||(strcmp(meshString{1},'FlexFlatPlateTh02F')==1)...
                ||(strcmp(meshString{1},'FlexFlatPlateTh005F')==1)
                % vorticity
            b1=fU_gradp(DeEF,ElemF,Sol.U(:,1,1));
            b2=fU_gradp(DeEF,ElemF,Sol.U(:,2,1));
            vort=GM.fM\(b2(:,1)-b1(:,2));
         end
        
        if strcmp(meshString{1}(1:4),'FPCf')==1 || strcmp(meshString{1},'BenchMarkF')==1 || (strcmp(meshString{1},'FlexFlatPlateF')==1)

    %         Sol.ad(nn,:)=Sol.eta(175,:,1); %%% why only 175... or from where did we get the 175..

            [Sol.ldc(nn,:)]=ldc(ElemF,Phy,1,Sol.U(:,:,1),Sol.p,DeEF,meshString);
            nn=nn+1;
        end
        
        if (strcmp(meshString{1},'BenchMarkF')==1 | strcmp(meshString{1}(1:4),'FPCf')==1) & mod(nit,50)==0
                save('SOL.mat','Sol');
        end
        
        if strcmp(meshString{1}(1:4),'FPCf')==1
            DataEta(:,:,Counter) = Sol.eta(:,:,1);
            Counter = Counter+1;
        end
            
%         if abs(Sol.t-1)<1e-6
%             eta = Sol.eta;
%             save(cell2mat(titleat1(loop)),'eta');
%             save('ElemS.mat','ElemS');
%         end
        
        if (strcmp(meshString{1},'FlexFlatPlateF')==1)...
                ||(strcmp(meshString{1},'FlexFlatPlateTh02F')==1)...
                ||(strcmp(meshString{1},'FlexFlatPlateTh005F')==1)
%             if mod(nit,25)==0
%                 ALE = MeshS(ALE.Elem,DeEF.dof);
%     %             ALE = SVKFluidMeshUpdateLinear(ALE.Elem,DeEF.dof);
%                 Sol.dof = DeEF.dof;
%                 Sol.dis=zeros(ElemS.nDof,2); % displacement
%             end
            if mod(nit,20)==0
                DataU(:,:,Counter) = Sol.U(:,:,1);
                DataEta(:,:,Counter) = Sol.eta(:,:,1);
                DataEtat(:,:,Counter) = Sol.etat(:,:,1);
        %         DataP(:,:,Counter) = Sol.p;
                DataME(:,:,Counter) = Sol.me(:,:,1);
        %         DataVort(:,:,Counter) = vort;
                Counter = Counter+1;
            end

            if mod(nit,100)==0
                save(['DataU',int2str(DataCount),'.mat'],'DataU');
                save(['DataEta',int2str(DataCount),'.mat'],'DataEta');
                save(['DataEtat',int2str(DataCount),'.mat'],'DataEtat');
%                 save(['DataP',int2str(DataCount),'.mat'],'DataP');
                save(['DataME',int2str(DataCount),'.mat'],'DataME');
%                 save(['DataVort',int2str(DataCount),'.mat'],'DataVort');
                save('SolverInputs.mat','ALE','DeEF','ElemF','ElemS','Phy','meshString','nit','GM')
                MU = Sol.MU;
                save('MU.mat','MU');
    %             
    %             ALE.eta=zeros(ElemS.nDof,2); % displacement
    %             ALE.dof = DeEF.dof;
    %             ALE = SVKFluidMeshUpdateLinear(ALE);
                DataCount=DataCount+1;
                Counter = 1;
                DataU = [];
                DataEta = [];
                DataEtat = [];
                DataP = [];
                DataDeF = [];
                DataVort = [];
                t
            end
        end           
   
    
%     if mod(nit,round(plotDt/Sol.delt))==0
%         figure(3); clf;
%         subplot(3,1,1);
%         quiver(DeEF.dof(:,1),DeEF.dof(:,2),Sol.U(:,1,1),Sol.U(:,2,1));
%         hold on; 
%         showmesh(ElemS.dof(1:ElemS.nVertex,:)+Sol.eta(1:ElemS.nVertex,:,1),...
%             ElemS.elem2dof(:,1:3));
%         hold on; 
%         text(1.3,0.2,num2str(Sol.t));
%         axis equal;  % drawnow;  
% 
%         % stream function
%         set0index=1;
%         U_qt=uGraduAtQuadPts(Sol.U(:,:,1),DeEF,ElemF);
%         % fq_divu calculate <U,grad phi>. We now let u->v, v->-u;
%         fq=fq_divu(U_qt(:,:,2),-U_qt(:,:,1),DeEF,ElemF);
%         fq(set0index)=0;
%         GM.fS(set0index,:)=0; GM.fS(:,set0index)=0; GM.fS(set0index,set0index)=1;
%         strea=GM.fS\fq;
%         subplot(3,1,2);
%         FPCplot(DeEF,ElemF,ElemS,Sol,strea,plotFrame,50,[]);
%         
%         % vorticity
%         b1=fU_gradp(DeEF,ElemF,Sol.U(:,1,1));
%         b2=fU_gradp(DeEF,ElemF,Sol.U(:,2,1));
%         vort=GM.fM\(b2(:,1)-b1(:,2));
%         subplot(3,1,3);
%         % so far, the FPCplot is for vorticity only
%         FPCplot(DeEF,ElemF,ElemS,Sol,vort,plotFrame,100,[]);
%         drawnow;
%         if Sol.movie==1 && Sol.t<=8
% 
%             % vorticity
%             b1=fU_gradp(DeEF,ElemF,Sol.U(:,1,1));
%             b2=fU_gradp(DeEF,ElemF,Sol.U(:,2,1));
%             vort=GM.fM\(b2(:,1)-b1(:,2));
%             figure(1); clf;
%             % so far, the FPCplot is for vorticity only
%             FPCplot(DeEF,ElemF,ElemS,Sol,vort,plotFrame,200,'Vorticity');
%             drawnow;
% 
% %             % stream function
% %             set0index=1;
% %             U_qt=uGraduAtQuadPts(Sol.U(:,:,1),DeEF,ElemF);
% %             % fq_divu calculate <U,grad phi>. We now let u->v, v->-u;
% %             fq=fq_divu(U_qt(:,:,2),-U_qt(:,:,1),DeEF,ElemF);
% %             fq(set0index)=0;
% %             fS(set0index,:)=0; fS(:,set0index)=0; fS(set0index,set0index)=1;
% %             strea=fS\fq;
% %             figure(2); clf;
% %             FPCplot(DeEF,ElemF,ElemS,Sol,strea,plotFrame,'Stream Function');
% 
%             if Sol.t<=8
%                 figure(1);
%                 
%                 aviobj = addframe(aviobj,gcf);
%                 drawnow; 
% %                 MovieData(moviecount)=getframe(gcf);
% %                 moviecount=moviecount+1;
%             end
%         end
%     end % end of mod(nit,round(plotDt/Sol.delt))==0  
%     
% %     this step has been commented by me on 23/12/2012
% %     tt=[1,2,3,4,5,6,7,8,9,10,15,20];
% %     if min(abs(Sol.t-tt))<Sol.delt/3
% %         str=['flelueGCLT',num2str(Sol.t,3),meshString{1},...
% %             'NGR',num2str(Sol.nGlobalRefine),...
% %             'Dt',num2str(Sol.delt),'Deg',num2str(ElemF.deg)];
% %         save([str,'.mat']);
% %     end
end
toc
disp(['hmin=',num2str(ElemF.hm(1))]);
disp(['hmax=',num2str(ElemF.hm(2))]);
hold off;

% viobj = close(aviobj);
% 
% if Sol.movie==1
%    movie2avi(MovieData,'movie1.avi','compression', 'none','fps',1); % generate the .avi movie
% end

% res = Sol.ldc;
% res1 = Sol.ldc1;
% save('res.mat','res','res1');
if strcmp(meshString{1},'FS4f')==1
    eta = Sol.eta;
    U = Sol.U(:,:,1);
    P = Sol.p(:,:,1);
    V = Sol.etat(:,:,1);
    save(cell2mat(title(loop)),'eta');
    save(cell2mat(titleU(loop)),'U');
    save(cell2mat(titleP(loop)),'P');   
    save(cell2mat(titleV(loop)),'V'); 
end

if strcmp(meshString{1},'BenchMarkF')==1
    DelP = Sol.DelP;
    save('DelP.mat','DelP');
    save('SOL.mat','Sol');
end

if strcmp(meshString{1}(1:4),'FPCf')==1
    save('Eta.mat','DataEta');
    eta = Sol.eta;
    U = Sol.U(:,:,1);
    P = Sol.p(:,:,1);
    V = Sol.etat(:,:,1);
    save(cell2mat(title(loop)),'eta');
    save(cell2mat(titleU(loop)),'U');
    save(cell2mat(titleP(loop)),'P');   
    save(cell2mat(titleV(loop)),'V'); 
end

end

exit;
% figure(10)
% plo