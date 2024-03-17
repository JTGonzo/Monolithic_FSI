function Sol = elSVK(Elem,Phy,Sol,theta)
global GM
if Sol.t==0 % initialize eta
    Sol.eta=zeros(Elem.nDof,2,3); % eta(:,:,1:3) is the eta at t=n+1,n,n-1
    Sol.etaAtItf=Sol.eta(Elem.itfNode,:,1);
    if ~isempty(Sol.deg) % strong form at quad pts
        Sol.sTraAtItf=zeros(Elem.Quad.bdQuad.N*size(Elem.itfElem,1),4);
    else % weak form at interface dof
        Sol.sTraAtItf=zeros(2*size(Elem.itfNode,1),1);
    end
    return
end % end of t=0
%==========================================================================

% Note that sdof(sitfNode,:)-dof(itfNode,:)=0
temp=theta*(Sol.eta(Elem.itfNode,:,1))+...
    (1-theta)*(Sol.eta(Elem.itfNode,:,2)+Sol.delt*Sol.uAtItf);
% Sol.eta(:,:,1)=Sol.eta(:,:,2); % init guess. This should NOT go before temp=theta*......
Sol.eta(Elem.bdDirN,:,1)=0;
Sol.eta(Elem.itfNode,:,1)=temp; % itfNode is part of bdDirN

eta=ones(2*Elem.nFreeN,1);
freeN=[Elem.freeN;Elem.freeN+Elem.nDof];
while max(abs(eta))>1e-7
    [a,S,b]=initMatrixSVK(Elem,Phy,Sol.eta(:,:,1));
    % eta^n+1 - 2 eta^n + eta^n-1 - delt^2(gg - b) = 0
    temp=GM.sM*(Sol.eta(:,:,1)-2*Sol.eta(:,:,2)+Sol.eta(:,:,3)+...
        (Sol.delt^2/Phy.srho)*repmat([0,Phy.g],Elem.nDof,1));
    bb=temp(:)+(Sol.delt^2/Phy.srho)*b;
    E=[GM.sM,zeros(Elem.nDof,Elem.nDof);zeros(Elem.nDof,Elem.nDof),GM.sM]+Sol.delt^2*S;
    eta=E(freeN,freeN)\bb(freeN);
    Sol.eta(Elem.freeN,:,1)=Sol.eta(Elem.freeN,:,1)-reshape(eta,[],2);
end

Sol.etaAtItf=Sol.eta(Elem.itfNode,:,1);
% weak form at quad pts
Sol.sTraAtItf=(Phy.srho/Sol.delt^2)*bb([Elem.itfNode;Elem.itfNode+Elem.nDof]);