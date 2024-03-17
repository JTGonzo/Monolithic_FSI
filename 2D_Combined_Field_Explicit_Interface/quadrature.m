function Quad=quadrature(deg)
% first, generate eQuad
switch deg
    case 1
        [x,y,kappa]=Dunavant(3);
        x=x'; y=y'; kappa=kappa/2;
        Nquad=length(x);                
    case 2
        [x,y,kappa]=Dunavant(5);
        x=x'; y=y'; kappa=kappa/2;
        Nquad=length(x);
        a=sqrt(15);   
        
        % Not symmetric ? Should comment out the following lines ?
%         pos=[6-a,9+2*a,6+a,9-2*a,7]/21;
%         x=pos([1,2,1,3,3,4,5]);
%         y=pos([1,1,2,4,3,3,5]);
%         wts=[155-a,155+a,270]/2400;
%         kappa=wts([1,1,1,2,2,2,3])';
    case 3
        % the quad pts of dunavant(11) is not always inside the triang
        [x,y,kappa]=Dunavant(8);
        x=x'; y=y'; kappa=kappa/2;
        Nquad=length(x);        
    case 4
        % the quad pts of dunavant(11) is not always inside the triang
        [x,y,kappa]=Dunavant(12);
        x=x'; y=y'; kappa=kappa/2;
        Nquad=length(x);
    case 5
        % the quad points of dunavant(14) remain inside the triangle
        [x,y,kappa]=Dunavant(14);
        x=x'; y=y'; kappa=kappa/2;
        Nquad=length(x);        
    otherwise
        disp(['not implemented yet']);
        pause;
end

[eQuad.psi,eQuad.psi_x,eQuad.psi_y]=basefun(deg,x,y);
eQuad.N=Nquad;
eQuad.kappa=kappa;
eQuad.xy=[x;y];

bdQuad.N=deg+1;  %deg+2;
[x,delta] = legendre_compute(bdQuad.N);
bdQuad.delta=delta'/2;
bdQuad.x=(x+1)/2;
[bdQuad.psi,bdQuad.psi_x,bdQuad.psi_y]=basefun(deg,bdQuad.x,1-bdQuad.x); % y=1-x;

x=linspace(0,1,deg+1); y=1-x;
[bdGrid.psi,bdGrid.psi_x,bdGrid.psi_y]=basefun(deg,x,y);

Quad.eQuad=eQuad;
Quad.bdQuad=bdQuad;
Quad.bdGrid=bdGrid;