%% Implementation of RBF-FD method for solving Tumor growth model (Prostate Tumor)
clearvars; close all; clc; format shorte;
%% Free energy function & its derivative
M=16;
F=@(u) M*u.^2.*(1-u).^2;
syms u;
dF=inline(diff(F(u),u));
%% Spatial step
a=0; b=2; NN =128; h=(b-a)/NN;
%% Collocation points
[X,Xb,XI,N,Nb,NI]=points('Regular',h,a,b); % see the function "points"
%% Dimensionless parameters
tau1=0.01;
a1=0.16;  a2=100;             a3=600;          a4=600; 
b1=5;     b2=365*2.70;        b3=365*2.75;     b4=1000;
c1=0.16;  c2=6.25*10^(-3);    c3=15*c2;        c4=100;
%%
T=0.6; dt=1e-4; MT=T/dt;% Time step & final time (The user can be changed the final time)
tol=1e-6; maxit=100;     % Tolerance & the maxiumun number of iterations for BiCGSTAB method
ns=13; c=15;             % the number of points in each stencil & shape parameter
IDX = knnsearch([X(:,1),X(:,2)],[X(:,1),X(:,2)],'K',ns); % Use Knnsearch
[A,Axx]= weights_RBF_FD([X(:,1),X(:,2)],[X(:,1),X(:,2)],IDX,c); 
AII=A(Nb+1:N,:);
%% The required differentiation matrices (if it is necessary)
ak1=A;            A1=ak1(Nb+1:N,Nb+1:N);
ak2=Axx;          A2=ak2(Nb+1:N,Nb+1:N);
AT=[(1+dt*a4)*A1-dt*a1*(A2)  -dt*a3*A1  sparse(NI,NI)
    dt*b3*A1  (1+dt*b4)*A1-dt*b1*(A2)  sparse(NI,NI)
    (dt*c2-dt*c3)*A1  sparse(NI,NI)  (1+dt*c4)*A1-dt*c1*(A2)];
%% The initial conditions
u=zeros(N,1);
ind=find(((X(:,1)-1).^2/(0.1)^2)+((X(:,2)-1).^2/(0.15)^2)<=1);
u(ind)=1;
sigma=zeros(N,1);
p=zeros(N,1);
%% The solution at different times
sz=[120 220];   % surface grid parameters
[ll,tt]=meshgrid(linspace(a,b,sz(2)),linspace(a,b,sz(1)));
xx=[ll(:) tt(:)];
%% **************** Interpolate to the grid ************* %%
IDXI = knnsearch([X(:,1),X(:,2)],[xx(:,1),xx(:,2)],'K',ns);
[AI,~]= weights_RBF_FD([X(:,1),X(:,2)],[xx(:,1),xx(:,2)],IDXI,c);
yy=reshape(xx(:,2),sz);
xx=reshape(xx(:,1),sz);
%% Preconditioner for different times
setup.type='nofill';
[Lu,Uu] = ilu(AT,setup);
t=0; m=0;
kp=0; kS=3; % the counters for drawing subplot
while t<=T
    rhs=[(AII*u-dt*a2*dF(AII*u));(AII*sigma+dt*b2);(AII*p+dt*c2)];
    U=bicgstabl(AT,rhs,tol,maxit,Lu,Uu);
    uI=U(1:NI);
    sigmaI=U(NI+1:2*NI);
    pI=U(2*NI+1:3*NI);
    
    u=[zeros(Nb,1);uI];
    sigma=[zeros(Nb,1);sigmaI];
    p=[zeros(Nb,1);pI];
    m=m+1;
    t=m*dt
    
    if  (t==0.2 || t==0.4 || t==0.6 )
        uapp=AI*u;
        induapp1=find(uapp<0);
        uapp(induapp1)=0;
        uu=reshape(uapp,sz);
        
        sigmaapp=AI*sigma;
        induapp2=find(sigmaapp<0);
        sigmaapp(induapp2)=0;
        sigmasigma=reshape(sigmaapp,sz);
        
        papp=AI*p;
        induapp3=find(papp<0);
        papp(induapp3)=0;
        pp=reshape(papp,sz);
        
        kp=kp+1;
        subplot(kS,3,kp),pcolor(xx,yy,uu)
        shading interp;
        str = ['t=',num2str(t)];
        title(str);
        colormap Copper;
        
        kp=kp+1;
        subplot(kS,3,kp),pcolor(xx,yy,sigmasigma)
        shading interp;
        str = ['t=',num2str(t)];
        title(str);
        colormap Copper;
        
        kp=kp+1;
        subplot(kS,3,kp),pcolor(xx,yy,pp)
        shading interp;
        str = ['t=',num2str(t)];
        title(str);
        colormap Copper;
        hold all;
        
    end
     if isnan(u)
        break;
    end
end
% End of program