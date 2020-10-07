clearvars; clc; format shorte;
%% Spatial step
a=0; b=2; NN =128; h=(b-a)/NN;
%% Collocation points
[X,Xb,XI,N,Nb,NI]=points('Regular',h,a,b); % see the function "points"
ns=13; c=15;             % the number of points in each stencil & shape parameter
IDX = knnsearch([X(:,1),X(:,2)],[X(:,1),X(:,2)],'K',ns); % Use Knnsearch
[~,Axx]= weights_RBF_FD([X(:,1),X(:,2)],[X(:,1),X(:,2)],IDX,c); 
%% Drawing all the eigenvalues
a1=0.16; b1=5; % Special parameter for the first simulation
DDX=a1*Axx(Nb+1:N,Nb+1:N); 
Lambda1= eigs(DDX,NI,'LM');
Lambda2=(b1/a1)*Lambda1;
Lambda3=Lambda1;
Lambda=[Lambda1;Lambda2;Lambda3];
plot(real(Lambda),imag(Lambda),'+');
xlabel('$Re(\lambda)$','Interpreter','latex')
ylabel('$Im(\lambda)$','Interpreter','latex')
str = ['$N_{I}$=',num2str(NI)];
TitleH=title(str);
set(TitleH, 'Position')
hold on;