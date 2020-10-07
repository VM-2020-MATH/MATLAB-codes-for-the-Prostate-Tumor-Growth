function [X,Xb,XI,N,Nb,NI]=points(char,h,a,b)
% Inputs:
% char can be chosen by the user Regular or distmesh (it needs to be downloaded this package that is free for the user)
% h is the maximum distance between two distinct points
% a & b are chosen here 0 & 2, respectively

% Outputs:
% X is the total number of collocation points with N
% Xb is the number of boundary points with Nb
% XI is the number of interior points with NI

switch (char)
    case ('Regular')
        X1=[a:h:b];
        Y1=[a b];
        [Xb1,Yb1]=meshgrid(X1,Y1);
        X2=[a b];
        Y2=[a+h:h:b-h];
        [Xb2,Yb2]=meshgrid(X2,Y2);
        Xb=[Xb1(:) Yb1(:);Xb2(:) Yb2(:)];
        [Xi,Yi]=meshgrid(a+h:h:b-h);
        XI=[Xi(:) Yi(:)];
        X=[Xb;XI];
        [Nb,~]=size(Xb);
        [NI,~]=size(XI);
        N=Nb+NI;
    case ('Distmesh')
        X1=[a:h:b];
        Y1=[a b];
        [Xb1,Yb1]=meshgrid(X1,Y1);
        X2=[a b];
        Y2=[a+h:h:b-h];
        [Xb2,Yb2]=meshgrid(X2,Y2);
        Xb=[Xb1(:) Yb1(:);Xb2(:) Yb2(:)];
        Nb=length(Xb);
        
        fd=@(p) drectangle(p,a+h,b-h,a+h,b-h);
        fh=@(p) (0.01+0.3*abs(dcircle(p,0,0,0)))+(0.025+0.3*abs(dpoly(p,[2.5,2.5; 2.5,2.5])));% For quasi-uniform points is the best choice.
        hp=h;
        [p,~]=distmesh2d(fd,fh,hp,[a+h,a+h;b-h,b-h],[a+h,a+h;b-h,a+h;a+h,b-h;b-h,b-h]);
        XI=p;                NI=length(XI);
        X=[Xb;XI];           N=length(X);
    otherwise
        error('There is no method')
end





