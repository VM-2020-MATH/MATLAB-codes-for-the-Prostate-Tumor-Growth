function [D,Dxx]= weights_RBF_FD(x,y,IDX,c)
% Hint: This is borrowed from the written code by Natasha
% Flyer and co-authors 

% Inputs: 
% x is the collocation points
% y is the point which the approximation should be computed
% IDX the indices for each stencil
% c is the shape parameter

% Outputs:
% the collocation matrix, D
% the differentiation matrix, Dxx

% Use MQ function and the second derivative in two-dimensional space
phi=@(r,c) (c.^2*r.^2+1).^(1/2);
d2phi=@(r,c) c.^2.*(c.^2.*r.^2+2)./(c.^2*r.^2+1).^(3/2);
[N,n] = size(IDX);  M=length(x);
w=zeros(N,n);  wxx=zeros(N,n);
for k= 1:N
    X = x(IDX(k,:) ,:); % Nodes in kth stencil
    R=sqrt((X(:,1)-X(:,1)').^2+(X(:,2)-X(:,2)').^2);
    A=phi(R,c);
    x1=y(k,:);         % Node, which the approximation is computed in it.
    r=sqrt((x1(:,1)-X(:,1)).^2+(x1(:,2)-X(:,2)).^2);
    L=phi(r,c);
    W = A\L;
    w(k,:)=W;
    
    d2xL=d2phi(r,c);
    Wxx = A\d2xL;
    wxx(k,:)=Wxx;
end
it = (1:N)'; it = it(:,ones(1,n));
D = sparse(it(:),IDX(:),w(:),N,M);
Dxx = sparse(it(:),IDX(:),wxx(:),N,M);