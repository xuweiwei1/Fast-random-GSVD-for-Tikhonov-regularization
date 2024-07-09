%%首先建立随机GSVD，针对A是低秩矩阵
function [x_k,Q,miu]=TIK_drgsvd1(A,L,tol,b)
[~,n]=size(A);
[~,n1]=size(L);
if (n1 ~= n)
    error('No. columns in A and L must be the same')
end
Q1=random_qb(full(A'),tol);
B=A*Q1;
Q=random_qb(B,tol);
C1=Q'*A*Q1;
[U,~,X,C,D]=fgsvd(C1,full(L)*Q1);
UU=Q*U;
[p,n2]=size(C1);
w=n2-p+1;
X1=Q1/(X);
    c=diag(C(:,w:end));
    s=diag(D(w:end,w:end));
    s=s(1:p);
    sm=[c,s];
    miu = gcv(UU,sm,b);
    W=c./(c.^2+miu^2*s.^2);
    beta=UU'*b;
    x_k=X1(:,w:end)*(W.*beta);
end