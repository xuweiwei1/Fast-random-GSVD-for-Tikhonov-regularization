function [Q] = random_qb(A,tol)
[m,n] = size(A);
eps=tol;
Q = zeros(m,n);
R = Q(1:n,:);
% count=10;
% k = round(n/count); %% 设置每次取n/20列进行QR计算
k=50;
count=round(n/k);
a = sqrt(3);
for j = 1:count+1
    ost = (j-1)*k+1; oen = min(j*k,n); oenk = ost-1;
    if oenk>size(Q,2)
        break
    end
    v = A*((-a) + 2*a.*rand(n,oen-ost+1));% (a*rands(n,oen-ost+1)); %% 使用的随机矩阵是均匀矩阵
    v = v - Q(:,1:oenk) * (Q(:,1:oenk)' * v);
    [Q(:,ost:oen),R(ost:oen,ost:oen)] = qr(v,0);
    i = find(abs(diag(R(ost:oen,ost:oen)))<eps,1);
    if ~isempty(i)
        Q = Q(:,1:oen-k+i-1);
        return
    end
end
end