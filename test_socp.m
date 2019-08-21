
% test psd to socp

clear;clc;

Num = 1500;
A = [];
X = [];
Z = [];
for i = 1:Num
    temp = rand(2);
    At = temp + temp';
    Xt = eye(2);
    Zt = eye(2);    
    X = [X;Xt(:)];
    Z = [Z;Zt(:)];
    A = [A, At(:)'];
end
b = A*X(:);
y = rand(1);
c = y*A' + Z(:);
K.s = 2*ones(Num,1);

% original SDP with all 2 by 2 psd cones
[x,y,info] = sedumi(A',b,c,K);

% transform into second-order cones
[An,bn,cn,Kn,indsocp] = psd2socp(A,b,c,K);
[xn,yn,infon] = sedumi(An,bn,cn,Kn);

[info.cpusec,infon.cpusec]
[y - yn, cn'*xn - c'*x, norm(indsocp*xn - x)]

