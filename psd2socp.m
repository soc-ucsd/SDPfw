function [An,bn,cn,Kn,indsocp] = psd2socp(A,b,c,K)
% Transform all 2 by 2 psd cone into second-order cone constraints
% Assume K.s only contains 2 by 2 psd cones or 2 by 2 psd cones come first in K.s

%% Input check

 if length(find(K.s == 2)) ~= length(K.s)
     error('psd2socp can only deal with all 2 by 2 psd cones ...')
 end

if size(A,1) ~= length(b) 
    A = A';
end
if ~isfield(K,'f') || isempty(K.f) 
    K.f = 0;
end
if ~isfield(K,'l') || isempty(K.l) 
    K.l = 0;
end
if ~isfield(K,'q') || isempty(K.q) 
    K.q = 0;
end

% Non PSD part
Anonpsd = A(:,1:K.f+K.l+K.q);
cnonpsd = c(1:K.f+K.l+K.q);

% psd data
Apsd = A(:,K.f+K.l+K.q+1:end);   % PSD data 
cpsd = c(K.f+K.l+K.q+1:end);     
   

%% 2-2 psd to socp
P = [0.5 0  0.5
     0  0.5  0;
     0  0.5  0;
     0.5 0 -0.5];
Pt = P';
P  = sparse(P);
Pt = sparse(Pt);
I  = speye(length(K.s));
P  = kron(I,P);
Pt = kron(I,Pt);
 
An = [Anonpsd,Apsd*P];
cn = [cnonpsd;Pt*cpsd];
bn = b;

Kn = K;
Kn.s = [];
if K.q(1) ~= 0
    Kn.q = [Kn.q,3*ones(1,length(K.s))];
else
    Kn.q = 3*ones(1,length(K.s));
end

indsocp = P;

end

