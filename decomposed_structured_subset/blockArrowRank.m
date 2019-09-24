function [At,b,c,K, X_star] = blockArrowRank(m,nBlk,BlkSize,ArrowHead, r)
% Is a modification of blockArrowMultCones
% Setup problem for block-arrow sdp with multiple cones. Inputs:
% m        : number of equality constraints
% nBlk     : number of diagonal blocks in block-arrow matrix data (vector of length nCones)
% BlkSize  : size of each diagonal block (vector of length nCones)
% ArrowHead: size of head of arrow pattern (vector of length nCones)
% r        : rank of optimal X

nCones = 1;

% cone
K.f = 0;
K.l = 0;
K.q = 0;
K.s = zeros(1,nCones);

% Sparsity pattern of each cone
Spa = cell(nCones,1);
for k = 1:nCones
    n = nBlk(k)*BlkSize(k)+ArrowHead(k);
    Spa{k} = zeros(n);
    for i = 1:nBlk(k)
        Spa{k}((i-1)*BlkSize(k)+1:BlkSize(k)*i, ...
               (i-1)*BlkSize(k)+1:BlkSize(k)*i) = ones(BlkSize(k));
    end
    Spa{k}(nBlk(k)*BlkSize(k)+1:n,:) = 1;
    Spa{k}(:,nBlk(k)*BlkSize(k)+1:n) = 1;
    
    % Set cone size
    K.s(k) = n;
end



% Data
%At = [];
At = sparse(n^2, m);
%At = zeros(
for i = 1:m
%     Ai = [];
    %for k = 1:nCones
        M = 100*sprandsym(Spa{k});  % random symmetric data with given sparsity pattern
%         Ai = [Ai; M(:)];        % concatenate
%     end
    At(:, i) = M(:);
end

%I suggest that you can generate random sparse SDPs that admit low-rank 
%solutions, and then try the reweighted algorithm again. Basically, you can
%generate random sparse Ai, choose a low-rank X >= 0, and determine 
%bi = <Ai,X> that guarantees primal feasible. For dual feasibility, you can
%generate a sparse S > 0 and vector y, and determine C = S + \sum yi Ai. 
%In this case, we know there exists a low-rank solution. Then the algorithm
%has the potential to reduce the rank of obtained solution. 

%Optimal X:
R = randn(n, r);
X_star_mat = R*R';
X_star = X_star_mat(:);

% stictly feasible primal point
% X = cell(nCones,1);
% for k = 1:nCones
%     Temp = 10*sprandsym(Spa{k});
%     Temp = Temp + (-min(eig(full(Temp)))+1)*speye(size(Spa{k}));
%     X{k} = Temp(:);
% end
b = At'*X_star;

% stictly feasible dual point
y = rand(m,1);
S = cell(nCones,1);
for k = 1:nCones
    Temp = 10*sprandsym(Spa{k});
    Temp = Temp + (-min(eig(full(Temp)))+1)*speye(size(Spa{k}));
    S{k} = Temp(:);
end
c = vertcat(S{:}) + At*y;

end
