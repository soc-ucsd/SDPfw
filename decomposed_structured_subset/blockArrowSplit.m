function [model, model_split] = blockArrowSplit(m,nBlk,BlkSize,ArrowHead)
% Is a modification of blockArrowMultCones
% Setup problem for block-arrow sdp with multiple cones. Inputs:
% m        : number of equality constraints
% nBlk     : number of diagonal blocks in block-arrow matrix data (vector of length nCones)
% BlkSize  : size of each diagonal block (vector of length nCones)
% ArrowHead: size of head of arrow pattern (vector of length nCones)
% r        : rank of optimal X
% num_c    : number of cost functions desired

nCones = 1;
num_c = 1;
% cone
K.f = 0;
K.l = 0;
K.q = [];
K.s = zeros(1,nCones);

% Sparsity pattern of each cone
Spa = cell(nCones,1);

%modify this later
cliques = cell(nBlk, 1);
%idx
%n = nBlk*BlkSize + ArrowHead;
%idx_corner = reshape(n*((nBlk(k)*BlkSize(k)+1:n))-1)' + (nBlk(k)*BlkSize(k)+1:n), [], 1);

for k = 1:nCones
    n = nBlk(k)*BlkSize(k)+ArrowHead(k);
    Spa{k} = zeros(n);
    for i = 1:nBlk(k)
        Spa{k}((i-1)*BlkSize(k)+1:BlkSize(k)*i, ...
               (i-1)*BlkSize(k)+1:BlkSize(k)*i) = ones(BlkSize(k));
        cliques{i} = [(i-1)*BlkSize(k)+1:BlkSize(k)*i, nBlk(k)*BlkSize(k)+1:n];
    end
    Spa{k}(nBlk(k)*BlkSize(k)+1:n,:) = 1;
    Spa{k}(:,nBlk(k)*BlkSize(k)+1:n) = 1;
    
    
    
    % Set cone size
    K.s(k) = n;
end



% Data
%At = [];
At = sparse(n^2, m);
for i = 1:m
%     Ai = [];
    %for k = 1:nCones
        M = 10*sprandsym(Spa{k});  % random symmetric data with given sparsity pattern
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
Temp = 10*sprandsym(Spa{k});
Temp = Temp + (-min(eig(full(Temp)))+1)*speye(size(Spa{k}));
X_star = Temp(:);

b = At'*X_star;

% stictly feasible dual point

S = cell(nCones,1);
for k = 1:nCones
    Temp = 10*sprandsym(Spa{k});
    Temp = Temp + (-min(eig(full(Temp)))+1)*speye(size(Spa{k}));
    S{k} = Temp(:);
end

S_cat = vertcat(S{:});
c = sparse(size(S_cat, 1), num_c);

y = rand(m,1);

if (nargin >= 6) && num_c > 1
    for i = 1:num_c
        y = rand(m,1);
        c(:, i) = S_cat + At*y;
    end
else
    y = rand(m,1);
    c = S_cat + At*y;
end


model.At = At;
model.b = b;
model.c = c;
model.K = K;

%now do the split
model_split.At = [];
model_split.b  = model.b;
model_split.c  = [];
model_split.K  = model.K;
model_split.K.s = [];
n = nBlk*BlkSize+ArrowHead;

%clique consistency equality constraints
At_eq = [];
b_eq = [];
Count = 0;

for i = 1:nBlk
    cli = cliques{i};
    cli_corner = cli(end-ArrowHead+1:end);
    lc = BlkSize + ArrowHead;
    ind = reshape(n*(cli-1)' + cli, [], 1);
    
    ind_corner = lc*((lc-ArrowHead+1:lc)-1)' + (lc-ArrowHead+1:lc);
    
    Ati = At(ind, :);
    ci  = c(ind);
    
    %clique overlap constraints
    if i == 1        
        ind_corner_v0 = triu_vec(ind_corner);
    else
        %don't double-count cliques
        Ati(ind_corner, :) = 0;
        ci(ind_corner) = 0;
        
        %clique consistency equality constraints
        ind_corner_v = triu_vec(ind_corner) + Count;
        v_len = length(ind_corner_v);
                
        At_eq_curr = sparse([ind_corner_v0; ind_corner_v], [1:v_len, 1:v_len], [ones(v_len, 1), -ones(v_len, 1)], nBlk*lc^2, v_len);
        b_eq_curr = sparse(v_len, 1);
        
        At_eq = [At_eq At_eq_curr];
        b_eq  = [b_eq;  b_eq_curr];
    end
    
    model_split.At = [model_split.At; Ati];
    model_split.c  = [model_split.c; ci];
    
    model_split.K.s = [model_split.K.s lc];
    Count = Count + lc^2;
end

model_split.At = [model_split.At At_eq];
model_split.b  = [model_split.b; b_eq];

end

function v = triu_vec(K)
    %get upper triangular elements of matrix K
    %is inefficient, but I just want this to work.
    Kt = K.';
    m  = (1:size(K,1)).' >= (1:size(Kt,2));
    v  = Kt(m);
end
