function [model, model_sym] = sym_block_arrow_sdp(m, nBlk, BlkSize, ArrowHead, orbits, splits)
%number of matrices

 if nargin < 6
     splits = zeros(size(orbits));
 end

N = nBlk*BlkSize + ArrowHead;

orbit_ind = cellfun(@(x) (x-1)*BlkSize +  (1:BlkSize)', orbits , 'UniformOutput', false);
arrow_ind = nBlk*BlkSize + (1:ArrowHead);

[A, A_sym, U, block_perm] = sym_block_arrow(m+1, nBlk, BlkSize, ArrowHead, orbits, splits);



%% Unstructured Program
%Dual variable S >=0 and permutation invariant
S =  10*sprandsym(A{1});
S = (-min(eig(full(S)))+1)*speye(size(S));
S = reynolds_perm(S, orbit_ind, BlkSize, arrow_ind);
S_blk = U'*S*U;
S_sym = S_blk(block_perm, block_perm);
S_sym(abs(S_sym) < 1e-8) = 0; 
    

%Sf = S(:);

model = struct;
model.K.f = 0;
model.K.l = 0;
model.K.q = [];
model.K.s = N;

model.At = sparse(N^2, m);
for i = 1:m
    model.At(:,  i) = A{i}(:);    
end



%initial feasible X
X_star = 10*sprandsym(A{1});
X_star = X_star + (-min(eig(full(X_star)))+1)*speye(size(A{1}));
%X_star = Temp(:);
y = rand(m,1);
model.b = model.At'*X_star(:);


model.pars.fid = 1;


model.c = S(:) + model.At*y;

%% Symmetric Program
model_sym = struct;
orbit_weight = cellfun(@(x) length(x)-1, orbits);
orbit_split = orbits(logical(splits));
orbit_split_weight = cellfun(@(x) length(x), orbit_split);

% %size remaining in the arrowhead block
%orbit_rem = nBlk - sum(orbit_weight);
%sym_blocks = nnz(orbit_weight)
% 
% model_sym.At 
%model_sym.At;
%At_block = [];
%At_head = [];

%Next step: 
%I need to reynolds/diagonalize 

%USE S_SYM!



sub_block =  [];
w_block = [];
count = 0;
K = [];
Kw = [];
for i = 1:length(orbit_weight)
    if orbit_weight(i) >  0
        curr_size = BlkSize/(1+splits(i));         
        curr_ind = count + (1:curr_size);
        [meshx, meshy] = meshgrid(curr_ind);
        sub_block = [sub_block; meshx(:) meshy(:)];
        w_block = [w_block; ones(length(meshx)^2, 1)*orbit_weight(i)];
        count = count + orbit_weight(i)*curr_size;
        K = [K curr_size];
        Kw = [Kw curr_size*orbit_weight(i)];
    end
end

for i = 1:length(orbit_split_weight)
    %if orbit_split_weight(i) >  0
        curr_size = BlkSize/2;         
        curr_ind = count + (1:curr_size);
        [meshx, meshy] = meshgrid(curr_ind);
        sub_block = [sub_block; meshx(:) meshy(:)];
        w_block = [w_block; ones(length(meshx)^2, 1)*orbit_split_weight(i)];
        count = count + curr_size*orbit_split_weight(i);
        K = [K curr_size];
        Kw = [Kw curr_size*orbit_split_weight(i)];
    %end
end
    
%now the arrowhead
head_start = sum(Kw)+1;
head_ind = head_start:N;
[meshx, meshy] = meshgrid(head_ind);
sub_block = [sub_block; meshx(:) meshy(:)];
w_block = [w_block; ones(length(meshx)^2, 1)];
    

ind_block = sub2ind([N, N], sub_block (:, 1), sub_block(:, 2));


model.At_sym = [];
for k=1:m
    model_sym.At(:, k) =  A_sym{k}(ind_block).*w_block;
end
S_sym_block = S_sym(ind_block).*w_block;
model_sym.b = model.b;

model_sym.c = model_sym.At*y + S_sym_block;

model_sym.K =  model.K;
model_sym.K.s =  [K, N - head_start + 1];

%model_sym.At(:, i) = A_sym{i}(:);


model_sym.pars = model.pars;
%now generate the SDP
%y = 5*randn(m, 1);
%S = s

end
