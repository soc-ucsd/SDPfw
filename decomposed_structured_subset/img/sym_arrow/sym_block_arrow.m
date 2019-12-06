function [MG, MG_diag, U, block_perm] = sym_block_arrow(m, nBlk, BlkSize, ArrowHead, orbits)


%total size
N = nBlk*BlkSize + ArrowHead;

visualize = 1;

orbit_ind = cellfun(@(x) (x-1)*BlkSize +  (1:BlkSize)', orbits , 'UniformOutput', false);
arrow_ind = nBlk*BlkSize + (1:ArrowHead);

nCones = 1;
num_c = 1;
% % cone
% K.f = 0;
% K.l = 0;
% K.q = [];
% K.s = zeros(1,nCones);

% Sparsity pattern of each cone


%modify this later
cliques = cell(nBlk, 1);

%sparsity pattern
for k = 1:nCones
    n = nBlk(k)*BlkSize(k)+ArrowHead(k);
    Spa = zeros(n);
    for i = 1:nBlk(k)
        Spa((i-1)*BlkSize(k)+1:BlkSize(k)*i, ...
               (i-1)*BlkSize(k)+1:BlkSize(k)*i) = ones(BlkSize(k));
        cliques{i} = [(i-1)*BlkSize(k)+1:BlkSize(k)*i, nBlk(k)*BlkSize(k)+1:n];
    end
    Spa(nBlk(k)*BlkSize(k)+1:n,:) = 1;
    Spa(:,nBlk(k)*BlkSize(k)+1:n) = 1;           
    % Set cone size
    %K.s(k) = n;
end

M = cell(m, 1);

%M = 10*sprandsym(Spa);  % random symmetric data with given sparsity pattern

%now perform the reynolds operator and figure out block diagonalizer
MG = M;


%first and subsequent indices for each block, useful in block
%diagonalization
for k = 1:m
    M{k} = 10*sprandsym(Spa);
    MG{k} = M{k};
    for i = 1:length(orbit_ind)
        curr_orbit = orbit_ind{i};
        num_orbits = size(curr_orbit, 2);

        Diag_block = zeros(BlkSize);
        Arrow_block = zeros(BlkSize, ArrowHead);

        %average over orbits
        %w =  num_orbits;
        %w = 1;
        w = sqrt(num_orbits);
        for j = 1:num_orbits
            curr_ind = curr_orbit(:, j);
            Diag_block = Diag_block + M{k}(curr_ind, curr_ind)/w;
            Arrow_block = Arrow_block + M{k}(curr_ind, arrow_ind)/w;
        end

        %assign back into matrix
        for j = 1:num_orbits
            curr_ind = curr_orbit(:, j);
            MG{k}(curr_ind,  curr_ind) = Diag_block;
            MG{k}(curr_ind, arrow_ind) = Arrow_block;
            MG{k}(arrow_ind, curr_ind) = Arrow_block';
        end

        %Block Diagonalization
        %raw_orbit = sparse(orbits{i}, orbits{i}, ones(size(orbits{i})), nBlk, nBlk);
        %U_dct(orbits{i}, orbits{i}) = dctmtx(num_orbits);    
    end
end

%DCT to compactify blocks
U_dct = sparse(nBlk);
%block_weight = ones(length(orbits));
ind_first = [];
ind_second = [];
ind_last = [];
for i = 1:length(orbits)
    curr_orbit = orbit_ind{i};
    num_orbits = size (curr_orbit, 2);
    ind_first = [ind_first curr_orbit(:, 1)];    
    if num_orbits > 1
        ind_second = [ind_second curr_orbit(:, 2)];
        ind_last  = [ind_last  curr_orbit(:, 2:end )]; 
    end
    U_dct(orbits{i}, orbits{i}) = dctmtx(num_orbits);    
end
U_dct_kron = kron(U_dct', speye(BlkSize));
U = [U_dct_kron, sparse(nBlk*BlkSize, ArrowHead);
    sparse(ArrowHead, nBlk*BlkSize), speye(ArrowHead)];

%MG_blk = U'*MG*U;
%MG_blk = cellfun(@(x) U'*x*U, MG);

%now shuffle indices around to form a block diagonal matrix
block_perm = [ind_last(:); ind_first(:); arrow_ind'];
%MG_diag = MG_blk(block_perm, block_perm);
%MG_diag(abs(MG_diag) < 1e-8) = 0;
%MG_blk = cellfun(@(x) U'*x*U, MG);

MG_diag = cell(size(MG));
for k = 1:m
    MG_blk = U'*MG{k}*U;
    MG_perm = MG_blk(block_perm, block_perm);
    MG_perm(abs(MG_perm) < 1e-8) = 0; 
    
    MG_diag{k} = MG_perm;
end

end
