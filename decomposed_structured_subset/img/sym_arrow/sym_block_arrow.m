function [MG, MG_diag, U, block_perm] = sym_block_arrow(m, nBlk, BlkSize, ArrowHead, orbits, splits)
%forms a set of matrices that are symmetrized block-arrow
% Inputs:
%   m:                          Number of matrices
%   nBlk, BlkSize, ArrowHead:   Block Arrow parameters
%   orbits:                     Orbits of the blocks (symmetry structure)
%   splits:                     Which blocks should have additional
%                               symmetry (Z2)

%total size
N = nBlk*BlkSize + ArrowHead;

if nargin < 5
    splits = zeros(size(orbits));
end

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
        curr_split = splits(i);
        num_orbits = size(curr_orbit, 2);

        Diag_block = zeros(BlkSize);
        Arrow_block = zeros(BlkSize, ArrowHead);

        %average over orbits
        %w =  num_orbits;
        %w = 1;
        w = sqrt(num_orbits);
        for j = 1:num_orbits
            curr_ind = curr_orbit(:, j);
            N_ind = length(curr_ind);
            if curr_split
                split_switch = [(N_ind/2+1):N_ind, 1:(N_ind/2)];
                split_ind = curr_ind(split_switch);
                Diag_block = Diag_block + (M{k}(curr_ind, curr_ind)/w + M{k}(split_ind, split_ind)/w)/2;
                Arrow_block = Arrow_block + (M{k}(curr_ind, arrow_ind)/w + M{k}(split_ind, arrow_ind)/w)/2; 
            else
                Diag_block = Diag_block + M{k}(curr_ind, curr_ind)/w;
                Arrow_block = Arrow_block + M{k}(curr_ind, arrow_ind)/w;
            end
            
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
ind_split = [];
ind_last = [];
for i = 1:length(orbits)
    curr_orbit = orbit_ind{i};    
    num_orbits = size (curr_orbit, 2);
    if splits(i)
        ind_first = [ind_first; curr_orbit(1:BlkSize/2, 1)]; 
        
        ind_split = [ind_split; reshape(curr_orbit((BlkSize/2+1):end, :), [], 1)];
        if num_orbits > 1        
            ind_last  = [ind_last; reshape(curr_orbit(1:(BlkSize/2), 2:end ), [], 1)]; 
        end
    else 
        ind_first = [ind_first; curr_orbit(:, 1)];          

        if num_orbits > 1        
            ind_last  = [ind_last; reshape(curr_orbit(:, 2:end ), [], 1)]; 
        end

    end

    
    
    
    
    U_dct(orbits{i}, orbits{i}) = dctmtx(num_orbits);
end
if any(splits) 
    U_dct_part = full(kron(U_dct', speye(2)));
    F = dctmtx(2);
    for i = 1:length(orbits)
        if splits(i)
            for j = 1:length(orbits{i})
                orb_curr = orbits{i}(j);
                curr_ind = 2*orb_curr +  [-1, 0] ;
                U_dct_part(curr_ind, :) = F*U_dct_part(curr_ind, :);
            end
        end        
    end
    U_dct_kron =  kron(U_dct_part, speye(BlkSize/2));
else
    U_dct_kron = kron(U_dct', speye(BlkSize));
end
U = [U_dct_kron, sparse(nBlk*BlkSize, ArrowHead);
    sparse(ArrowHead, nBlk*BlkSize), speye(ArrowHead)];

%MG_blk = U'*MG*U;
%MG_blk = cellfun(@(x) U'*x*U, MG);

%now shuffle indices around to form a block diagonal matrix
block_perm = [ind_last; ind_split; ind_first; arrow_ind'];
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
