rng(10, 'twister')

%number of matrices
m =  80;
ArrowHead = 10;
BlkSize = 10;

% nBlk = 6;
% orbits = {[1,2,3,4],5, 6};
% 
% visualize = 1;
% 

nBlk = 8;
orbits = {[1,2,3,4], 5, 6, [7,8]};
% nBlk = 10;
% orbits = {[1,2,3], 4, [5, 6], [7,  8,9, 10]};
%orbits = {[1,3,4,9], [2,10], [5,8],  6, 7};
%orbits = {1:10};

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


if visualize
    %plot orbit structure
for k = 1:1
    figure(2)
    clf
    subplot(1, 3, 1)
    hold on
    Colors = linspecer(length(orbits));
    for i = 1:length(orbits)
        curr_orbit =  orbit_ind{i};
        curr_color = Colors(i, :);
        for j = 1:size(curr_orbit, 2)
            rectangle('Position', [curr_orbit(1, j)-1, curr_orbit(1, j)-1, BlkSize, BlkSize],...
                'FaceColor', curr_color)
            rectangle('Position', [curr_orbit(1, j)-1, nBlk*BlkSize, BlkSize, ArrowHead],...
                'FaceColor', curr_color)
            rectangle('Position', [nBlk*BlkSize,  curr_orbit(1, j)-1, ArrowHead, BlkSize],...
                'FaceColor', curr_color)            
        end
    end
    xlim([0,N]);
    ylim([0,N]);
    axis square
    
    rectangle('Position', [nBlk*BlkSize, nBlk*BlkSize, ArrowHead, ArrowHead], 'FaceColor', 'k')
    set(gca, 'YDir','reverse')
    title('Permutation Orbits', 'FontSize', 22)
%    plot elements of *-algebra
    %figure(1)    
    
        %figure('units','normalized','outerposition',[0 0 1 1])
        %clf
%         subplot(1,3,1)
%         M_display = M{k};
%         M_display(M_display == 0) = NaN;
%         imagesc(M_display)
%         axis square
%         title('Sparse Random')
        subplot(1,3,2)
        %imagesc(MG)
        M_display = MG{k};
        M_display(M_display == 0) = NaN;
        imagesc(M_display)
        axis square
        title('Group-Invariant Matrix', 'FontSize', 22)
        subplot(1,3,3)
        %imagesc(MG_diag)
        M_display = MG_diag{k};
        M_display(M_display == 0) = NaN;
        imagesc(M_display)
        axis square
        title('Block Diagonalization', 'FontSize', 22)
        %pause(0.5);
        fname = strcat('sym_block_arrow-', num2str(k));
        %export_fig(fname, '-png', '-m', '4');
       %keyboard
        
    end
end