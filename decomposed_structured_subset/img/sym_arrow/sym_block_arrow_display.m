rng(10, 'twister')
visualize = 1;
%number of matrices
m =  1;
ArrowHead = 10;
BlkSize = 10;

% nBlk = 3;
% orbits = {[1, 2], 3};
%  splits = [1 1];
 
visualize = 1;

% nBlk = 6;
% orbits = {[1,2,3,4],5, 6};
% splits =  [1, 0, 0];
% visualize = 1;
% 

nBlk = 8;
orbits = {[1,2,3,4], 5, 6, [7,8]};
%splits = [1, 0, 1, 0 ];
splits = [ 0 0 0 0];


N = nBlk*BlkSize + ArrowHead;
orbit_ind = cellfun(@(x) (x-1)*BlkSize +  (1:BlkSize)', orbits , 'UniformOutput', false);
arrow_ind = nBlk*BlkSize + (1:ArrowHead);

% nBlk = 10;
% orbits = {[1,2,3], 4, [5, 6], [7,  8,9, 10]};
%orbits = {[1,3,4,9], [2,10], [5,8],  6, 7};
%orbits = {1:10};

[MG, MG_diag, U, block_perm] = sym_block_arrow(m, nBlk, BlkSize, ArrowHead, orbits, splits);


if visualize
    %plot orbit structure
for k = 1:m
    f2 = figure(2);
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    clf
    subplot(2, 2, 1)
    hold on
    Colors = linspecer(length(orbits) + nnz(splits));
    color_count = 1;
    for i = 1:length(orbits)
        curr_orbit =  orbit_ind{i};
        curr_color = Colors(color_count, :);
        curr_split = splits(i);
        for j = 1:size(curr_orbit, 2)            
            rectangle('Position', [curr_orbit(1, j)-1, curr_orbit(1, j)-1, BlkSize, BlkSize],...
                'FaceColor', curr_color)
            rectangle('Position', [curr_orbit(1, j)-1, nBlk*BlkSize, BlkSize, ArrowHead],...
                'FaceColor', curr_color)
            rectangle('Position', [nBlk*BlkSize,  curr_orbit(1, j)-1, ArrowHead, BlkSize],...
                'FaceColor', curr_color)
            if curr_split
                split_color = Colors(color_count+1, :);
                rectangle('Position', [curr_orbit(1 +BlkSize/2, j)-1, curr_orbit(1, j)-1, BlkSize/2, BlkSize/2],...
                    'FaceColor', split_color)
                rectangle('Position', [curr_orbit(1, j)-1, curr_orbit(1 +BlkSize/2, j)-1, BlkSize/2, BlkSize/2],...
                    'FaceColor', split_color)
                rectangle('Position', [curr_orbit(1+BlkSize/2, j)-1, nBlk*BlkSize, BlkSize, ArrowHead],...
                    'FaceColor', curr_color)
                rectangle('Position', [nBlk*BlkSize,  curr_orbit(1+BlkSize/2, j)-1, ArrowHead, BlkSize],...
                    'FaceColor', curr_color)
            end
        end
        color_count = color_count + 1 + splits(i);
    end
    xlim([0,N]);
    ylim([0,N]);
    axis square
    
    rectangle('Position', [nBlk*BlkSize, nBlk*BlkSize, ArrowHead, ArrowHead], 'FaceColor', 'k')
    set(gca, 'YDir','reverse')
    title('Permutation Orbits', 'FontSize', 22, 'interpreter', 'latex')
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
        subplot(2,2,2)
        %imagesc(MG)
        M_display = MG{k};
        M_display(M_display == 0) = NaN;
        imagesc(M_display)
        axis square
        title('Group-Invariant Matrix $X$', 'FontSize', 22, 'interpreter', 'latex')
        
        subplot(2,2,3)
        U_display = U;
        imagesc(U_display)
        axis square
        title('Block-Diagonalizer $P$', 'FontSize', 22, 'interpreter', 'latex')
        
        
        subplot(2,2,4)
        
        M_display = MG_diag{k};
        M_display(M_display == 0) = NaN;
        imagesc(M_display)
        axis square
        title('Block Diagonalization $X_k$', 'FontSize', 22,'interpreter', 'latex')
        pause(0.5);
        fname = strcat('sym_block_arrow-', num2str(k));
        %export_fig(fname, '-png', '-m', '4');
        %keyboard        
    end
end

