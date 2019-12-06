function [MG] = reynolds_perm(M, orbit_ind, BlkSize, arrow_ind)
%REYNOLDS_PERM Performs the reynolds operator of projection onto an
% permutation invariant subspace. orbit_ind are the orbits of indices. 
%   Detailed explanation goes here

    MG = M;
    ArrowHead = length(arrow_ind);
    for i = 1:length(orbit_ind)
        curr_orbit = orbit_ind{i};
        num_orbits = size(curr_orbit, 2);

        Diag_block = zeros(BlkSize);
        Arrow_block = zeros(BlkSize, ArrowHead);

        %average over orbits
        %w =  num_orbits;
        w = 1;
        %w = sqrt(num_orbits);
        for j = 1:num_orbits
            curr_ind = curr_orbit(:, j);
            Diag_block = Diag_block + M(curr_ind, curr_ind)/w;
            Arrow_block = Arrow_block + M(curr_ind, arrow_ind)/w;
        end

        %assign back into matrix
        for j = 1:num_orbits
            curr_ind = curr_orbit(:, j);
            MG(curr_ind,  curr_ind) = Diag_block;
            MG(curr_ind, arrow_ind) = Arrow_block;
            MG(arrow_ind, curr_ind) = Arrow_block';
        end

        %Block Diagonalization
        %raw_orbit = sparse(orbits{i}, orbits{i}, ones(size(orbits{i})), nBlk, nBlk);
        %U_dct(orbits{i}, orbits{i}) = dctmtx(num_orbits);    
    end


end

