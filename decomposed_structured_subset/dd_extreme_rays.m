function rays =  dd_extreme_rays(N)

L = N + 2*N*(N-1)/2;
M = reshape(1:(N^2) , N, N);

%setup

%diagonal
i = diag(M)';
j = 1:N;
v = ones(size(i));


curr_ind = N+1;
for k1 = 1:N
    for k2 = (k1+1):N
        j_pos = [1 1 1 1]*curr_ind;
        i_pos = reshape(M([k1, k2], [k1,k2]),1, 4);
        v_pos = [1 1 1 1];
        v_neg = [1 -1 -1 1];
        
        i = [i i_pos i_pos];
        j = [j j_pos (j_pos + 1)];       
        v = [v v_pos v_neg];
        
        curr_ind = curr_ind + 2;
    end    
end


rays = sparse(i, j, v , N^2, L);
end