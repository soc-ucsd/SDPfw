function A_out = tri_matrix(A, two_scale)

[m, N2] = size(A);
N = sqrt( N2 );
N_out =  N* (N+1)/2;
ai = [];
aj = [];
av = [];
for i = 1:m
    a_curr = A(i, :);
    A_curr = reshape(a_curr, N, N);
    v_curr = tri_vector(A_curr, two_scale);
    
    
    I = find(v_curr);
    count = length(I);
    ai = [ai; i*ones(count, 1)];
    aj = [aj; I];
    av = [av; v_curr(I)];           
end

A_out =  sparse(ai, aj, av, m, N_out);