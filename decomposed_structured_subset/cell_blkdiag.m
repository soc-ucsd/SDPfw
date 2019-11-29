function [A_diag] = cell_blkdiag(A)
%CELL_BLKDIAG Given a cell filled with matrices A, generate a sparse matrix
%with A along the main diagonal.
%this should really be default

m = cellfun(@(x) size(x,1), A);
n = cellfun(@(x) size(x,2), A);

%mc = cumsum([0, m]+1);
%nc = cumsum([0, n]+1);

i = [];
j = [];
v = [];

m_count = 0;
n_count = 0;
for B = 1:length(A)
    
    %row_ind = (1:m(B))'*ones(1, n(B)) + m_count;
    %col_ind = ones(m(B),1)*(1:n(B))   + n_count;
    
    [ii, jj, vv] = find(A{B});
    i = [i; ii + m_count];
    j = [j; jj + n_count];
    v = [v; vv];
    
    m_count = m_count + m(B);
    n_count = n_count + n(B);
end

A_diag = sparse(i, j, v, sum(m), sum(n));
end

