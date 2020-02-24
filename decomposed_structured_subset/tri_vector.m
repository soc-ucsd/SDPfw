function v_out = tri_vector(A, two_scale)
    %vectorize upper triangular part of matrix    

    %scale off-diagonal elements by 2
    if nargin < 2
        two_scale = 0;
    end
    
    if two_scale
        %A = 2*A - diag(diag(A));
        A = spdiags(diag(A), 0, 2*A);
    end                            
    
    %mask = triu(true(size(A)));
    mask = tril(true(size(A)));
    v_out = A(mask);
end
        