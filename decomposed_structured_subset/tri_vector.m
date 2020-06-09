function v_out = tri_vector(A, two_scale, lower)
    %vectorize upper triangular part of matrix    

    %scale off-diagonal elements by 2
    if nargin < 2
        two_scale = 0;
    end
    
    if nargin < 3
        lower = 1;
    end
    
    if two_scale
        %A = 2*A - diag(diag(A));
        A = spdiags(diag(A), 0, 2*A);
    end                            
    
    %mask = triu(true(size(A)));
    if lower
        mask = tril(true(size(A)));
    else
        mask = triu(true(size(A)));
    end
    
    
    v_out = A(mask);
end
        