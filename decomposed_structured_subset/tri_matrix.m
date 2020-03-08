function A_out = tri_matrix(A, two_scale)

[m, N2] = size(A);
N = sqrt( N2 );
N_out =  N* (N+1)/2;
mask = tril(true(N));
mask_shape = reshape(mask, [], 1) ;

A_out =  A(:, mask_shape);

if two_scale
    maskw = 2*mask - eye(N);
    maskw_shape = reshape(maskw(mask), 1, []); 
    
    %Aw  = A_out .* maskw_shape;
    A_out = bsxfun(@times, A_out, maskw_shape);
end
