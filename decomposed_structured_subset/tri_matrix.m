function A_out = tri_matrix(A, two_scale, lower)

if nargin < 3
    lower = 1;
end

if nargin < 2
    two_scale = 0;
end

[m, N2] = size(A);
N = sqrt( N2 );
N_out =  N* (N+1)/2;
if lower
    mask = tril(true(N));
else
    mask = triu(true(N));
end
mask_shape = reshape(mask, [], 1) ;

A_out =  A(:, mask_shape);

if two_scale
    maskw = 2*mask - eye(N);
    maskw_shape = reshape(maskw(mask), 1, []); 
    
    %Aw  = A_out .* maskw_shape;
    A_out = bsxfun(@times, A_out, maskw_shape);
end
