function D_new  = dd_shuffle(D, L)
    D_new  = zeros(size(D));
    for i = 1: size(D, 2)
         di = D(:, i);
         Di = v2m(di);
         Di_out = L * Di * L';
         D_new(:, i) = m2v(Di_out)
    end
end