function ec  = eig_cdd(CDD)
    ec = zeros(4,  size(CDD, 2));
    for i = 1:size(CDD, 2)
        X_cdd = v2m(CDD(: , i));    
        X1 = X_cdd([1,2], [1,2]);
        X2 = X_cdd([2,3], [2,3]);
        ec([1, 2], i) = eig(X1);
        ec([3, 4], i) = eig(X1);
    end
end