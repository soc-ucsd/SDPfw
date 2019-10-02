function L_out = basis_shift(X)
    delta = 1e-5;
    X_inv = inv(delta*eye(length(X)) + X);
    X_shift = X + 0.1*delta*X_inv;
    %x_dd_shift = m2v(X_dd_shift);
    L_out = chol(X_shift);
end