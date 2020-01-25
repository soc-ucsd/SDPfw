function x_out = mon_basis(x, alpha)
    %x^alpha because I dont know how yalmip does it
    %initialize into vector of sdpvar
    N_alpha = size(alpha, 2);
    x_out = x(1)*ones(N_alpha, 1);
    
    for i = 1:N_alpha
        alpha_curr = alpha(:, i);
        x_out(i) = prod(x.^alpha_curr);
    end    
end