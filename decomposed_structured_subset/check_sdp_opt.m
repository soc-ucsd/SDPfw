function [sdp_opt, cone_valid] = check_sdp_opt(x_opt, y_opt, A, b, c, K, cones, dual)
%CHECK_SDP_OPT Does the decomposed structured subset solve the SDP?
%   Calls check_opt_dual after forming the dual optimal solution from KKT
%Input:
%   x_opt:   primal optimal solution from Sedumi, could be over K or K*
%            Assumes that x_opt is indeed optimal over decomposed structured subsets           
%   y_opt:   multipliers on equality constraints from Sedumi
%   A, c, K: sedumi cones from original SDP
%   cones:   list of structured subsets for each entry of K.s 
%   dual:    Optimization over K or K* (0 or 1)?
%
%Output:
%   sdp_optimal:    Does z solve the SDP?
%   cone_valid:     Which cones are valid, and which cones are violated
%                   (not PSD)?


if dual
    z_opt = x_opt;
else
    %recover the dual solution
    if size(A, 1) == length(c)
        A = A';
    end    

    z_opt = c - A'*y_opt(1:length(b));
end 

[sdp_opt,cone_valid] = check_opt_dual(z_opt,K, cones);

end

