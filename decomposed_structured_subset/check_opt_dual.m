function [sdp_optimal,cone_valid] = check_opt_dual(z_opt,K, cones)
%CHECK_OPT_DUAL Is a decomposed structured subset solution dual feasible?
%               Does it have the same cost as the true SDP solution?
%
%Input:
%   z:     dual optimal solution: from K* optimization, or output from
%          recover_opt_dual
%   K:     sedumi cones from original SDP
%   cones: list of structured subsets
%
%Output:
%   sdp_optimal:    Does z solve the SDP?
%   cone_valid:     Which cones are valid, and which cones are violated
%                   (not PSD)?
%tolerance of 1e-6, edit later?

    numPSD = length(K.s);
    %set up cones array
    if ~iscell(cones)
        cone_str = cones;
        cones = cell(numPSD, 1);
        for i = 1:numPSD
            cones{i} = cone_str;
        end
    end
    
    cone_valid = zeros(numPSD, 1);
    
    Count = K.f + K.l + sum(K.q);
    for PSDind = 1:numPSD
        Ksi = K.s(PSDind);
        curr_range = Count + (1:Ksi^2);
        
        if strcmp('psd', cones{PSDind})
            cone_valid(PSDind) = 1;
        else
            Z_curr = reshape(z_opt(curr_range), Ksi, Ksi);
            d = eig(Z_curr);
            cone_valid(PSDind) = all(d >= -1e-6);
        end
            
        
        Count = Count + Ksi^2;
    end
       
    sdp_optimal = all(cone_valid);
    
end

