function [x_out] = decomposed_recover(x, info)
%DECOMPOSED_RECOVER Recover the optimum solution of a decomposed structured
%subset.
%
%Input:
%   x: vector of input that was previously 


%maybe its sparse? I don't know, future work.
x_out = zeros(info.prev_var, 1);

count_dd = info.count_dd;
%free and linear variables before conversion
x_out(1:count_dd) = x(1:count_dd);


%diagonally dominant cone
for i = 1:length(info.dd)
    info_curr = info.dd{i};
    
    range_curr = count_dd + (1:info_curr.num_var);
    
    x_dd = x(range_curr);
    
    x_out(info_curr.ind_orig) = info_curr.rays * x_dd;
    count_dd = count_dd + info_curr.num_var;
end

%second order cone variables from original problem
if info.count_non_dd > 0
    num_q = info.count_non_dd - count_dd;
    x_out(info.count_dd + (1:num_q)) = x(count_dd:info.count_non_dd);
end

%Factor Width Cones
count_non_dd = info.count_non_dd;
for i = 1:length(info.non_dd)
    
    info_curr = info.non_dd{i};
    range_curr = count_non_dd + (1:info_curr.num_var);
    if isfield(info_curr, 'Ech')
        %factor width
        
        x_reconstructed = accumarray(info_curr.Ech, x(range_curr));
        
        x_out(info_curr.ind_orig) = x_reconstructed;
    else
        %psd
        x_out(info_curr.ind_orig) = x(range_curr);
    end    

    count_non_dd = count_non_dd + info_curr.num_var;
end
