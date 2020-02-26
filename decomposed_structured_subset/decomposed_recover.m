function [x_out] = decomposed_recover(x, info)
%DECOMPOSED_RECOVER Recover the optimum solution of a decomposed structured
%subset.
%
%Input:
%   x: vector of input that was previously 

dual = info.dual;

%maybe its sparse? I don't know, future work.
x_out = zeros(info.prev_var, 1);

%count_dd = info.count_dd;
%free and linear variables before conversion
Kold = info.prev_cone;
%x_out(1:count_dd) = x(1:(Kold.f + Kold.l));
x_out(1:(Kold.f + Kold.l)) = x([1:Kold.f, info.count_free + (1:Kold.l)]);

%diagonally dominant cone
count_dd = Kold.f + Kold.l;
count_non_dd = info.count_non_dd;
for i = 1:length(info.dd)    
    info_curr = info.dd{i};
    
    if dual
        %DD*
        x_out(info_curr.ind_orig) =  x(info_curr.ind_dual);
    else
        %DD
        range_curr = count_dd + (1:info_curr.num_var);

        x_dd = x(range_curr);

        x_out(info_curr.ind_orig) = info_curr.rays * x_dd;
        
    end
    count_dd = count_dd + info_curr.num_var;
end
count_dd = count_dd + info.num_rel;
%second order cone variables from original problem
if info.count_non_dd > 0    
    num_q = sum(Kold.q);
    x_out(Kold.f + Kold.l + (1:num_q)) = x(count_dd:count_non_dd);
end

%Factor Width Cones
count_non_dd = info.count_non_dd;
for i = 1:length(info.non_dd)
    
    info_curr = info.non_dd{i};
    range_curr = count_non_dd + (1:info_curr.num_var);
    if isfield(info_curr, 'ind_dual')
        %FW*
        x_out(info_curr.ind_orig) =  x(info_curr.ind_dual);
    elseif isfield(info_curr, 'Ech')
        %factor width
        
        x_reconstructed = accumarray(info_curr.Ech, x(range_curr));
        
        x_out(info_curr.ind_orig) = x_reconstructed;
    else
        %psd
        x_out(info_curr.ind_orig) = x(range_curr);
    end    

    count_non_dd = count_non_dd + info_curr.num_var;
end
