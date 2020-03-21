function [x, y] = convert_mosek2sedumi_var(sol,bardim)
%Convert the mosek-format variable into a sedumi-format variable
%input is res.sol.itr (solution) and prob.bardim (sizes of PSD blocks)
%x = sol.xx;
%does nott handle second order cones (K.q), add later
num_var = length(sol.xx) + sum(bardim.^2);
x = zeros(num_var, 1);

N_free_lin = length(sol.xx);
x(1:N_free_lin) = sol.xx;

count_sedumi = N_free_lin;
count_mosek = 0;
for i = 1:length(bardim)
    Ksi = bardim(i);
    
    vsi = (Ksi/2)*(Ksi+1);
    
    %precompute and cache, future improvement
    Mi = tri_indexer(Ksi);
    x_mosek_curr = sol.barx(count_mosek + (1:vsi));
    X_curr = x_mosek_curr(Mi);
    
    x(count_sedumi + (1:Ksi^2)) = reshape(X_curr, [], 1);
    count_sedumi = count_sedumi + Ksi^2;
    count_mosek = count_mosek + vsi;
end

y = sol.y;

end

