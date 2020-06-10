rng(20, 'twister');

% n_list = 20*(1:50);

n_list = 20*(1:50);

%Hout_list = zeros(length(n_list), 1);
% time_list_solve = zeros(length(n_list), 1);
% time_list_conv = zeros(length(n_list), 1);
G = 1;

outname = 'reference_dense3.mat';

QUIET = 0;
SOLVE = 0;
PLOT = 1;

if SOLVE

opts = sdpsettings('solver', 'mosek');

for k = 11:length(n_list)
% for k = 1:1
    n = n_list(k);
    m = ceil(n/10);
    d = ceil(n/10);
    Flag = 3;
    Sys = GenerateDynamics(G, [], [n], [m], [d], Flag);

    
    epsilon = 0.01;
    P = sdpvar(n);
    Constraint = [P - epsilon*eye(sum(n)) >= 0];
    gamma2 = sdpvar(1);
    Constraint2 = [[P*Sys.globalA+Sys.globalA'*P + Sys.globalC'*Sys.globalC, P*Sys.globalB + Sys.globalC'*Sys.globalD; 
                        Sys.globalB'*P + Sys.globalD'*Sys.globalC, Sys.globalD'*Sys.globalD-gamma2*eye(sum(m))] + epsilon*eye(sum(n)+sum(m)) <= 0];
    Constraint = [Constraint, Constraint2, gamma2 >= 0];
    Cost = gamma2;
    
    sol = optimize(Constraint, Cost, opts);
%     flag_str = 'Hinf0';
    
%     [Hout, time_solve, sdp_opt] = sea_star_reference(model, QUIET);
    Hout = sqrt(value(gamma2));
    Hout_list(k) = Hout;
    time_list_solve(k) = sol.solvertime;
    time_list_conv(k) = sol.yalmiptime;
    
    save(outname, 'Hout_list', 'time_list_conv', 'time_list_solve', 'k', 'n')
    
    [k, Hout, sol.solvertime]
    
%     keyboard
    
end

end

if PLOT
   figure(1)
   clf
   k_crit = max(find(time_list_solve ~= 0));
   
   plot(n_list(1:k_crit), time_list_solve(1:k_crit)/60, 'linewidth', 4);
   
   title('H infinity Size of System vs. Computation Time')
   xlabel('Number of States')
   ylabel('Time (minutes)')
   
end