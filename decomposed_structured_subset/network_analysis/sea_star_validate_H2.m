load('sea_star_H2_small.mat')

%gamma = [    2.4331    1.6710    0.2955    0.9732];
%gamma_med = [5.9351    2.7972    1.4517    0.0974    0.9592];

epsilon = 0.01;
P = [];
N = length(Sys.A);
for i = 1:N
    P = blkdiag(P,sdpvar(n(i)));
end


i = 4;
%gamma = gamma_med(i);
P = [];
for i = 1:N
    P = blkdiag(P,sdpvar(n(i)));
end

Constraint = [P - epsilon*eye(sum(n)) >= 0];
Constraint = [Constraint, Sys.globalA*P + P*Sys.globalA' + Sys.globalB*Sys.globalB' + epsilon*eye(sum(n)) <= 0];
Cost = trace(Sys.globalC*P*Sys.globalC');


opts          = sdpsettings('verbose',1,'solver','sparsecolo','sparsecolo.domain',1,'sparsecolo.range',0,'sparsecolo.EQorLMI',1,'sparsecolo.SDPsolver','sedumi');
solSparseCoLO = optimize(Constraint,Cost,opts);
%H2colo = sqrt(value(Cost));
PV = value(P);
check(Constraint)

plot(log10(eig(PV)))