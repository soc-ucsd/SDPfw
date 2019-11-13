%load('sea_star_H2_tiny.mat', 'Sys', 'n', 'd', 'm')
load('sea_star_H2_micro.mat', 'Sys', 'n', 'd', 'm')
%load('sea_star_Hinf0_small.mat', 'Sys', 'n', 'd', 'm')


epsilon = 0.01;
P = sdpvar(sum(n), sum(n));
P_reg = P - epsilon*eye(sum(n));
N = sum(n);

num_input  = size(Sys.globalB, 2);
num_output = size(Sys.globalC, 1);
num_state  = size(Sys.globalA, 2);

Sys.globalD = sparse(num_output, num_input);

gamma2 = sdpvar(1);
Bounded_Real = -[P*Sys.globalA+Sys.globalA'*P + Sys.globalC'*Sys.globalC, P*Sys.globalB + Sys.globalC'*Sys.globalD; 
                Sys.globalB'*P + Sys.globalD'*Sys.globalC, Sys.globalD'*Sys.globalD-gamma2*eye(sum(m))]  - epsilon*eye(sum(n)+sum(m));
NB = length(Bounded_Real);

Constraint = [P_reg >= 0, Bounded_Real + epsilon*eye(NB) >= 0];

Cost = gamma2;

opts = sdpsettings('solver', 'mosek', 'verbose', 1);

%sol = optimize(Constraint, Cost, opts);

%gamma = sqrt(value(gamma2));

model = export(Constraint, Cost, sdpsettings('solver', 'sedumi'));
% 

% Couldn't get dual optimization to work on the LMI.
% Think I'll try doing some stuff manually.
% parCoLO.domain    = 2;  % dConvCliqueTree  ---> equalities 
% parCoLO.range     = 1;   % rConvMatDecomp   ---> equalities 
% parCoLO.EQorLMI   = 1; % CoLOtoEQform     ---> LMI standard form
% parCoLO.SDPsolver = []; % CoLOtoEQform     ---> LMI standard form       
% %parCoLO.SDPsolver = 'sedumi'; % CoLOtoEQform     ---> LMI standard form       
% 
% parCoLO.quiet     = 1; % Some peace and quiet       
% J.f = length(model.b);
% 
% [~,~,~,cliqueDomain,cliqueRange,LOP] = sparseCoLO(model.A',model.b,model.C,model.K,J,parCoLO); 
% 
% 
% [model_dd_star.At, model_dd_star.b, model_dd_star.c, model_dd_star.K, model_dd_star.info] = ...
% dd_star_convert(LOP.A,LOP.b,LOP.c,LOP.K);    
% %dd_convert(model.At,model.b,model.c,model.K);
% 
% [x, y, info] = sedumi(model_dd_star.At, model_dd_star.b, model_dd_star.c, model_dd_star.K);