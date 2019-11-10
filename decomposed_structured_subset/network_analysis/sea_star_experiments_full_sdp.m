%load('sea_star_Hinf0_verylarge.mat', 'model')
load('sea_star_H2_tiny.mat', 'model')

%[Houtdd, Resdd, tsdd, tcdd] = run_model_star(model, 'dd', 1);
%cost_table= latex(vpa(sym(time_solve),4));
USE_DD = 0;

%yalmip outputs LMI
parCoLO.domain    = 0;  % dConvCliqueTree  ---> equalities 
parCoLO.range     = 0;   % rConvMatDecomp   ---> equalities 
parCoLO.EQorLMI   = 2; % CoLOtoEQform     ---> LMI standard form
parCoLO.SDPsolver = []; % CoLOtoEQform     ---> LMI standard form       
%parCoLO.SDPsolver = 'sedumi'; % CoLOtoEQform     ---> LMI standard form       

parCoLO.quiet     = 1; % Some peace and quiet       
J.f = length(model.b);

[~,~,~,cliqueDomain,cliqueRange,LOP] = sparseCoLO(model.A',model.b,model.C,model.K,J,parCoLO); 


tic;
if USE_DD
    [A, b, c, K, ~] = decomposed_subset(-LOP.A,-LOP.C,-LOP.b,LOP.J, 'dd');
    %[A, b, c, K, ~] = decomposed_subset(model.A,model.b,model.C,model.K, 'dd');
    prob = sedumi2mosek(A, b, c, K);
else
    prob = sedumi2mosek(-LOP.A, -LOP.c, -LOP.b, LOP.J);    
end
time_convert = toc;
    
prob2 = sedumi2mosek(model.A, model.b, model.C, model.K);
% Set log level (integer parameter)
%param.MSK_IPAR_LOG = 1;
% Select interior-point optimizer... (integer parameter)
param.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-6;
if USE_DD
    param.MSK_IPAR_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
    % ... without basis identification (integer parameter)
param.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';            
end
tic;
[r,res] = mosekopt('minimize' ,prob, param);
time_solve = toc;
%[r,res] = mosekopt('minimize echo(0)',prob, param);

if  strcmp(res.sol.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE')        
    cost = res.sol.itr.pobjval;
else
    cost = NaN;
end


Hout = sqrt(cost);