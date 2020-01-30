% [model.A, model.b, model.c, model.K] = convert_mosek2sedumi(prob0);
% prob1 = convert_sedumi2mosek(model.A, model.b, model.c, model.K);
% [model1.A, model1.b, model1.c, model1.K] = convert_mosek2sedumi(prob1);

CONSTRAINED = 0;

if CONSTRAINED
    %Constrained (optimum of 535.98)
    read_file = 'read(LR18_cons_4.task.gz)';
else
    %Unconstrained (optimum of 12.47)
    read_file = 'read(LR18_uncons_3.task.gz)';
end

[r, res0] = mosekopt(read_file);


%Initial Mosek
prob0 = res0.prob;
[r, res] = mosekopt('maximize', prob0);
cost_mosek0 = res.sol.itr.pobjval;

%Mosek -> Sedumi
[model.A, model.b, model.c, model.K] = ...
 convert_mosek2sedumi(res0.prob);

[x, y, info] = sedumi(model.A, model.b, -model.c, model.K);
cost_sedumi = model.c'*x;


%Mosek -> Sedumi -> Mosek
prob1 = convert_sedumi2mosek(model.A, model.b, model.c, model.K);
[r, res] = mosekopt('maximize', prob1);
cost_mosek1 = res.sol.itr.pobjval;

%Mosek -> Sedumi -> Mosek -> Sedumi
[model1.A, model1.b, model1.c, model1.K] = convert_mosek2sedumi(prob1);
[x, y, info] = sedumi(model.A, model.b, -model.c, model.K);
cost_sedumi1 = model.c'*x;
%

% cost_mosek2 = res.sol.itr.pobjval;
