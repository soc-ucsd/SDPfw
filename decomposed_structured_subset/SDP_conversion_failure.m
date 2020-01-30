% [model.A, model.b, model.c, model.K] = convert_mosek2sedumi(prob0);
% prob1 = convert_sedumi2mosek(model.A, model.b, model.c, model.K);
% [model1.A, model1.b, model1.c, model1.K] = convert_mosek2sedumi(prob1);

CONSTRAINED = 1;
OPTIMIZE = 0;

% if CONSTRAINED
%     %Constrained (optimum of 535.98)
%     read_file = 'read(LR18_cons_4.task.gz)';
% else
%     %Unconstrained (optimum of 12.47)
%     read_file = 'read(LR18_uncons_3.task.gz)';
% end

%[r, res0] = mosekopt(read_file);

%res_cons_trans = load('LR_120.mat', 'res_cons_trans');
res_cons_trans = load('LR120_box_1_2.mat', 'res_cons_trans');
res0 = res_cons_trans.res_cons_trans;
%Initial Mosek
prob0 = res0.prob;

prob0.c = -prob0.c;
prob0.barc.val = -prob0.barc.val;
prob0 = rmfield(prob0, 'cfix');
prob0 = rmfield(prob0, 'names');
prob0 = rmfield(prob0, 'objsense');

prob0.barc.subj = prob0.barc.subj';
prob0.barc.subk = prob0.barc.subk';
prob0.barc.subl = prob0.barc.subl';
prob0.barc.val  = prob0.barc.val';
prob0.bara.subi = prob0.bara.subi';
prob0.bara.subj = prob0.bara.subj';
prob0.bara.subk = prob0.bara.subk';
prob0.bara.subl = prob0.bara.subl';
prob0.bara.val  = prob0.bara.val';
prob0.blx = prob0.blx';
prob0.bux = prob0.bux';

prob0.a = [prob0.a(:, end), prob0.a(:, 1:(end - 1))];
prob0.blx(1) = -inf;
prob0.blx(end) = 0;
cnew = prob0.c(1);
prob0.c(1) = prob0.c(end);
prob0.c(end) = prob0.c(1);

res_cons_trans = prob0;
model_cons_trans = model;
save('LR120_box_1_2.mat', 'res_cons_trans', 'model_cons_trans');

%Mosek -> Sedumi
[model.A, model.b, model.c, model.K] = ...
convert_mosek2sedumi(prob0);


%Mosek -> Sedumi -> Mosek
prob1 = convert_sedumi2mosek(model.A, model.b, model.c, model.K);


%Mosek -> Sedumi -> Mosek -> Sedumi
[model1.A, model1.b, model1.c, model1.K] = convert_mosek2sedumi(prob1);

%

% cost_mosek2 = res.sol.itr.pobjval;

if OPTIMIZE
    [r, res] = mosekopt('minimize', prob0);
    cost_mosek0 = res.sol.itr.pobjval;

    [x, y, info] = sedumi(model.A, model.b, model.c, model.K);
    cost_sedumi = model.c'*x;

    [r, res] = mosekopt('minimize', prob1);
    cost_mosek1 = res.sol.itr.pobjval;

    [x, y, info] = sedumi(model1.A, model1.b, model1.c, model1.K);
    cost_sedumi1 = model1.c'*x;
end