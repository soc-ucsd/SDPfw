%load('quartic_data_3.mat');
%load('quartic_data_3m.mat');
load('quartic_data_3m.mat');

pars.fid = 0;

model.At = At;
model.b  = b;
model.c  = c;
model.K  = K;
model.pars = pars;

%original constraints
%cons = [1, 3,  5, 6, 7, 8, 9, 11, 12, 15];

%cliques = {[1,2,5] ,[1,2,3], [2,3,4] , [3,4,6]};

%i missed 7
%cliques = {[2, 3,4, 6],[1,2,3],[1,2,5],[2,4,5]}
%cliques = {[2, 3,4, 6],[2,3,4,5], [1,3,5]};

%cons = [1, 2,     3, 5, 6, 9, 12, 15]  ;
%cliques = {[1, 3], [2, 3], [ 3 , 6], [4, 6] , [5, 6] };
%cliques_11 = {[1 , 3, 5], [3, 5, 6], [2, 3, 6] , [2, 6, 4]};
%Ats = At(K.f+1:end, :);
%cons = [1:13 15];
%cons = 1:15;

%b nonzero at 1, 5, 11, 12, 15
%quartic_data_3m
cons = [1, 2, 3, 4, 5:12 15];


%zero sparsity:
%2:  [32, 37]
%4:  [19, 29]
%13: [6, 21]
%14: [4,9]
%strip out 2 and 4?
%cons = [1, 3, 5:12, 15];
%cons = [1, 3, 5:13, 15];
%cons = [1, 2, 3, 4, 5:12, 15];

%add in {1,2}, {3,5}
%cons = [1, 3, 4, 5:12, 13, 15];
%zero forcing constraints


K.l = 0;
K.q = [];

%I found the bug. Wrong cliques being found, didn't index out K.f
[cl, Ech, Jch, usedvars, s] = chordalDecomposition(model.At(K.f+1:end, cons), model.c(K.f+1:end), K);

cliques = {};
cl_ind = 0;
idx = (1:K.f)';

for i = 1:(cl{1,1}.NoC)
    cli = cl{1 ,1}.Elem (cl_ind + (1:cl{1,1}.NoElem(i)));
    cl_ind = cl_ind + cl{1,1}.NoElem(i);
    cliques{i} =  cli;
    idxi = reshape(cli +  K.s*(cli-1)' + K.f, [],1 );
    idx = [idx; idxi];
    if i ==1
        idx0 = idxi;
    end
end

Ats = At(K.f+1:end, cons);
cs = c(K.f+1:end);        
SP = spones(spones(cs) + sparse(sum(spones(Ats),2)));  % vector of 1s and 0s
mask = reshape(SP, K.s, K.s);
% 
% 
% spy(mask)
% G = graph(mask);
% plot(G)


model_agler.K.f = K.f;
model_agler.K.l = 0;
model_agler.K.q = [];
model_agler.K.s = cellfun(@length, cliques);
model_agler.pars = pars;
 
numvar = K.f + sum(model_agler.K.s .^2);


model_agler.At = At(idx, cons );
model_agler.b = b(cons);
model_agler.c  = c(idx);


%Now do the decompositions


%SDD
Partition = Inf;
opts.bfw  = 1;
opts.nop  = Partition;
opts.socp = 1;   % second-order cone constraints

model_sdd.pars = pars;
[model_sdd.At, model_sdd.b, model_sdd.c, model_sdd.K, info_fw_decomp] = ...
    factorwidth(model.At',model.b,model.c,model.K,opts);

%DD
model_dd.pars = pars;
[model_dd.At, model_dd.b, model_dd.c, model_dd.K, model_dd.info] = ...
    dd_convert(model.At',model.b,model.c,model.K);



%block factor width over cliques
model_csdd.pars = pars;
[model_csdd.At, model_csdd.b, model_csdd.c, model_csdd.K, info_cl_fw_decomp] = ...
    factorwidth(model_agler.At',model_agler.b,model_agler.c,model_agler.K,opts);

%SDD on 5-size clique
At_large = At(idxi, cons);
c_large = c(idxi);
K_large.f  = 0;
K_large.l  = 0;
K_large.q  = [];
K_large.s  = model_agler.K.s(end);

[At_large_fw, ~, c_large_fw, K_large_fw, ~] =  ...
    factorwidth(At_large, model_agler.b, c_large, K_large, opts);

model_large.At = [model_agler.At(1:K.f, :); At_large_fw'; model_agler.At(K.f+(1:sum(model_agler.K.s(1:end-1).^2)), :)];
model_large.c  = [model_agler.c(1:K.f, :); c_large_fw; model_agler.c(K.f+(1:sum(model_agler.K.s(1:end-1).^2)), :)];
model_large.b = model_agler.b;
model_large.K = model_agler.K;
model_large.K.q = K_large_fw.q;
model_large.K.s = model_agler.K.s(1:end-1);
model_large.pars = pars;

%DD on 5-size clique
[At_large_dd, ~, c_large_dd, K_large_dd, ~] = ...
    dd_convert(At_large, model_agler.b, c_large, K_large);

model_large_dd.At = [model_agler.At(1:K.f, :); At_large_dd'; model_agler.At(K.f+(1:sum(model_agler.K.s(1:end-1).^2)), :)];
model_large_dd.c  = [model_agler.c(1:K.f, :); c_large_dd; model_agler.c(K.f+(1:sum(model_agler.K.s(1:end-1).^2)), :)];
model_large_dd.b = model_agler.b;
model_large_dd.K = model_agler.K;
model_large_dd.K.l = model_agler.K.l + K_large_dd.l;
model_large_dd.K.s = model_agler.K.s(1:end-1);
model_large_dd.pars = pars;

%SDD on 4-size clique
At_small = At(idx0, cons);
c_small = c(idx0);
K_small.f  = 0;
K_small.l  = 0;
K_small.q  = [];
K_small.s  = model_agler.K.s(1);

[At_small_fw, ~, c_small_fw, K_small_fw, ~] =  ...
    factorwidth(At_small, model_agler.b, c_small, K_small, opts);

%Kbsmall = sum(model_agler.K.s(2:end).^2);
Kbsmall = K.f + model_agler.K.s(1)^2 + 1;
model_small.At = [model_agler.At(1:K.f, :); At_small_fw'; model_agler.At(Kbsmall:end, :)];
model_small.c  = [model_agler.c(1:K.f, :); c_small_fw; model_agler.c(Kbsmall:end, :)];
model_small.b = model_agler.b;
model_small.K = model_agler.K;
model_small.K.q = K_small_fw.q;
model_small.K.s = model_agler.K.s(2:end);

model_small.pars = pars;
%
save('quartic_process_3m.mat', 'model', 'model_agler', 'model_sdd', 'model_csdd',  'model_large', 'model_small');