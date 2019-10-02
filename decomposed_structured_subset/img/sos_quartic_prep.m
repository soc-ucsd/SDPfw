load('quartic_data_sparse.mat');

pars.fid = 0;

model.At = At;
model.b  = b;
model.c  = c;
model.K  = K;
model.pars = pars;

K2 = K;
cliques = {[1,3], [2,3]};
K2.s = cellfun(@length, cliques);

%i1 = K.f + [1, 3, 7, 9];
%i2 = K.f + [5, 6, 8, 9];
%idx = [1:K.f, i1, i2];
cons = [1, 2, 3, 5];
idx = 1:K.f;

numvar = K.f + sum(model2.K.s .^2);
for i = 1:length(cliques)
    cli = cliques{i};
    idxi = reshape(cli +  K.s*(cli-1)' + K.f, 1, []);
    idx = [idx idxi];
end

% kills  number 4
At2 =  At(idx, cons);
c2 =  c(idx);
b2  =  b(cons);

model2.At = At2;
model2.b  = b2;
model2.c  = c2;
model2.K  = K2;
model2.pars = pars;

save('quartic_process0.mat', 'model', 'model2');