%load('quartic_data_3.mat');
load('quartic_data_3m.mat');

pars.fid = 0;

model.At = At;
model.b  = b;
model.c  = c;
model.K  = K;
model.pars = pars;

%important constraints, i'm not sure about 8
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
%cons = [1:13, 15];
cons = [1:12 15]

K.l = 0;
K.q = [];

[cl, Ech, Jch, usedvars, s] = chordalDecomposition(model.At(:, cons), model.c, K);

cliques = {};
cl_ind = 0;
for i = 1:(cl{1,1}.NoC)
    cli = cl{1 ,1}.Elem (cl_ind + (1:cl{1,1}.NoElem(i)));
    cl_ind = cl_ind + cl{1,1}.NoElem(i);
    cliques{i} =  cli;
end

Ats = At(K.f+1:end, cons);
cs = c(K.f+1:end);        
SP = spones(spones(cs) + sparse(sum(spones(Ats),2)));  % vector of 1s and 0s
mask = reshape(SP, K.s, K.s);


%spy(mask)
%G = graph(mask);
%plot(G)


model2.K.f = K.f;
model2.K.l = 0;
model2.K.q = [];
model2.K.s = cellfun(@length, cliques);
model2.pars = pars;
 
numvar = K.f + sum(model2.K.s .^2);

idx = (1:K.f)';
for i = 1:length(cliques)
    cli = cliques{i};
    idxi = reshape(cli +  K.s*(cli-1)' + K.f, [],1 );
    idx = [idx; idxi];
end


model2.At = At(idx, cons );
model2.b = b(cons);
model2.c  = c(idx);
% K2 = K;
% K2.s = [2, 2];
% i1 = K.f + [1, 3, 7, 9];
% i2 = K.f + [5, 6, 8, 9];
% is = [1:K.f, i1, i2];
% js = [1, 2, 3, 5];
%  
% % kills  number 4
% At2 =  At(is, js);
% c2 =  c(is);
% b2  =  b(js);
% 
% model2.At = At2;
% model2.b  = b2;
% model2.c  = c2;
% model2.K  = K2;
% model2.pars = pars;
% 
save('quartic_process_3.mat', 'model', 'model2');