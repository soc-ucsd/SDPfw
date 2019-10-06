m = 20;
nBlk = 20;
BlkSize = 1;
ArrowHead = 1;

num_c = 2;
n = nBlk*BlkSize+ArrowHead;
%r = 7;
r = n;
rng(400, 'twister')

pars.fid = 0;
model.pars = pars;
[model.At,model.b,model.c,model.K, X_star] = blockArrowRank(m,nBlk,BlkSize,ArrowHead, r, num_c);

Ats = model.At(model.K.f+1:end, :);
cs = model.c(model.K.f+1:end, 1);        
SP = spones(spones(cs) + sparse(sum(spones(Ats),2)));  % vector of 1s and 0s
mask = reshape(SP, model.K.s, model.K.s);
% figure(2)
% spy(mask)

%chordal decomposition
%[cl, Ech, Jch, usedvars, s] = chordalDecomposition(model.At, model.c(:, 1), model.K);

cliques = {};
cl_ind = 0;
idx = [];

for i = 1:nBlk
    %cli = cl{1 ,1}.Elem (cl_ind + (1:cl{1,1}.NoElem(i)));
    %cl_ind = cl_ind + cl{1,1}.NoElem(i);
    cli =  [BlkSize*(i-1)+(1:BlkSize), ((-ArrowHead + 1):0) + model.K.s];
    cliques{i} = cli;
    idxi = reshape(cli +  model.K.s*(cli-1)', [],1 );
    idx = [idx; idxi];
    if i ==1
        idx0 = idxi;
    end
end

model_agler.K.f = model.K.f;
model_agler.K.l = 0;
model_agler.K.q = [];
model_agler.K.s = cellfun(@length, cliques);
model_agler.pars = pars;
numvar = sum(model_agler.K.s .^2);

model_agler.At = model.At(idx, :);
model_agler.b = model.b;
model_agler.c = model.c(idx, :);
% for i = 1:num_c
%     model_agler.c(:, i)  = model.c(idx, i);
% end


cone = 'dd';
[model_dd.At,model_dd.b,model_dd.c,model_dd.K, ~] = ...
    decomposed_subset(model.At, model.b, model.c, model.K, cone);
model_dd.pars = pars;

[model_cdd.At,model_cdd.b,model_cdd.c,model_cdd.K, ~] = ...
    decomposed_subset(model_agler.At, model_agler.b, model_agler.c, model_agler.K, cone);
model_cdd.pars = pars;

cone = 'sdd';
[model_sdd.At,model_sdd.b,model_sdd.c,model_sdd.K, ~] = ...
    decomposed_subset(model.At, model.b, model.c, model.K, cone);
model_sdd.pars = pars;

[model_csdd.At,model_csdd.b,model_csdd.c,model_csdd.K, ~] = ...
    decomposed_subset(model_agler.At, model_agler.b, model_agler.c, model_agler.K, cone);
model_csdd.pars = pars;

N = 30;
th = linspace(0, 2*pi, N);
PSD = draw_feasibility(model, th);
% CSDD = draw_feasibility(model_csdd, th);
% SDD = draw_feasibility(model_sdd, th);
% CDD = draw_feasibility(model_cdd, th);
% DD = draw_feasibility(model_dd, th);


%Plot the figures
C = linspecer(6);
figure(1)
clf
hold on
plot(PSD.a(PSD.conv), PSD.b(PSD.conv), 'k', 'linewidth', 2)
%plot(DD.a(DD.conv), DD.b(DD.conv),  'color', C (4, :), 'linewidth', 2)
%plot(SDD.a(SDD.conv), SDD.b(SDD.conv),  'color', C (1, :), 'linewidth', 2)
%plot(CDD.a(CDD.conv), CDD.b(CDD.conv),  'color', C (5, :), 'linewidth', 1)
%plot(CSDD.a(CSDD.conv), CSDD.b(CSDD.conv),  'color', C (5, :), 'linewidth', 1)
hold off

function out = draw_feasibility(model, th)
    N = length(th);
    out.a  = zeros(N, 1);
    out.b  = zeros(N, 1);
    out.Q  = cell(N, 1);
    
    c = model.c;
    c1 = model.c(:, 1);
    c2 = model.c(:, 2);
    
    for i = 1:N
        theta = th(i);
        
        
        
        c = c1*cos(theta) + c2*sin(theta);
        %c(2) = sin(theta);
        
        
        %there's a bug here, the chordal decomposition did not add new
        %equality constraints. How to add these in?
        [x,y,info] = sedumi(model.At, model.b, c, model.K, model.pars);
        K = model.K;
        out.a(i) = c1'*x;
        out.b(i) = c2'*x;
        out.Q{i} = x(K.f+1:end);
    end
    
    out.conv = convhull(out.a, out.b);
end