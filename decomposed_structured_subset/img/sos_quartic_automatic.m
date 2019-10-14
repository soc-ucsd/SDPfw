load('quartic_data_3m.mat');

pars.fid = 0;

model.At = At;
model.b  = b;
model.c  = c;
model.K  = K;
model.pars = pars;

cons = [1, 2, 3, 4, 5:12 15];

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

model_agler.K.f = K.f;
model_agler.K.l = 0;
model_agler.K.q = [];
model_agler.K.s = cellfun(@length, cliques);
model_agler.pars = pars;
 
numvar = K.f + sum(model_agler.K.s .^2);


model_agler.At = At(idx, cons );
model_agler.b = b(cons);
model_agler.c  = c(idx);

%Initial feasible point
model_feas_0 = model_agler;
model_feas_0.At = [model_feas_0.At, [eye(2); zeros(size(model_feas_0.At,1)-2, 2)]];
model_feas_0.b = [model_feas_0.b; 0; 0];

[x0, y0, info0] = sedumi(model_feas_0.At, model_feas_0.b, model_feas_0.c, model_feas_0.K, model.pars);

model_change = basis_change(x0, model_feas_0);


%Now do the decompositions

%same decompositions
%model_new.pars = pars;
%cone = 'sdd';

cones = {'dd', 'sdd', 'psd'};
LC = length(cones);
OUT = cell(LC, LC);

N = 200;
th = linspace(0, 2*pi, N);

OUT_0 = draw_feasibility(model, 0, th);

for i = 1:LC
    for j = 1:LC
        cone_curr = {cones{i}, cones{j}};
        OUT{i, j} = draw_feasibility(model_change, cone_curr, th);
    end
end


figure(1)
clf
C = linspecer(6);
for i = 1:LC
    for j= 1:LC
        subplot(LC, LC, i + LC*(j-1))
        hold on
        
        plot(OUT_0.a(OUT_0.conv), OUT_0.b(OUT_0.conv), 'color', 'k', 'linewidth', 2)
        
        OUT_curr = OUT{i, j};
        plot(OUT_curr.a(OUT_curr.conv), OUT_curr.b(OUT_curr.conv), 'color', C(1, :), 'linewidth', 2)
        
        plot(0, 0, 'xk', 'markersize', 12)
        
        model_out = basis_change(OOUT_CURR.x, model);
        
        title(strcat('$K_4=', upper(cones{i}), ',$  $K_5=', upper(cones{j}), '$'), 'Interpreter', 'latex', 'FontSize', 18)
        xlabel('a')
        ylabel('b')
        hold off

    end
end


function out = draw_feasibility(model, cone, th)
    N = length(th);
    out.a  = zeros(N, 1);
    out.b  = zeros(N, 1);
    out.Q  = cell(N, 1);
    
    [model_new.At, model_new.b, model_new.c, model_new.K, model_new.info] = ...
        decomposed_subset(model.At',model.b,model.c,model.K,cone);
        model_new.pars = model.pars;

    
    c = model_new.c;
    
    for i = 1:N
        theta = th(i);
        
        c(1) = cos(theta);
        c(2) = sin(theta);
        
        [x,y,info] = sedumi(model_new.At, model_new.b, c, model_new.K, model_new.pars);
        K = model.K;
        out.a(i) = x(1);
        out.b(i) = x(2);
        out.Q{i} = x(K.f+1:end);
        out.x = decomposed_recover(x, model_new.info);
    end
    
    out.conv = convhull(out.a, out.b);
end

function model_out = basis_change(x, model)
    K = model.K;
    
    basis = cell(length(K.s), 1);
    
    Count = K.f + K.l + sum(K.q);
    
    %delta = 1e-5;
    
    %NO_C = (nargin < 3) 
    
    model_out = model;
    
    for i = 1:length(K.s)
        Ksi = K.s(i);
        ind = Count + (1:Ksi^2);
        
        %Find new basis
        X = reshape(x(ind), Ksi, Ksi);        
                
        L = chol(X);
        basis{i} = L;
        
        %perform basis change
        C = reshape(model.c(ind), Ksi, Ksi);
        C_new = L'*C*L;
        model_out.c(ind) = reshape(C_new, [], 1);
        
        At_curr = model.At(ind, :);
        for j = 1:size(model.At, 2)
            At_j = At_curr(:, j);
            At_j_mat = reshape(At_j, Ksi, Ksi);
            
            At_j_mat_new = L'*At_j_mat*L;
            At_curr(:, j) = reshape(At_j_mat_new, [], 1);
        end
        
        model_out.At(ind, :) = At_curr;
        
        Count = Count + Ksi^2;
    end
    
    
    if isfield(model, 'basis')
        for i = 1:length(model.basis)
            basis{i} = model.basis{i}*basis{i};
        end
    else
        model_out.basis = basis;
    end
end