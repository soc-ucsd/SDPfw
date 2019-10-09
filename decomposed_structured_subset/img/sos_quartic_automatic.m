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
        OUT{i, j} = draw_feasibility(model_agler, cone_curr, th);
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
        
        title(strcat('$K_1=', upper(cones{i}), ',$  $K_2=', upper(cones{j}), '$'), 'Interpreter', 'latex', 'FontSize', 18)
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
    
    if iscell(cone)
        [model_new.At, model_new.b, model_new.c, model_new.K, ~] = ...
        decomposed_subset(model.At',model.b,model.c,model.K,cone);
        model_new.pars = model.pars;
    else
        model_new = model;
    end
    
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
    end
    
    out.conv = convhull(out.a, out.b);
end