% load('quartic_data_3m.mat');
load('quartic_agler.mat')

model_agler.At = At;
model_agler.b = b;
model_agler.c = c;
model_agler.K = K;
model_agler.pars.fid = 0;

% %% Set up the chordal decomposition
% pars.fid = 0;
% 
% model.At = At;
% model.b  = b;
% model.c  = c;
% model.K  = K;
% model.pars = pars;
% 
% cons = [1, 2, 3, 4, 5:12 15];
% 
% K.l = 0;
% K.q = [];
% 
% %I found the bug. Wrong cliques being found, didn't index out K.f
% [cl, Ech, Jch, usedvars, s] = chordalDecomposition(model.At(K.f+1:end, cons), model.c(K.f+1:end), K);
% 
% cliques = {};
% cl_ind = 0;
% idx = (1:K.f)';
% 
% for i = 1:(cl{1,1}.NoC)
%     cli = cl{1 ,1}.Elem (cl_ind + (1:cl{1,1}.NoElem(i)));
%     cl_ind = cl_ind + cl{1,1}.NoElem(i);
%     cliques{i} =  cli;
%     idxi = reshape(cli +  K.s*(cli-1)' + K.f, [],1 );
%     idx = [idx; idxi];
%     if i ==1
%         idx0 = idxi;
%     end
% end
% 
% model_agler.K.f = K.f;
% model_agler.K.l = 0;
% model_agler.K.q = [];
% model_agler.K.s = cellfun(@length, cliques);
% model_agler.pars = pars;
%  
% numvar = K.f + sum(model_agler.K.s .^2);
% 
% 
% model_agler.At = At(idx, cons );
% model_agler.b = b(cons);
% model_agler.c  = c(idx);


%% Start the change of basis process
%Initial feasible point
model_feas_0 = model_agler;
model_feas_0.At = [model_feas_0.At, [eye(2); zeros(size(model_feas_0.At,1)-2, 2)]];
model_feas_0.b = [model_feas_0.b; 0; 0];

[x0, y0, info0] = sedumi(model_feas_0.At, model_feas_0.b, model_feas_0.c, model_feas_0.K, model_feas_0.pars);


%basis change
N = 200;
iter_change = 4;
eps_change = 0.05;

%model_change = cell(iter_reweight, 1);

%Travel in the direction of the corner:
theta_corner = 0.762159001611428 + pi;

model_change = basis_change(x0, model_agler);
%[x0c, y0c, info0c] = sedumi(model_change.At, model_change.b, model_change.c, model_change.K);

%Now do the decompositions

%same decompositions
%model_new.pars = pars;
%cone = 'sdd';


%% Draw feasibility and run
cones = {'dd', 'sdd', 'psd'};
LC = length(cones);
%LC = 1;
%OUT = cell(LC, LC);

%set up reweighting cycle
OUT = cell(iter_change, LC, LC);
MDL= cell(iter_change+1, LC, LC);

for i = 1:LC
    for j = 1:LC
        MDL{1,i,j} = model_change;
    end
end


%N = 200;
th = linspace(0, 2*pi, N);

OUT_0 = draw_feasibility(model_change, 0, th, theta_corner);

for w = 1:iter_change
    for i = 1:LC
        for j = 1:LC
            cone_curr = {cones{i}, cones{j}};
            model_curr = MDL{w,i,j};
            
            if (i == LC)&& (j== LC) && (w > 1)
                OUT{w, i, j} = OUT{1, i, j};
            else                

                %draw feasibility set
                OUT{w, i, j} = draw_feasibility(model_curr, cone_curr, th, theta_corner);



                %update basis for next iteration
                x_star = OUT{w,i,j}.x_star;    

                %SDPs will land on the boundary, move a little bit inside
                %before doing the basis change.
                x_change = x_star*(1-eps_change) + x0*eps_change;

                [MDL{w+1, i, j}, x_fake] = basis_change(x_change, model_curr);
            end
        end
    end
end


figure(1)
clf
C = linspecer(6);

for i = 1:LC
    for j= 1:LC
        subplot(LC, LC, i + LC*(j-1))
        hold on
        
        plot(OUT_0.a, OUT_0.b, 'color', 'k', 'linewidth', 2)
        
        for w = 1:iter_change
            OUT_curr = OUT{w, i, j};
            plot(OUT_curr.a, OUT_curr.b, 'color', C(w, :), 'linewidth', 2)

            plot(0, 0, 'xk', 'markersize', 12)
            quiver(0.5, -1.5, -cos(theta_corner)*0.8, -sin(theta_corner)*0.8, 'k', 'MaxHeadSize', 2, 'LineWidth', 2)
            plot(OUT_curr.a_star, OUT_curr.b_star, '*', 'color', C(w, :), 'linewidth', 2)
        end
        title(strcat('$K_4=', upper(cones{i}), ',$  $K_5=', upper(cones{j}), '$'), 'Interpreter', 'latex', 'FontSize', 18)
        xlabel('a')
        ylabel('b')
        hold off
        axis square

    end
end


function out = draw_feasibility(model, cone, th, theta_star)
    %Draw the feasibility plot
    %model: Original sdp
    %cone:  decomposed structured subset of interest in the problem
    %th:    list of angles
    %th_star: special cost angle, used for next iteration of change of
    %basis


    if nargin < 4
        th_star = 0;
    end
    
    N = length(th);
    out.a  = zeros(N+1, 1);
    out.b  = zeros(N+1, 1);
    out.Q  = cell(N+1, 1);
    
    [model_new.At, model_new.b, model_new.c, model_new.K, model_new.info] = ...
        decomposed_subset(model.At',model.b,model.c,model.K,cone);
        model_new.pars = model.pars;
    
    out.model_dec = model_new;

    

    
    for i = 1:N
        theta = th(i);
        
        [out.a(i), out.b(i), ~] = run_model(model_new, theta);
    end
    
    
    [out.a_star, out.b_star, out.x_star] = run_model(model_new, theta_star);
    
    %periodic theta
    out.a(end) = out.a(1);
    out.b(end) = out.b(1);
    
    %out.conv = convhull(out.a, out.b);
end

function [a, b,  x] = run_model(model_new, theta)
    %runs the feasibility problem where a and b are in direction theta
        c = model_new.c;
        c(1) = cos(theta);
        c(2) = sin(theta);
        
        [x,y,info] = sedumi(model_new.At, model_new.b, c, model_new.K, model_new.pars);
        a = x(1);
        b = x(2);       
        x = decomposed_recover(x, model_new.info);    
end

