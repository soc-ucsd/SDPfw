N = 200;
th = linspace(0,2*pi, N);

%PSD matrix  X(4, 4)


%i = [1; 2; 3; 4];
% i = [1;2;3;4] ;
% j = [3;8;13;18 ];
% v = [1; 1; 1; 1];


%diagonal [1; 2-b; 5; 2]
i = [1;2;3;4] ;
j = [3; 8; 13;18 ];
v = [1; 1; 1; 1];
b = [1; 2; 5; 2];
% 
% i = [i; 2; 2];
% j = [j; 8; 2];
% v = [v; 1; -1/2];


%X(1, 2) = 1+a
i =  [i; 5; 5; 5];
j =  [j; 1; 4; 7];
v =  [v; -2; 1; 1];
b = [b; 1];

%X(2,4) = a + b
i =  [i; 6; 6; 6 ; 6];
j =  [j; 1; 2;  10; 16];
v =  [v; -2; -2; 1; 1];
b = [b; 0];


%X(2,3) =  -2a
i =  [i;  7; 7; 7];
j =  [j; 1;  9; 12];
v =  [v; 4; 1;  1];
b =  [b; 0];

%X(3,4) =  b/2
i =  [i;  8; 8; 8];
j =  [j; 2;  14; 17];
v =  [v; -1; 1;  1];
b =  [b; 0];

model.A = sparse(i, j, v);
model.b = b;
model.C = sparse(18, 1);
model.K.f = 2;
model.K.s = 4;
model.K.l = 0;
model.K.q = [];
model.pars.fid = 0;

%% manually do the split because sparsecolo doesn't want to.

%start with clique overlap


is = [1; 1];
js = [3; 15];
vs = [1; -1];
bs = 0;
%diagonal [1; 1; 5; 2]
is = [is; 2;  3; 4; 5];
js = [js; 12; 3; 7; 11];
vs = [vs; 1;  1; 1; 1];
bs = [bs; 1;  2; 5; 2];
%join is 1-b/2
% is = [is; 3; 3];
% js = [js; 3; 2];
% vs = [vs; 1; -1/2];

%bs = [bs; 1; 1; 1; 1];

%X(1, 2) = 1+a
is =  [is; 6; 6; 6];
js =  [js; 1; 13; 14];
vs =  [vs; -2; 1; 1];
bs =  [bs; 1];

%X(2,4) = a + b
is =  [is; 7;7;7;7];
js =  [js; 1; 2;  5; 9];
vs =  [vs; -2; -2; 1; 1];
bs =  [bs; 0];
% 

%X(2,3) =  -2a
is =  [is;  8;8;8];
js =  [js; 1;  4; 6];
vs =  [vs; 4; 1;  1];
bs =  [bs; 0];
% 
%X(3,4) =  b/2
is =  [is;  9;9;9];
js =  [js; 2;  8; 10];
vs =  [vs; -1; 1;  1];
bs =  [bs; 0];

model_split.A = sparse(is, js, vs);
%model_split.b = sparse(9, 1) ;
model_split.b = bs;
model_split.C = sparse(15, 1);

model_split.K.f = 2;
model_split.K.s = [3 2];                                                                                  
model_split.K.l = 0;
model_split.K.q = [];
model_split.pars.fid = 0;
cliques = {{2,3,4}, {1,2}}; 

% Cs = model_split.C;
% 
% Cs(1) = -1;
% Cs(2) = -2;
% [xs, ys, infos] = sedumi(model_split.A, model_split.b, Cs, model_split.K, model_split.pars);
% 
% as = xs(1);
% bs = xs(2);
% x1s = reshape(xs(3:11), 3, 3);
% x2s = reshape(xs(12:15), 2, 2);
% 
% C = model.C;
% C(1) = -1;
% C(2) = -2;
% 
% 
% [x, y, info] = sedumi(model.A, model.b, C, model.K, model.pars);
% a0 = x(1);
% b0 = x(2);
% x0 = reshape(x(3:18), 4, 4);
% 
% da = a0 - as;
% db = b0 - bs;



out    = draw_feasibility(model, 'psd', th);
out_dd = draw_feasibility(model, 'dd', th);
out_sdd = draw_feasibility(model, 'sdd', th);

out_c = draw_feasibility(model_split,    'psd', th);
out_cdd = draw_feasibility(model_split,  'dd', th);
out_csdd = draw_feasibility(model_split, 'sdd', th);
out_cdd_sdd = draw_feasibility(model_split, {'dd', 'sdd'}, th);
out_csdd_dd = draw_feasibility(model_split, {'sdd', 'dd'}, th);



C = linspecer(6);

%homogenous cones
    figure(30)
    clf
    subplot(2, 2, 1)
    hold on
    % plot(out.a, out.b, 'k', 'linewidth', 2)
    % plot(out_dd.a, out_dd.b, 'linewidth', 2, 'color', C(1, :) )
    % plot(out_cdd.a, out_cdd.b, 'linewidth', 2, 'color', C(2, :))
    plot_region(out, [0, 0, 0])
    plot_region(out_cdd, C(2, :))
    plot_region(out_dd, C(1, :))
    hold off
    xlabel('$a$', 'interpreter', 'latex', 'Fontsize', 16)
    ylabel('$b$', 'interpreter', 'latex', 'Fontsize', 16)
    axis square
    title('$\mathcal{D}\mathcal{D}^4 \subset \mathcal{D}\mathcal{D}^4(E, ?) \subset S_+^4$', 'interpreter', 'latex', 'Fontsize', 18)
    l1 = {'$S_+^4$', ...
        '$\mathcal{D}\mathcal{D}^4(E, 0)$', ...
        '$\mathcal{D}\mathcal{D}^4$'  };    
    legend(l1, 'interpreter', 'latex', 'location', 'SouthWest', 'fontsize', 12);

    subplot(2, 2, 2)
    hold on
    % plot(out.a, out.b, 'k', 'linewidth', 2)
    % plot(out_cdd.a, out_cdd.b, 'linewidth', 2, 'color', C(2, :) )
    % plot(out_csdd.a, out_csdd.b, 'linewidth', 2, 'color', C(6, :))
    plot_region(out, [0, 0, 0])
    plot_region(out_csdd, C(6, :))
    plot_region(out_cdd, C(2, :))
    hold off
    xlabel('$a$', 'interpreter', 'latex', 'Fontsize', 16)
    ylabel('$b$', 'interpreter', 'latex', 'Fontsize', 16)
    axis square
    title('$\mathcal{D}\mathcal{D}^4(E, 0) \subset \mathcal{S}\mathcal{D}\mathcal{D}^4(E, ?) \subset S_+^4$', 'interpreter', 'latex', 'Fontsize', 18)

    l1 = {'$S_+^4$', ...
        '$\mathcal{S}\mathcal{D}\mathcal{D}^4(E, 0)$', ...
        '$\mathcal{D}\mathcal{D}^4(E, 0)$'};    
    legend(l1, 'interpreter', 'latex', 'location', 'SouthWest', 'fontsize', 12);


    subplot(2, 2, 3)
    hold on
    % plot(out.a, out.b, 'k', 'linewidth', 2)
    % plot(out_dd.a, out_dd.b, 'linewidth', 2, 'color', C(1, :))
    % plot(out_sdd.a, out_sdd.b, 'linewidth', 2, 'color', C(5, :))
    plot_region(out, [0, 0, 0])
    plot_region(out_sdd, C(5, :))
    plot_region(out_dd, C(1, :))
    hold off
    xlabel('$a$', 'interpreter', 'latex', 'Fontsize', 16)
    ylabel('$b$', 'interpreter', 'latex', 'Fontsize', 16)
    axis square
    title('$\mathcal{D}\mathcal{D}^4 \subset \mathcal{S}\mathcal{D}\mathcal{D}^4 \subset S_+^4$', 'interpreter', 'latex', 'Fontsize', 18)
    l1 = {'$S_+^4$', ...         
         '$\mathcal{S}\mathcal{D}\mathcal{D}^4$',...
         '$\mathcal{D}\mathcal{D}^4$',};    
    legend(l1, 'interpreter', 'latex', 'location', 'SouthWest', 'fontsize', 12);



    subplot(2, 2, 4)
    hold on
    %plot(out.a, out.b, 'k', 'linewidth', 2)
    %plot(out_sdd.a, out_sdd.b, 'linewidth', 2, 'color', C(5, :))
    %plot(out_csdd.a, out_csdd.b, 'linewidth', 2, 'color', C(6, :))
    plot_region(out, [0, 0, 0])
    plot_region(out_csdd, C(6, :))
    plot_region(out_sdd, C(5, :))

    hold off
    xlabel('$a$', 'interpreter', 'latex', 'Fontsize', 16)
    ylabel('$b$', 'interpreter', 'latex', 'Fontsize', 16)
    axis square
    title('$\mathcal{S}\mathcal{D}\mathcal{D}^4 \subset \mathcal{S}\mathcal{D}\mathcal{D}^4(E, ?) \subset S_+^4$', 'interpreter', 'latex', 'Fontsize', 18)
    l1 = {'$S_+^4$', ...
         '$\mathcal{S}\mathcal{D}\mathcal{D}^4$', ...
         '$\mathcal{S}\mathcal{D}\mathcal{D}^4(E, ?)$'};
    legend(l1, 'interpreter', 'latex', 'location', 'SouthWest', 'fontsize', 12);

%mixed cones
figure(31)
clf
subplot(1, 2, 1)
hold on
%plot(out.a, out.b, 'k', 'linewidth', 2)
%plot(out_cdd.a, out_cdd.b, 'linewidth', 2, 'color', C(2, :) )
%plot(out_csdd_dd.a, out_csdd_dd.b, 'linewidth', 2, 'color', C(4, :))
%plot(out_csdd.a, out_csdd.b, 'linewidth', 2, 'color', C(6, :))
plot_region(out_csdd, C(6, :))
plot_region(out_csdd_dd, C(4, :))
plot_region(out_cdd, C(2, :))
hold off

xlabel('$a$', 'interpreter', 'latex', 'Fontsize', 16)
ylabel('$b$', 'interpreter', 'latex', 'Fontsize', 16)
axis square
title('$(\mathcal{D}\mathcal{D}^3, \mathcal{D}\mathcal{D}^2) \subset (\mathcal{S}\mathcal{D}\mathcal{D}^3, \mathcal{D}\mathcal{D}^2) \subset (\mathcal{S}\mathcal{D}\mathcal{D}^3, \mathcal{S}\mathcal{D}\mathcal{D}^2)$', 'interpreter', 'latex', 'Fontsize', 18)
l1 = {'$(\mathcal{S}\mathcal{D}\mathcal{D}^3, \mathcal{S}\mathcal{D}\mathcal{D}^2)$', ...
    '$(\mathcal{S}\mathcal{D}\mathcal{D}^3, \mathcal{D}\mathcal{D}^2)$', ...
    '$(\mathcal{D}\mathcal{D}^3, \mathcal{D}\mathcal{D}^2)$'};
legend(l1, 'interpreter', 'latex', 'location', 'SouthWest', 'fontsize', 14);

subplot(1, 2, 2)
hold on
%plot(out.a, out.b, 'k', 'linewidth', 2)
% plot(out_cdd.a, out_cdd.b, 'linewidth', 2, 'color', C(2, :) )
% plot(out_cdd_sdd.a, out_cdd_sdd.b, 'linewidth', 2, 'color', C(5, :))
% plot(out_csdd.a, out_csdd.b, 'linewidth', 2, 'color', C(6, :))

plot_region(out_csdd, C(6, :))
plot_region(out_cdd_sdd, C(3, :))
plot_region(out_cdd, C(2, :))


hold off
xlabel('$a$', 'interpreter', 'latex', 'Fontsize', 16)
ylabel('$b$', 'interpreter', 'latex', 'Fontsize', 16)
axis square
title('$(\mathcal{D}\mathcal{D}^3, \mathcal{D}\mathcal{D}^2) \subset (\mathcal{D}\mathcal{D}^3, \mathcal{S}\mathcal{D}\mathcal{D}^2) \subset (\mathcal{S}\mathcal{D}\mathcal{D}^3, \mathcal{S}\mathcal{D}\mathcal{D}^2) $', 'interpreter', 'latex', 'Fontsize', 18)
l1 = {'$(\mathcal{S}\mathcal{D}\mathcal{D}^3, \mathcal{S}\mathcal{D}\mathcal{D}^2)$', ...
    '$(\mathcal{D}\mathcal{D}^3, \mathcal{S}\mathcal{D}\mathcal{D}^2)$', ...
    '$(\mathcal{D}\mathcal{D}^3, \mathcal{D}\mathcal{D}^2)$'};
legend(l1, 'interpreter', 'latex', 'location', 'SouthWest', 'fontsize', 14);
function out = draw_feasibility(model, cone, th)
    N = length(th);
    out.a  = zeros(N, 1);
    out.b  = zeros(N, 1);
    out.x  = cell(N, 1);
    
    %c = %model.C;
    %c1 = c(1);
    %c2 = c(2);
    
    [model_new.A,model_new.b,model_new.C,model_new.K, model_new.info] = ...
        decomposed_subset(model.A, model.b, model.C, model.K, cone);
    model_new.pars = model.pars;
    
    
    %yalmip poses this as a dual optimization problem for some reason
    
    for i = 1:N
        theta = th(i);                
        
        C = zeros(size(model_new.C));
        
        C(1) = -cos(theta);
        C(2) = -sin(theta);                
        %c = c1*cos(theta) + c2*sin(theta);
        
        
        %there's a bug here, the chordal decomposition did not add new
        %equality constraints. How to add these in?
        [x,y,info] = sedumi(model_new.A, model_new.b, C, model_new.K, model_new.pars);
        %K = model.K;
        x = decomposed_recover(x, model_new.info);  
        out.a(i) = x(1);
        out.b(i) = x(2);
        
        if length(model.K.s) > 1
            out.x{i} = {reshape(x(3:11), 3, 3), reshape(x(12:15), 2, 2)};
        else
            out.x{i} = reshape(x(3:end), 4, 4);
        end
    end
    
    %out.conv = convhull(out.a, out.b);
end

