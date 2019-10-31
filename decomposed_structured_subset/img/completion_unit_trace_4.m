N = 200;
th = linspace(0,2*pi, N);


a = sdpvar(1,1);
b = sdpvar(1,1);
c = sdpvar(1,1);
X = sdpvar(4,4, 'symmetric');

C = [X(1,1) == 1; X(2,2) == 1; X(3,3) == 1; X(4, 4) == 1; X>=0];

%C = [X == eye(3)];
C = [C; X(1,2) == a; X(2,3) == -2*a; X(2, 4) == a + b; X(3, 4) == b/2];

obj = 100*a + 200*b;

%put this in primal form tomorrow

%,'sparsecolo','sparsecolo.domain',1,'sparsecolo.range',0,'sparsecolo.EQorLMI',1,'sparsecolo.SDPsolver','sedumi'


opts = sdpsettings('solver', 'sedumi', 'verbose', 0);


[model, recovery_model] = export(C, obj, opts);
model.pars.fid = 0;

out    = draw_feasibility(model, 'psd', th);
out_dd = draw_feasibility(model, 'dd', th);

% sola = optimize(C, obj_a, opts);
% va_a = value(a);
% va_b = value(b);
% solb = optimize(C, obj_b, opts);
% vb_a = value(a);
% vb_b = value(b);
% [x,y,info] = sedumi(model.A, model.b, model.C, model.K, model.pars);
%value([a,b])

figure(30)
clf
hold on
plot(out.a, out.b, 'k', 'linewidth', 2)
plot(out_dd.a, out_dd.b, 'k', 'linewidth', 2)
hold off
function out = draw_feasibility(model, cone, th)
    N = length(th);
    out.a  = zeros(N, 1);
    out.b  = zeros(N, 1);
    out.x  = cell(N, 1);
    
    %c = %model.C;
    %c1 = c(1);
    %c2 = c(2);
    
    [model_new.A,model_new.b,model_new.C,model_new.K, ~] = ...
        decomposed_subset(model.A, model.b, model.C, model.K, cone);
    model_new.pars = model.pars;
    
    
    %yalmip poses this as a dual optimization problem for some reason
    
    for i = 1:N
        theta = th(i);                
        
        b = zeros(size(model_new.b));
        
        b(1) = -cos(theta);
        b(2) = -sin(theta);                
        %c = c1*cos(theta) + c2*sin(theta);
        
        
        %there's a bug here, the chordal decomposition did not add new
        %equality constraints. How to add these in?
        [x,y,info] = sedumi(model_new.A, b, model_new.C, model_new.K, model_new.pars);
        %K = model.K;
        out.a(i) = y(1);
        out.b(i) = y(2);
        %out.x{i} = x;
        %out.a(i) = c1'*x;
        %out.b(i) = c2'*x;
        %out.Q{i} = x(K.f+1:end);
    end
    
    %out.conv = convhull(out.a, out.b);
end