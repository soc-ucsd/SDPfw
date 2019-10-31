N = 30;
th = linspace(0,2*pi, N);


a = sdpvar(1,1);
b = sdpvar(1,1);
c = sdpvar(1,1);
X = sdpvar(3,3, 'symmetric');

C = [X(1,1) == 1; X(2,2) == 1; X(3,3) == 1; X>=0];

%C = [X == eye(3)];
C = [C; X(1,2) == 0.25-a; X(2,3) == 0.5+a-b];

obj_a = a;
obj_b = b;

%,'sparsecolo','sparsecolo.domain',1,'sparsecolo.range',0,'sparsecolo.EQorLMI',1,'sparsecolo.SDPsolver','sedumi'


opts = sdpsettings('solver', 'sedumi', 'verbose', 0);


[model_a, recovery_model] = export(C, obj_a, opts);
[model_b, recovery_model] = export(C, obj_b, opts);
ba = model_a.b;
bb = model_b.b;


%model.pars.fid = 1;
out = draw_feasibility(model_a, th);

% sola = optimize(C, obj_a, opts);
% va_a = value(a);
% va_b = value(b);
% solb = optimize(C, obj_b, opts);
% vb_a = value(a);
% vb_b = value(b);
% [x,y,info] = sedumi(model.A, model.b, model.C, model.K, model.pars);
%value([a,b])

figure(30)
plot(out.a, out.b, 'k', 'linewidth', 2)

function out = draw_feasibility(model, th)
    N = length(th);
    out.a  = zeros(N, 1);
    out.b  = zeros(N, 1);
    out.x  = cell(N, 1);
    
    %c = %model.C;
    %c1 = c(1);
    %c2 = c(2);
    
    
    %yalmip poses this as a dual optimization problem for some reason
    
    for i = 1:N
        theta = th(i);                
        
        b = zeros(size(model.b));
        
        b(1) = -cos(theta);
        b(2) = -sin(theta);                
        %c = c1*cos(theta) + c2*sin(theta);
        
        
        %there's a bug here, the chordal decomposition did not add new
        %equality constraints. How to add these in?
        [x,y,info] = sedumi(model.A, b, model.C, model.K, model.pars);
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