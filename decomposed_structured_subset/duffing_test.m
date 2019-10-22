N = 20;
b = (rand(N) + 0.5)/N;
a = (rand(N,1)+0.5)/N/N;
y = sdpvar(N, 1);

Va = a'*(y.^2 / 2 + y.^4/4);
%Vb = sum(sum(b.*((y*ones(1, N) - ones(N, 1)*y').^4)));
Vb = sum(sum(b.*((y*ones(1, N) - ones(N, 1)*y').^4)))/4;
%V = sum(y.^2 - y.^4);
%V = Va + Vb;
V = Va + Vb;

lam = sdpvar(N, 1);
g = sdpvar(1, 1);
g = 100;
%g = 2;

%mult_term = y.^2;
%mult_term = y.^2 .* (g-y.^2);
mult_term = ones(size(y));
%p = V - lam'*mult_term;
p = V  - sum(lam.* (g-y.^2) .* (y.^2));
%dVdy = jacobian(V, y);
%[sL,cL] = polynomial(x,2);

%p = (x'*x)*(V - rho) - sL*Vdot;

%F = [sos(p); sos(sL)];
F = [sos(p), lam>=0];
obj = -g;
%obj = 0;
opts_csp = sdpsettings('solver','sedumi');

opts_csp.sos.csp = 1;
[model_csp, recoverymodel_csp] = export(F , obj, opts_csp);
Kcsps = model_csp.K.s;

sol = solvesos(F,obj,opts_csp,[g; lam]);
rho_g = value(g)
opts_csp.sos.csp = 0;
[model, recoverymodel] = export(F , obj, opts_csp);
K0s = model.K.s;

%opts_csp.sos.csp = 1;
%[model_csp, recoverymodel_csp] = export(F , obj, opts_csp);
%Kcsps = model_csp.K.s;

%sol = optimize(F, obj, opts_csp) 