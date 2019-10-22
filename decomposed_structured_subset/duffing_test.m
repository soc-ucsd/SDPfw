N = 5;
b = (rand(N) + 0.5)/N;
a = (rand(N,1)+0.5)/N/N;
y = sdpvar(N, 1);

%Va = a'*(y.^2 / 2 + y.^4/4);
%Vb = sum(sum(b.*((y*ones(1, N) - ones(N, 1)*y').^4)));
V = sum(y.^4);
%V = Va + Vb;

lam = sdpvar(N, 1);
%g = sdpvar(1, 1);
g = 2;

mult_term = y.^2;
%mult_term = y.^2 .* (g-y.^2);
p = V - lam'*mult_term;



F = sos(p);
%obj = -g;
obj = 0;
opts_csp = sdpsettings('solver','moment');
opts_csp.sos.csp = 0;
[model, recoverymodel] = export(F , obj, opts_csp);
K0s = model.K.s;

opts_csp.sos.csp = 1;
[model_csp, recoverymodel_csp] = export(F , obj, opts_csp);
Kcsps = model_csp.K.s;