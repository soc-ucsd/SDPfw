%feasibility quartic example through YALMIP

x = sdpvar(3, 1);
x1 = x(1);
x2 = x(2);
x3 = x(3);

a = sdpvar(1, 1);
b = sdpvar(1, 1);

f = x1^4 + (1+a)*x2^4 + (2+b)*x3^4 + (2-a-b)*x1^2*x2^2 / 2 + 2*b*x1*x2^3 ...
    + (a-b)*x2^2*x3^2 + 3*a*x1*x3^3 + 2*x1^2*x2*x3;

F = sos(f);
obj =  2*a + 3*b;

opts_csp = sdpsettings('solver','SEDUMI');

opts_csp.sos.csp = 1;
[model_csp, recoverymodel_csp] = export(F , obj, opts_csp);
