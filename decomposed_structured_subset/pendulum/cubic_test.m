x = sdpvar(1, 1);
rho = sdpvar(1,1);

f = x^3 - x;
V = x^2;
dVdx = jacobian(V, x);
Vdot = dVdx' * f;

%L = 1;
%degL = 2;
degL = 4;
[sL,cL] = polynomial(x,2);


p = (x'*x)*(V - rho) - sL*Vdot;

cons = [sos(p); sos(sL)]; 
obj  = -rho;

opts = sdpsettings('solver', 'sedumi');

%sol = optimize(cons, obj, opts)
sol = solvesos(cons,-rho,opts,[rho; cL]);
rho_V = value(rho)